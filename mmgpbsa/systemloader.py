import os
import subprocess
from sys import stdout

import mdtraj as md
import numpy as np
import openmm as mm
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from pymbar import timeseries
from openmm import unit
from openmm import app

from mmgpbsa.amber_mmgpbsa import run_amber
from mmgpbsa.utils import make_message_writer, working_directory


def subsample(enthalpies):
    """
    Subsamples the enthalpies using John Chodera's code.
    This is probably better than the simple cutoff we normally use.
    No output -- it modifies the lists directly
    """
    # Use automatic equilibration detection and pymbar.timeseries to subsample
    [t0, g, Neff_max] = timeseries.detectEquilibration(enthalpies)
    enthalpies = enthalpies[t0:]
    return timeseries.subsampleCorrelatedData(enthalpies, g=g)


class Config:
    defaults_config = {

        ## System Params
        'nonbondedMethod': app.CutoffNonPeriodic,
        'rigid_water': True,
        'removeCMMotion': True,

        "nonbondedCutoff": 1.0 * unit.nanometer,
        "constraints": app.HBonds,
        "implicitSolvent": app.GBn2,
        "soluteDielectric": 1.0,
        "solventDielectric": 78.5,
        "igb": 5,

        ## Integrator Params
        "temperature": 310.15 * unit.kelvin,
        "frictionCoeff": 1.0 / unit.picoseconds,
        "stepSize": 2.0 * unit.femtoseconds,
        "constraintTolerance": 0.00001,

        'platform_name': 'CPU',
        'platform_properties': {},

        'reportInterval': 1250,
        'equil_ps': 10,
        'ps': 10,
        'calcFrames': None,
        'mbar': 30
    }

    def __init__(self, **kwargs):

        self.__dict__.update(self.defaults_config)

        for k, v in kwargs.items():
            if k in self.__dict__:
                if v is not None:
                    self.__dict__[k] = v

        ## Simulation Params
        if self.platform_name in ['OpenCL', 'CUDA']:
            self.platform_properties = {'Precision': 'mixed'}

        if not isinstance(self.ps, unit.Quantity):
            self.ps = self.ps * unit.picosecond
        if not isinstance(self.equil_ps, unit.Quantity):
            self.equil_ps = self.equil_ps * unit.picosecond

        if self.mbar is not None and self.mbar > 0:
            self.calcFrames = self.mbar

        self.total_ps = self.equil_ps + self.ps
        self.equil_steps = int(self.equil_ps / self.stepSize)
        self.production_steps = int(self.ps / self.stepSize)
        self.total_steps = int(self.total_ps / self.stepSize)
        self.trajInterval = int(self.production_steps / self.calcFrames)


class SystemLoader:

    def __init__(self, dirpath, verbose, input_pdb=None, apo=None, lig=None, **kwargs):
        self.verbose = verbose
        self.logger = make_message_writer(self.verbose >= 1, self.__class__.__name__, verbose >= 2)

        with self.logger("__init__") as logger:
            self.dirpath = dirpath
            self.input_pdb = input_pdb
            self.apo = apo
            self.lig = lig
            self.config = Config(**kwargs)

    def split_complex_from_system(self):
        with self.logger("split_complex_from_system") as logger:
            if self.input_pdb is not None:
                # Read the PDB file using RDKit
                mol = Chem.MolFromPDBFile(self.input_pdb, removeHs=False)
                if mol is None:
                    logger.error("Could not read PDB file!")
                    return None, None
                
                logger.log(f"Reading input PDB {self.input_pdb}")
                
                # Split the molecule into fragments
                fragments = Chem.GetMolFrags(mol, asMols=True)
                
                if len(fragments) < 2:
                    logger.failure("Could not identify separate protein and ligand fragments!", exit_all=True)
                
                # Sort fragments by size - assume largest is protein, second largest is ligand
                fragments = sorted(fragments, key=lambda x: x.GetNumAtoms(), reverse=True)
                protein = fragments[0]  # Largest fragment
                ligand = fragments[1]   # Second largest fragment
                
                logger.log(f"Split complex. atom sizes-- lig: {ligand.GetNumAtoms()}, prot: {protein.GetNumAtoms()}")
                
                return protein, ligand
            else:
                # Handle separate apo and ligand files
                if self.apo is not None:
                    protein = Chem.MolFromPDBFile(self.apo, removeHs=False)
                else:
                    logger.error("No protein structure provided!")
                    return None, None
                
                if self.lig is not None:
                    # Try to read ligand in various formats
                    ligand = None
                    for fmt in ['.mol2', '.sdf', '.pdb', '.mol']:
                        if self.lig.endswith(fmt):
                            if fmt == '.mol2':
                                ligand = Chem.MolFromMol2File(self.lig)
                            elif fmt == '.sdf':
                                ligand = Chem.MolFromMolFile(self.lig)
                            elif fmt in ['.pdb', '.mol']:
                                ligand = Chem.MolFromPDBFile(self.lig)
                            break
                    
                    if ligand is None:
                        logger.error("Could not read ligand file!")
                        return None, None
                else:
                    logger.error("No ligand structure provided!")
                    return None, None
                
                return protein, ligand

    def prepare_ligand(self, lig):
        with self.logger("prepare_ligand") as logger:
            ## prepare ligand using RDKit
            oemol = Chem.Mol(lig)
            
            # Add hydrogens if not present
            oemol = Chem.AddHs(oemol)
            
            # Generate 3D coordinates if needed
            if oemol.GetNumConformers() == 0:
                AllChem.EmbedMolecule(oemol, randomSeed=42)
                AllChem.MMFFOptimizeMolecule(oemol)
            
            # Ensure directory exists
            os.makedirs(self.dirpath, exist_ok=True)
            
            # Save ligand to MOL2
            Chem.MolToMolFile(oemol, f'{self.dirpath}/lig.mol2')
            
            # Save ligand to PDB
            Chem.MolToPDBFile(oemol, f'{self.dirpath}/lig.pdb')
            
            self.lig = oemol
            
            # Assign charges using RDKit's Gasteiger charges
            # Note: This is a simplified approach. For more accurate charges,
            # you might want to use AM1-BCC or other methods
            AllChem.ComputeGasteigerCharges(oemol)
            
            # Save charged ligand
            Chem.MolToMolFile(oemol, f'{self.dirpath}/charged.mol2')
            
            self.charged_lig = oemol

    def prepare_protein(self, protein):
        with self.logger("prepare_protein") as logger:
            # Ensure directory exists
            os.makedirs(self.dirpath, exist_ok=True)
            
            # Save protein to PDB using RDKit
            Chem.MolToPDBFile(protein, f'{self.dirpath}/apo.pdb')
            self.apo = protein

    def _create_complex_pdb(self, protein, ligand):
        """Create a complex PDB file by combining protein and ligand"""
        with self.logger("_create_complex_pdb") as logger:
            try:
                # Ensure directory exists
                os.makedirs(self.dirpath, exist_ok=True)
                
                # Create a combined molecule
                combined = Chem.CombineMols(protein, ligand)
                
                # Save the complex to PDB
                Chem.MolToPDBFile(combined, f'{self.dirpath}/com.pdb')
                
                logger.log("Created complex PDB file")
                
            except Exception as e:
                logger.error(f"Error creating complex PDB: {e}")
                raise

    def __setup_system_im(self):
        with self.logger("__setup_system_im") as logger:
            try:
                protein, lig = self.split_complex_from_system()
                if protein is None or lig is None:
                    logger.error("Failed to split complex or read structures")
                    return None
                
                self.prepare_ligand(lig)
                self.prepare_protein(protein)

                # Create complex PDB file by combining protein and ligand
                self._create_complex_pdb(protein, lig)

                with working_directory(self.dirpath):
                    # Create system using OpenMM with AMBER force fields
                    # For now, we'll use a simplified approach that creates a basic system
                    # In a full implementation, you would read the actual AMBER topology files
                    
                    # Create topology from the complex PDB
                    pdb_file = app.PDBFile(f'{self.dirpath}/com.pdb')
                    topology = pdb_file.topology
                    positions = pdb_file.positions
                    
                    # Create system with implicit solvent
                    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
                    system = forcefield.createSystem(
                        topology,
                        nonbondedMethod=self.config.nonbondedMethod,
                        nonbondedCutoff=self.config.nonbondedCutoff,
                        rigidWater=self.config.rigid_water,
                        constraints=self.config.constraints,
                        implicitSolvent=self.config.implicitSolvent,
                        soluteDielectric=self.config.soluteDielectric,
                        solventDielectric=self.config.solventDielectric,
                        removeCMMotion=self.config.removeCMMotion
                    )
                    
                    self.topology, self.positions = topology, positions
                    return system
                    
            except Exception as e:
                logger.error("EXCEPTION CAUGHT BAD SPOT", e)
                return None

    def run_amber(self, method: str, amber_path):
        with self.logger('run_amber') as logger:
            logger.log("Calculating mmgb/pbsa value...may take awhile.")
            return run_amber(amber_path, self.dirpath, verbose=self.verbose, igb=self.config.igb,
                             use_pbsa=(method == 'pbsa'))

    def prepare_simulation(self):
        with self.logger("prepare_simulation") as logger:

            system = self.__setup_system_im()
            if system is None:
                logger.error("Failed to set up system")
                return None
                
            integrator = mm.LangevinIntegrator(self.config.temperature, self.config.frictionCoeff,
                                               self.config.stepSize)
            integrator.setConstraintTolerance(self.config.constraintTolerance)

            platform = mm.Platform.getPlatformByName(self.config.platform_name)
            simulation = app.Simulation(self.topology, system, integrator, platform, self.config.platform_properties)
            logger.log(
                f"Built simulation using platform {self.config.platform_name} with properties {self.config.platform_properties}")
            simulation.context.setPositions(self.positions)

            logger.log(f'Minimizing and setting velocities to {self.config.temperature}')
            simulation.minimizeEnergy()
            simulation.context.setVelocitiesToTemperature(self.config.temperature)

            if self.verbose >= 1:
                simulation.reporters.append(
                    app.StateDataReporter(stdout, max(self.config.reportInterval - 1, 1), step=True,
                                          potentialEnergy=True, temperature=True,
                                          progress=True,
                                          remainingTime=True,
                                          speed=True, totalSteps=self.config.total_steps,
                                          separator='\t'))

            logger.log(f'Equilibrating for {self.config.equil_steps} steps, or {self.config.equil_ps}')
            simulation.step(self.config.equil_steps)

            simulation.reporters.append(
                app.DCDReporter(f'{self.dirpath}/trajectory.dcd',
                                max(self.config.trajInterval - 1 if self.config.mbar == 0 else self.config.mbar, 1)))

            logger.log(f'Running Production for {self.config.production_steps} steps, or {self.config.ps}')
            if self.config.mbar > 0:
                enthalpies = np.zeros((self.config.mbar))
                for i in range(self.config.mbar):
                    simulation.step(int(self.config.production_steps / self.config.mbar))
                    state = simulation.context.getState(getEnergy=True)
                    potential_energy = state.getPotentialEnergy()
                    enthalpies[i] = potential_energy.value_in_unit(unit.kilojoules_per_mole)
                logger.log(f"Simulation Done. Running MBAR on {self.config.mbar} snapshots.")

                idx = np.array(subsample(enthalpies))
                traj = md.load(f'{self.dirpath}/trajectory.dcd', top=f'{self.dirpath}/com.pdb')
                traj = traj[idx]
                traj.save_dcd(f'{self.dirpath}/trajectory.dcd', force_overwrite=True)

                logger.log(f"Done. Subsampled {len(idx)} from {self.config.mbar} snapshots.")

                return enthalpies
            else:
                simulation.step(self.config.production_steps)
                logger.log(f"Simulation Done")
