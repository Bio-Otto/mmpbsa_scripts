"""
XTC-based MMGBSA Calculator
Works with existing trajectory files and complex PDB structures
"""

import os
import numpy as np
import mdtraj as md
from openmm import unit
from openmm import app


class XTCMMGBSACalculator:
    """
    Calculate MMGBSA using XTC trajectory and complex PDB
    """
    
    def __init__(self, verbose=1):
        self.verbose = verbose
        
    def calculate_mmgbsa(self, complex_pdb, trajectory_file, method='gbsa', 
                        igb=5, start_frame=0, end_frame=None, stride=1):
        """
        Calculate MMGBSA using XTC trajectory and complex PDB
        
        Args:
            complex_pdb: Path to complex PDB file
            trajectory_file: Path to trajectory file (XTC, DCD, TRR)
            method: 'gbsa' or 'pbsa'
            igb: GB model (5 for GBn, 7 for GBn2, 8 for OBC)
            start_frame: Starting frame (0-based)
            end_frame: Ending frame (None for all frames)
            stride: Frame stride for analysis
            
        Returns:
            tuple: (delta_g, std_error) in kcal/mol
        """
        
        if self.verbose >= 1:
            print(f"Loading trajectory: {trajectory_file}")
            print(f"Using complex PDB: {complex_pdb}")
        
        try:
            # Load trajectory
            traj = md.load(trajectory_file, top=complex_pdb)
            
            if self.verbose >= 1:
                print(f"Loaded trajectory with {traj.n_frames} frames")
            
            # Apply frame selection
            if end_frame is None:
                end_frame = traj.n_frames
            
            traj = traj[start_frame:end_frame:stride]
            
            if self.verbose >= 1:
                print(f"Analyzing {traj.n_frames} frames (frames {start_frame}:{end_frame}:{stride})")
            
            # Split trajectory into protein and ligand
            protein_traj, ligand_traj = self._split_trajectory(traj)
            
            # Calculate energies
            complex_energies = self._calculate_trajectory_energies(traj, "complex", igb)
            protein_energies = self._calculate_trajectory_energies(protein_traj, "protein", igb)
            ligand_energies = self._calculate_trajectory_energies(ligand_traj, "ligand", igb)
            
            # Calculate MMGBSA
            if len(complex_energies) == 0 or len(protein_energies) == 0 or len(ligand_energies) == 0:
                print("Error: Could not calculate energies for all systems")
                return None, None
            
            # Calculate binding free energy
            # ΔG = G_complex - G_protein - G_ligand
            binding_energies = complex_energies - protein_energies - ligand_energies
            
            # Convert to kcal/mol
            binding_energies_kcal = binding_energies * 0.593  # Convert kJ/mol to kcal/mol
            
            # Calculate mean and standard error
            delta_g = np.mean(binding_energies_kcal)
            std_error = np.std(binding_energies_kcal) / np.sqrt(len(binding_energies_kcal))
            
            if self.verbose >= 1:
                print(f"MMGBSA calculation complete:")
                print(f"  ΔG = {delta_g:.4f} ± {std_error:.4f} kcal/mol")
                print(f"  Number of frames: {len(binding_energies_kcal)}")
            
            return delta_g, std_error
            
        except Exception as e:
            print(f"Error in MMGBSA calculation: {e}")
            return None, None
    
    def _split_trajectory(self, traj):
        """
        Split trajectory into protein and ligand parts
        Assumes ligand is the smallest molecule/residue
        """
        if self.verbose >= 2:
            print("Splitting trajectory into protein and ligand...")
        
        # Get residue information
        residues = list(traj.topology.residues)
        
        if len(residues) < 2:
            print("Error: Need at least 2 residues to split protein and ligand")
            return traj, traj
        
        # Find the smallest residue (assumed to be ligand)
        residue_sizes = [(res, len(res.atoms)) for res in residues]
        residue_sizes.sort(key=lambda x: x[1])
        
        ligand_residue = residue_sizes[0][0]  # Smallest residue
        protein_residues = [res for res in residues if res != ligand_residue]
        
        if self.verbose >= 2:
            print(f"Ligand residue: {ligand_residue.name} ({len(ligand_residue.atoms)} atoms)")
            print(f"Protein residues: {[res.name for res in protein_residues]}")
        
        # Create protein and ligand trajectories
        protein_atoms = []
        ligand_atoms = []
        
        for atom in traj.topology.atoms:
            if atom.residue == ligand_residue:
                ligand_atoms.append(atom.index)
            else:
                protein_atoms.append(atom.index)
        
        protein_traj = traj.atom_slice(protein_atoms)
        ligand_traj = traj.atom_slice(ligand_atoms)
        
        return protein_traj, ligand_traj
    
    def _calculate_trajectory_energies(self, traj, system_type, igb):
        """
        Calculate energies for all frames in a trajectory
        """
        energies = []
        
        for i, frame in enumerate(traj):
            if self.verbose >= 2 and i % 10 == 0:
                print(f"Processing {system_type} frame {i+1}/{len(traj)}")
            
            energy = self._calculate_frame_energy(frame, system_type, igb)
            if energy is not None:
                energies.append(energy)
        
        return np.array(energies)
    
    def _calculate_frame_energy(self, frame, system_type, igb):
        """
        Calculate energy for a single frame using OpenMM
        """
        try:
            # Convert frame to OpenMM system
            topology = frame.topology.to_openmm()
            positions = frame.xyz[0] * unit.nanometers
            
            # Create force field
            forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
            
            # Set up GB model
            if igb == 5:
                gb_model = app.GBn
            elif igb == 7:
                gb_model = app.GBn2
            elif igb == 8:
                gb_model = app.OBC2
            else:
                gb_model = app.GBn2  # Default
            
            # Create system
            system = forcefield.createSystem(
                topology,
                nonbondedMethod=app.CutoffNonPeriodic,
                nonbondedCutoff=1.0*unit.nanometer,
                constraints=app.HBonds,
                implicitSolvent=gb_model,
                soluteDielectric=1.0,
                solventDielectric=78.5,
                removeCMMotion=True
            )
            
            # Create context
            integrator = app.LangevinIntegrator(
                310.15*unit.kelvin,
                1.0/unit.picoseconds,
                2.0*unit.femtoseconds
            )
            
            platform = app.Platform.getPlatformByName('CPU')
            context = app.Context(system, integrator, platform)
            context.setPositions(positions)
            
            # Calculate energy
            state = context.getState(getEnergy=True)
            potential_energy = state.getPotentialEnergy()
            
            return potential_energy.value_in_unit(unit.kilojoules_per_mole)
            
        except Exception as e:
            if self.verbose >= 2:
                print(f"Error calculating energy for {system_type}: {e}")
            return None


def run_xtc_mmgbsa(complex_pdb, trajectory_file, method='gbsa', igb=5, 
                   start_frame=0, end_frame=None, stride=1, verbose=1):
    """
    Main function to run XTC-based MMGBSA calculation
    
    Args:
        complex_pdb: Path to complex PDB file
        trajectory_file: Path to trajectory file
        method: 'gbsa' or 'pbsa'
        igb: GB model for GBSA
        start_frame: Starting frame
        end_frame: Ending frame
        stride: Frame stride
        verbose: Verbosity level
        
    Returns:
        tuple: (delta_g, std_error) in kcal/mol
    """
    
    # Check if files exist
    if not os.path.exists(complex_pdb):
        print(f"Error: Complex PDB file not found: {complex_pdb}")
        return None, None
    
    if not os.path.exists(trajectory_file):
        print(f"Error: Trajectory file not found: {trajectory_file}")
        return None, None
    
    # Create calculator
    calculator = XTCMMGBSACalculator(verbose)
    
    # Run calculation
    if method == 'pbsa':
        print("Note: Using GBSA as approximation for PBSA (full PBSA not implemented)")
        return calculator.calculate_mmgbsa(
            complex_pdb, trajectory_file, method='gbsa', igb=igb,
            start_frame=start_frame, end_frame=end_frame, stride=stride
        )
    else:
        return calculator.calculate_mmgbsa(
            complex_pdb, trajectory_file, method='gbsa', igb=igb,
            start_frame=start_frame, end_frame=end_frame, stride=stride
        ) 