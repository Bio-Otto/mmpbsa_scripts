"""
OpenMM-based MMGBSA Calculator
Replaces proprietary AMBER software for MMGBSA calculations
"""

import os
import numpy as np
from simtk import unit
from simtk.openmm import app
import mdtraj as md


class OpenMMGBSACalculator:
    """
    Calculate MMGBSA using OpenMM
    """
    
    def __init__(self, work_dir, verbose=1):
        self.work_dir = work_dir
        self.verbose = verbose
        
    def calculate_mmgbsa(self, complex_traj, protein_traj, ligand_traj, 
                        complex_top, protein_top, ligand_top, 
                        temperature=310.15, igb=5):
        """
        Calculate MMGBSA using OpenMM
        
        Args:
            complex_traj: Path to complex trajectory file
            protein_traj: Path to protein trajectory file  
            ligand_traj: Path to ligand trajectory file
            complex_top: Path to complex topology file
            protein_top: Path to protein topology file
            ligand_top: Path to ligand topology file
            temperature: Temperature in Kelvin
            igb: GB model (5 for GBn, 7 for GBn2, 8 for OBC)
            
        Returns:
            tuple: (delta_g, std_error) in kcal/mol
        """
        
        if self.verbose >= 1:
            print("Calculating MMGBSA using OpenMM...")
        
        # Load trajectories
        try:
            complex_traj = md.load(complex_traj, top=complex_top)
            protein_traj = md.load(protein_traj, top=protein_top)
            ligand_traj = md.load(ligand_traj, top=ligand_top)
        except Exception as e:
            print(f"Error loading trajectories: {e}")
            return None, None
        
        # Calculate energies for each frame
        complex_energies = []
        protein_energies = []
        ligand_energies = []
        
        # Process complex
        for i, frame in enumerate(complex_traj):
            if self.verbose >= 2 and i % 10 == 0:
                print(f"Processing complex frame {i+1}/{len(complex_traj)}")
            
            energy = self._calculate_frame_energy(frame, "complex", igb)
            if energy is not None:
                complex_energies.append(energy)
        
        # Process protein
        for i, frame in enumerate(protein_traj):
            if self.verbose >= 2 and i % 10 == 0:
                print(f"Processing protein frame {i+1}/{len(protein_traj)}")
            
            energy = self._calculate_frame_energy(frame, "protein", igb)
            if energy is not None:
                protein_energies.append(energy)
        
        # Process ligand
        for i, frame in enumerate(ligand_traj):
            if self.verbose >= 2 and i % 10 == 0:
                print(f"Processing ligand frame {i+1}/{len(ligand_traj)}")
            
            energy = self._calculate_frame_energy(frame, "ligand", igb)
            if energy is not None:
                ligand_energies.append(energy)
        
        # Calculate MMGBSA
        if len(complex_energies) == 0 or len(protein_energies) == 0 or len(ligand_energies) == 0:
            print("Error: Could not calculate energies for all systems")
            return None, None
        
        # Convert to numpy arrays
        complex_energies = np.array(complex_energies)
        protein_energies = np.array(protein_energies)
        ligand_energies = np.array(ligand_energies)
        
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
    
    def _calculate_frame_energy(self, frame, system_type, igb):
        """
        Calculate energy for a single frame using OpenMM
        
        Args:
            frame: MDTraj frame
            system_type: Type of system ("complex", "protein", "ligand")
            igb: GB model
            
        Returns:
            float: Energy in kJ/mol
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
    
    def calculate_pbsa(self, complex_traj, protein_traj, ligand_traj,
                      complex_top, protein_top, ligand_top,
                      temperature=310.15, salt_concentration=0.15):
        """
        Calculate MMPBSA using OpenMM (simplified)
        
        Note: Full PBSA requires solving the Poisson-Boltzmann equation,
        which is computationally expensive. This is a simplified implementation.
        """
        
        if self.verbose >= 1:
            print("Calculating MMPBSA using OpenMM (simplified)...")
            print("Note: This is a simplified PBSA implementation")
        
        # For now, use GB calculation as approximation
        # In a full implementation, you would need to solve the PB equation
        return self.calculate_mmgbsa(
            complex_traj, protein_traj, ligand_traj,
            complex_top, protein_top, ligand_top,
            temperature, igb=5
        )


def run_openmm_mmgbsa(work_dir, verbose=1, use_pbsa=False, igb=5):
    """
    Main function to run OpenMM-based MMGBSA calculation
    
    Args:
        work_dir: Working directory containing trajectory files
        verbose: Verbosity level
        use_pbsa: Whether to use PBSA (simplified) or GBSA
        igb: GB model for GBSA
        
    Returns:
        tuple: (delta_g, std_error) in kcal/mol
    """
    
    # Check for required files
    required_files = [
        'trajectory.dcd',  # Complex trajectory
        'com.pdb',         # Complex topology
        'apo.pdb',         # Protein topology  
        'lig.pdb'          # Ligand topology
    ]
    
    missing_files = []
    for file in required_files:
        if not os.path.exists(os.path.join(work_dir, file)):
            missing_files.append(file)
    
    if missing_files:
        print(f"Missing required files: {missing_files}")
        return None, None
    
    # Create calculator
    calculator = OpenMMGBSACalculator(work_dir, verbose)
    
    # Set up file paths
    complex_traj = os.path.join(work_dir, 'trajectory.dcd')
    complex_top = os.path.join(work_dir, 'com.pdb')
    
    # For protein and ligand, we'll use the same trajectory for now
    # In a full implementation, you would have separate trajectories
    protein_traj = complex_traj
    ligand_traj = complex_traj
    protein_top = os.path.join(work_dir, 'apo.pdb')
    ligand_top = os.path.join(work_dir, 'lig.pdb')
    
    # Run calculation
    if use_pbsa:
        return calculator.calculate_pbsa(
            complex_traj, protein_traj, ligand_traj,
            complex_top, protein_top, ligand_top
        )
    else:
        return calculator.calculate_mmgbsa(
            complex_traj, protein_traj, ligand_traj,
            complex_top, protein_top, ligand_top,
            igb=igb
        ) 