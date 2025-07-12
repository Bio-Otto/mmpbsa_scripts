"""
Simplified MMGBSA Calculator
Uses a more robust approach for handling unknown residues
"""

import os
import numpy as np
import mdtraj as md
from openmm import unit
from openmm import app


class SimpleMMGBSACalculator:
    """
    Simplified MMGBSA calculator that handles unknown residues better
    """
    
    def __init__(self, verbose=1):
        self.verbose = verbose
        
    def calculate_mmgbsa(self, complex_pdb, trajectory_file, method='gbsa', 
                        igb=5, start_frame=0, end_frame=None, stride=1):
        """
        Calculate MMGBSA using a simplified approach
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
            
            # Calculate energies using simplified approach
            complex_energies = self._calculate_simple_energies(traj, "complex")
            protein_energies = self._calculate_simple_energies(protein_traj, "protein")
            ligand_energies = self._calculate_simple_energies(ligand_traj, "ligand")
            
            # Debug information
            if self.verbose >= 1:
                print(f"Energy calculation results:")
                print(f"  Complex energies: {len(complex_energies)} frames")
                print(f"  Protein energies: {len(protein_energies)} frames")
                print(f"  Ligand energies: {len(ligand_energies)} frames")
            
            # Calculate MMGBSA
            if len(complex_energies) == 0 or len(protein_energies) == 0 or len(ligand_energies) == 0:
                print("Error: Could not calculate energies for all systems")
                if len(complex_energies) == 0:
                    print("  - Complex energy calculation failed")
                if len(protein_energies) == 0:
                    print("  - Protein energy calculation failed")
                if len(ligand_energies) == 0:
                    print("  - Ligand energy calculation failed")
                return None, None
            
            # Calculate binding free energy
            # ΔG = G_complex - G_protein - G_ligand
            binding_energies = complex_energies - protein_energies - ligand_energies
            
            # Convert to kcal/mol (assuming energies are in kJ/mol)
            binding_energies_kcal = binding_energies * 0.593
            
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
        """Split trajectory into protein and ligand parts"""
        if self.verbose >= 2:
            print("Splitting trajectory into protein and ligand...")
        
        # Get residue information
        residues = [res for res in traj.topology.residues]
        
        if len(residues) < 2:
            print("Error: Need at least 2 residues to split protein and ligand")
            return traj, traj
        
        # Find the smallest residue (assumed to be ligand)
        residue_sizes = [(res, len(list(res.atoms))) for res in residues]
        residue_sizes.sort(key=lambda x: x[1])
        
        ligand_residue = residue_sizes[0][0]  # Smallest residue
        protein_residues = [res for res in residues if res != ligand_residue]
        
        if self.verbose >= 2:
            print(f"Ligand residue: {ligand_residue.name} ({len(list(ligand_residue.atoms))} atoms)")
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
        
        if self.verbose >= 1:
            print(f"Trajectory splitting complete:")
            print(f"  Original trajectory: {len(list(traj.topology.atoms))} atoms, {len(list(traj.topology.residues))} residues")
            print(f"  Protein trajectory: {len(list(protein_traj.topology.atoms))} atoms, {len(list(protein_traj.topology.residues))} residues")
            print(f"  Ligand trajectory: {len(list(ligand_traj.topology.atoms))} atoms, {len(list(ligand_traj.topology.residues))} residues")
        
        return protein_traj, ligand_traj
    
    def _calculate_simple_energies(self, traj, system_type):
        """Calculate energies using a simplified approach"""
        energies = []
        
        for i, frame in enumerate(traj):
            if self.verbose >= 2 and i % 10 == 0:
                print(f"Processing {system_type} frame {i+1}/{len(traj)}")
            
            energy = self._calculate_simple_frame_energy(frame, system_type)
            if energy is not None:
                energies.append(energy)
        
        return np.array(energies)
    
    def _calculate_simple_frame_energy(self, frame, system_type):
        """Calculate energy using a very simple approach"""
        try:
            # Convert frame to OpenMM system
            topology = frame.topology.to_openmm()
            positions = frame.xyz[0] * unit.nanometers
            
            if self.verbose >= 2:
                print(f"  {system_type}: {len(list(topology.atoms()))} atoms, {len(list(topology.residues()))} residues")
            
            # Use a very simple force field approach
            # This avoids issues with unknown residues like "LIG"
            
            # Create system with minimal forces
            system = app.ForceField().createSystem(
                topology,
                nonbondedMethod=app.NoCutoff,
                constraints=None,
                removeCMMotion=False
            )
            
            # Add a simple GB force if possible
            try:
                gb_force = app.GBSAOBCForce()
                for atom in topology.atoms():
                    # Use default parameters for all atoms
                    gb_force.addParticle(0.0, 0.0, 0.0)
                system.addForce(gb_force)
            except:
                # If GB force fails, continue without it
                pass
            
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
            if self.verbose >= 1:
                print(f"Error calculating energy for {system_type}: {e}")
                print(f"  Frame has {len(list(frame.topology.atoms))} atoms, {len(list(frame.topology.residues))} residues")
            return None


def run_simple_mmgbsa(complex_pdb, trajectory_file, method='gbsa', igb=5, 
                      start_frame=0, end_frame=None, stride=1, verbose=1):
    """
    Main function to run simplified MMGBSA calculation
    """
    
    # Check if files exist
    if not os.path.exists(complex_pdb):
        print(f"Error: Complex PDB file not found: {complex_pdb}")
        return None, None
    
    if not os.path.exists(trajectory_file):
        print(f"Error: Trajectory file not found: {trajectory_file}")
        return None, None
    
    # Create calculator
    calculator = SimpleMMGBSACalculator(verbose)
    
    # Run calculation
    return calculator.calculate_mmgbsa(
        complex_pdb, trajectory_file, method=method, igb=igb,
        start_frame=start_frame, end_frame=end_frame, stride=stride
    ) 