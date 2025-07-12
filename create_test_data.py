#!/usr/bin/env python3
"""
Create test data for XTC-based MMGBSA testing
"""

import numpy as np
import mdtraj as md
from rdkit import Chem
from rdkit.Chem import AllChem

def create_test_complex():
    """Create a simple protein-ligand complex"""
    
    # Create a simple ligand (ethanol)
    ligand = Chem.MolFromSmiles("CCO")
    ligand = Chem.AddHs(ligand)
    AllChem.EmbedMolecule(ligand, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(ligand)
    
    # Create a simple protein (alanine dipeptide)
    protein = Chem.MolFromSmiles("CC(=O)NC(C)C(=O)N")
    protein = Chem.AddHs(protein)
    AllChem.EmbedMolecule(protein, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(protein)
    
    # Combine them
    complex_mol = Chem.CombineMols(protein, ligand)
    
    return complex_mol

def create_test_trajectory(complex_mol, n_frames=10):
    """Create a simple test trajectory"""
    
    # Convert RDKit molecule to MDTraj
    # First save as PDB
    Chem.MolToPDBFile(complex_mol, "temp_complex.pdb")
    
    # Load with MDTraj
    traj = md.load("temp_complex.pdb")
    
    # Create a simple trajectory by adding small random displacements
    frames = []
    for i in range(n_frames):
        # Add small random displacement to coordinates
        new_xyz = traj.xyz[0] + np.random.normal(0, 0.01, traj.xyz[0].shape)
        frames.append(new_xyz)
    
    # Create new trajectory
    new_traj = md.Trajectory(
        xyz=np.array(frames),
        topology=traj.topology
    )
    
    # Clean up temp file
    import os
    os.remove("temp_complex.pdb")
    
    return new_traj

def main():
    """Create test data"""
    print("Creating test data for MMGBSA...")
    
    # Create complex
    complex_mol = create_test_complex()
    print(f"Created complex with {complex_mol.GetNumAtoms()} atoms")
    
    # Save complex PDB
    Chem.MolToPDBFile(complex_mol, "test_complex.pdb")
    print("Saved test_complex.pdb")
    
    # Create trajectory
    traj = create_test_trajectory(complex_mol, n_frames=5)
    print(f"Created trajectory with {traj.n_frames} frames")
    
    # Save trajectory as XTC
    traj.save_xtc("test_trajectory.xtc")
    print("Saved test_trajectory.xtc")
    
    print("Test data created successfully!")
    print("You can now run:")
    print("python run_mmgpbsa.py --com test_complex.pdb --traj test_trajectory.xtc --odir results")

if __name__ == "__main__":
    main() 