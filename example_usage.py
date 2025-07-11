#!/usr/bin/env python3
"""
Example usage of the Open-Source MMGBSA Package
"""

import os
import tempfile
from pathlib import Path

def create_example_molecule():
    """
    Create a simple example molecule for testing
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    
    # Create a simple molecule (ethanol)
    mol = Chem.MolFromSmiles("CCO")
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    
    return mol

def create_example_protein():
    """
    Create a simple example protein structure
    """
    from rdkit import Chem
    
    # Create a simple peptide (alanine dipeptide)
    mol = Chem.MolFromSmiles("CC(=O)NC(C)C(=O)N")
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    
    return mol

def create_example_complex():
    """
    Create an example protein-ligand complex
    """
    protein = create_example_protein()
    ligand = create_example_molecule()
    
    # Combine them (simplified approach)
    from rdkit import Chem
    complex_mol = Chem.CombineMols(protein, ligand)
    
    return complex_mol

def run_basic_example():
    """
    Run a basic MMGBSA example
    """
    print("Running basic MMGBSA example...")
    
    # Create temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Create example structures
        complex_mol = create_example_complex()
        protein_mol = create_example_protein()
        ligand_mol = create_example_molecule()
        
        # Save structures
        complex_file = temp_path / "complex.pdb"
        protein_file = temp_path / "protein.pdb"
        ligand_file = temp_path / "ligand.mol2"
        
        from rdkit import Chem
        Chem.MolToPDBFile(complex_mol, str(complex_file))
        Chem.MolToPDBFile(protein_mol, str(protein_file))
        Chem.MolToMolFile(ligand_mol, str(ligand_file))
        
        print(f"Created example files in: {temp_dir}")
        print(f"  Complex: {complex_file}")
        print(f"  Protein: {protein_file}")
        print(f"  Ligand: {ligand_file}")
        
        # Run MMGBSA analysis
        try:
            from mmgpbsa.systemloader import SystemLoader
            
            # Initialize system loader
            sl = SystemLoader(
                dirpath=str(temp_path / "results"),
                verbose=1,
                input_pdb=str(complex_file),
                ps=5.0,  # Short simulation for example
                calcFrames=5,
                platform_name="CPU"
            )
            
            print("Preparing simulation...")
            sl.prepare_simulation()
            
            print("Running MMGBSA calculation...")
            delta_g, std = sl.run_amber("gbsa", amber_path=None)
            
            if delta_g is not None:
                print(f"MMGBSA Result: ΔG = {delta_g:.2f} ± {std:.2f} kcal/mol")
            else:
                print("MMGBSA calculation failed")
                
        except Exception as e:
            print(f"Error running MMGBSA: {e}")
            print("This might be due to missing dependencies or configuration issues.")

def run_component_scoring_example():
    """
    Run component scoring example
    """
    print("\nRunning component scoring example...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Create example complex
        complex_mol = create_example_complex()
        complex_file = temp_path / "complex.pdb"
        
        from rdkit import Chem
        Chem.MolToPDBFile(complex_mol, str(complex_file))
        
        try:
            from component_scoring import score_and_build
            
            scores = score_and_build(str(complex_file))
            print("Component scores:")
            for key, value in scores.items():
                print(f"  {key}: {value}")
                
        except Exception as e:
            print(f"Error running component scoring: {e}")

def main():
    """
    Main function to run examples
    """
    print("Open-Source MMGBSA Package - Example Usage")
    print("=" * 50)
    
    # Check if required packages are available
    try:
        import numpy as np
        import mdtraj as md
        from rdkit import Chem
        import openmm as mm
        from openmm import unit, app
        from pymbar import timeseries
        print("✓ All required packages are available")
    except ImportError as e:
        print(f"✗ Missing required package: {e}")
        print("Please install all dependencies: pip install -r requirements.txt")
        return
    
    # Run examples
    run_basic_example()
    run_component_scoring_example()
    
    print("\n" + "=" * 50)
    print("Example completed!")
    print("\nTo run your own analysis:")
    print("  python run_mmgpbsa.py --com your_complex.pdb --odir results")
    print("  python component_scoring.py --complex your_complex.pdb")

if __name__ == "__main__":
    main() 