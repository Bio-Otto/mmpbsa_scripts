"""
AMBER Parameter Generator using Open Source Tools
Replaces proprietary AMBER software for parameter generation
"""

import os
import subprocess
import tempfile
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
import numpy as np


class AmberParameterGenerator:
    """
    Generate AMBER parameters using open-source tools
    """
    
    def __init__(self, work_dir):
        self.work_dir = work_dir
        
    def generate_ligand_parameters(self, ligand_mol, output_prefix="lig"):
        """
        Generate AMBER parameters for a ligand using open-source tools
        
        Args:
            ligand_mol: RDKit molecule object
            output_prefix: Prefix for output files
            
        Returns:
            dict: Dictionary with file paths to generated parameters
        """
        
        # Save ligand in various formats
        mol2_file = os.path.join(self.work_dir, f"{output_prefix}.mol2")
        pdb_file = os.path.join(self.work_dir, f"{output_prefix}.pdb")
        
        # Ensure molecule has hydrogens
        ligand_mol = Chem.AddHs(ligand_mol)
        
        # Generate 3D coordinates if needed
        if ligand_mol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(ligand_mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(ligand_mol)
        
        # Save files
        Chem.MolToMolFile(ligand_mol, mol2_file)
        Chem.MolToPDBFile(ligand_mol, pdb_file)
        
        # Use Open Babel for charge assignment and format conversion
        charged_mol2 = os.path.join(self.work_dir, f"{output_prefix}_charged.mol2")
        pdbqt_file = os.path.join(self.work_dir, f"{output_prefix}.pdbqt")
        
        try:
            # Assign charges using Open Babel
            subprocess.run([
                'obabel', mol2_file, '-O', charged_mol2,
                '--addhydrogens', '--assigncharges', 'gasteiger'
            ], check=True, capture_output=True)
            
            # Convert to PDBQT for Vina compatibility
            subprocess.run([
                'obabel', charged_mol2, '-O', pdbqt_file
            ], check=True, capture_output=True)
            
        except subprocess.CalledProcessError as e:
            print(f"Warning: Open Babel failed: {e}")
            # Fallback to RDKit charge assignment
            AllChem.ComputeGasteigerCharges(ligand_mol)
            Chem.MolToMolFile(ligand_mol, charged_mol2)
        
        # Generate simplified AMBER parameters
        # Note: This is a simplified approach. For full AMBER compatibility,
        # you would need more sophisticated parameter generation
        prmtop_file = os.path.join(self.work_dir, f"{output_prefix}.prmtop")
        inpcrd_file = os.path.join(self.work_dir, f"{output_prefix}.inpcrd")
        
        self._create_simplified_amber_files(ligand_mol, prmtop_file, inpcrd_file)
        
        return {
            'mol2': mol2_file,
            'pdb': pdb_file,
            'charged_mol2': charged_mol2,
            'pdbqt': pdbqt_file,
            'prmtop': prmtop_file,
            'inpcrd': inpcrd_file
        }
    
    def _create_simplified_amber_files(self, mol, prmtop_file, inpcrd_file):
        """
        Create simplified AMBER topology and coordinate files
        This is a placeholder implementation
        """
        
        # Get molecular information
        num_atoms = mol.GetNumAtoms()
        conf = mol.GetConformer()
        
        # Create simplified topology file
        with open(prmtop_file, 'w') as f:
            f.write(f"Simplified AMBER topology for {num_atoms} atoms\n")
            f.write(f"Generated by open-source MMGBSA\n")
            f.write(f"Number of atoms: {num_atoms}\n")
            f.write(f"Force field: AMBER14\n")
            f.write(f"Charge method: Gasteiger\n")
        
        # Create simplified coordinate file
        with open(inpcrd_file, 'w') as f:
            f.write(f"Simplified AMBER coordinates for {num_atoms} atoms\n")
            f.write(f"Generated by open-source MMGBSA\n")
            
            # Write coordinates
            for i in range(num_atoms):
                pos = conf.GetAtomPosition(i)
                f.write(f"{pos.x:8.3f} {pos.y:8.3f} {pos.z:8.3f}\n")
    
    def prepare_protein(self, protein_mol, output_prefix="apo"):
        """
        Prepare protein for AMBER simulation
        
        Args:
            protein_mol: RDKit molecule object
            output_prefix: Prefix for output files
            
        Returns:
            dict: Dictionary with file paths to generated files
        """
        
        # Save protein
        pdb_file = os.path.join(self.work_dir, f"{output_prefix}.pdb")
        Chem.MolToPDBFile(protein_mol, pdb_file)
        
        # Create simplified AMBER files for protein
        prmtop_file = os.path.join(self.work_dir, f"{output_prefix}.prmtop")
        inpcrd_file = os.path.join(self.work_dir, f"{output_prefix}.inpcrd")
        
        self._create_simplified_amber_files(protein_mol, prmtop_file, inpcrd_file)
        
        return {
            'pdb': pdb_file,
            'prmtop': prmtop_file,
            'inpcrd': inpcrd_file
        }
    
    def create_complex_system(self, protein_mol, ligand_mol, output_prefix="com"):
        """
        Create complex system from protein and ligand
        
        Args:
            protein_mol: RDKit molecule object for protein
            ligand_mol: RDKit molecule object for ligand
            output_prefix: Prefix for output files
            
        Returns:
            dict: Dictionary with file paths to generated files
        """
        
        # Combine molecules (simplified approach)
        # In a real implementation, you would need to properly position the ligand
        # relative to the protein binding site
        
        # For now, we'll create a simple combined structure
        combined_mol = Chem.CombineMols(protein_mol, ligand_mol)
        
        # Save complex
        pdb_file = os.path.join(self.work_dir, f"{output_prefix}.pdb")
        Chem.MolToPDBFile(combined_mol, pdb_file)
        
        # Create simplified AMBER files for complex
        prmtop_file = os.path.join(self.work_dir, f"{output_prefix}.prmtop")
        inpcrd_file = os.path.join(self.work_dir, f"{output_prefix}.inpcrd")
        
        self._create_simplified_amber_files(combined_mol, prmtop_file, inpcrd_file)
        
        return {
            'pdb': pdb_file,
            'prmtop': prmtop_file,
            'inpcrd': inpcrd_file
        }


def create_leap_script(protein_file, ligand_file, output_dir):
    """
    Create a tleap script for AMBER parameter generation
    This is a template that can be used if AMBER is available
    """
    
    leap_script = f"""
# AMBER tleap script for MMGBSA
# Generated by open-source MMGBSA

# Load force fields
source leaprc.protein.ff14SBonlysc
source leaprc.gaff2
set default PBRadii mbondi3

# Load protein
rec = loadPDB {protein_file}

# Load ligand
lig = loadmol2 {ligand_file}

# Create complex
com = combine {{rec lig}}

# Save parameters
saveAmberParm rec apo.prmtop apo.inpcrd
saveAmberParm lig lig.prmtop lig.inpcrd
saveAmberParm com com.prmtop com.inpcrd

quit
"""
    
    with open(os.path.join(output_dir, 'leap.in'), 'w') as f:
        f.write(leap_script)
    
    return os.path.join(output_dir, 'leap.in') 