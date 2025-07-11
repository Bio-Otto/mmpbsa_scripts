import argparse
import os
import subprocess
import tempfile
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import numpy as np


def split_complex_from_system(complex_file):
    """
    Split complex PDB into protein and ligand using RDKit
    """
    # Read the PDB file
    mol = Chem.MolFromPDBFile(complex_file, removeHs=False)
    if mol is None:
        print("Could not read PDB file!")
        exit()
    
    # Split the molecule into fragments
    fragments = Chem.GetMolFrags(mol, asMols=True)
    
    if len(fragments) < 2:
        print("Could not identify separate protein and ligand fragments!")
        exit()
    
    # Assume the largest fragment is the protein and the smallest is the ligand
    # This is a simple heuristic - you might need to refine this logic
    fragments = sorted(fragments, key=lambda x: x.GetNumAtoms(), reverse=True)
    protein = fragments[0]  # Largest fragment
    ligand = fragments[1]   # Second largest fragment
    
    return protein, ligand


def prepare_receptor(protein, output_dir):
    """
    Prepare receptor for docking using Vina
    """
    # Save protein to PDB
    protein_file = os.path.join(output_dir, "receptor.pdb")
    Chem.MolToPDBFile(protein, protein_file)
    
    # Convert to PDBQT format for Vina (requires Open Babel)
    receptor_pdbqt = os.path.join(output_dir, "receptor.pdbqt")
    try:
        subprocess.run([
            "obabel", protein_file, "-O", receptor_pdbqt, 
            "--addhydrogens", "--assigncharges", "gasteiger"
        ], check=True, capture_output=True)
    except subprocess.CalledProcessError:
        print("Warning: Could not convert receptor to PDBQT format")
        return None
    
    return receptor_pdbqt


def prepare_ligand(ligand, output_dir):
    """
    Prepare ligand for docking using Vina
    """
    # Generate 3D coordinates if needed
    if ligand.GetNumConformers() == 0:
        AllChem.EmbedMolecule(ligand, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(ligand)
    
    # Save ligand to MOL2
    ligand_file = os.path.join(output_dir, "ligand.mol2")
    Chem.MolToMolFile(ligand, ligand_file)
    
    # Convert to PDBQT format for Vina
    ligand_pdbqt = os.path.join(output_dir, "ligand.pdbqt")
    try:
        subprocess.run([
            "obabel", ligand_file, "-O", ligand_pdbqt,
            "--addhydrogens", "--assigncharges", "gasteiger"
        ], check=True, capture_output=True)
    except subprocess.CalledProcessError:
        print("Warning: Could not convert ligand to PDBQT format")
        return None
    
    return ligand_pdbqt


def score_pose_vina(receptor_pdbqt, ligand_pdbqt, output_dir):
    """
    Score ligand pose using Vina
    """
    if receptor_pdbqt is None or ligand_pdbqt is None:
        return {"Vina Score": 0.0, "Error": "Could not prepare files for docking"}
    
    # Run Vina scoring
    output_file = os.path.join(output_dir, "vina_out.pdbqt")
    log_file = os.path.join(output_dir, "vina.log")
    
    try:
        result = subprocess.run([
            "vina", "--receptor", receptor_pdbqt, "--ligand", ligand_pdbqt,
            "--out", output_file, "--score_only", "--log", log_file
        ], check=True, capture_output=True, text=True)
        
        # Parse Vina score from log file
        vina_score = 0.0
        if os.path.exists(log_file):
            with open(log_file, 'r') as f:
                for line in f:
                    if "REMARK VINA RESULT" in line:
                        vina_score = float(line.split()[3])
                        break
        
        return {"Vina Score": vina_score}
        
    except subprocess.CalledProcessError as e:
        return {"Vina Score": 0.0, "Error": f"Vina failed: {e.stderr}"}


def calculate_molecular_properties(ligand):
    """
    Calculate molecular properties using RDKit
    """
    properties = {}
    
    # Basic properties
    properties['Molecular Weight'] = rdMolDescriptors.CalcExactMolWt(ligand)
    properties['LogP'] = rdMolDescriptors.CalcCrippenDescriptors(ligand)[0]
    properties['TPSA'] = rdMolDescriptors.CalcTPSA(ligand)
    properties['HBD'] = rdMolDescriptors.CalcNumHBD(ligand)
    properties['HBA'] = rdMolDescriptors.CalcNumHBA(ligand)
    properties['Rotatable Bonds'] = rdMolDescriptors.CalcNumRotatableBonds(ligand)
    properties['Aromatic Rings'] = rdMolDescriptors.CalcNumAromaticRings(ligand)
    
    # Lipinski's Rule of Five
    mw = properties['Molecular Weight']
    logp = properties['LogP']
    hbd = properties['HBD']
    hba = properties['HBA']
    
    lipinski_violations = 0
    if mw > 500: lipinski_violations += 1
    if logp > 5: lipinski_violations += 1
    if hbd > 5: lipinski_violations += 1
    if hba > 10: lipinski_violations += 1
    
    properties['Lipinski Violations'] = lipinski_violations
    
    return properties


def score_and_build(complex_pdb):
    """
    Main function to score complex and calculate properties
    """
    # Create temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as temp_dir:
        
        # Split complex
        protein, ligand = split_complex_from_system(complex_pdb)
        
        # Prepare files for docking
        receptor_pdbqt = prepare_receptor(protein, temp_dir)
        ligand_pdbqt = prepare_ligand(ligand, temp_dir)
        
        # Score with Vina
        vina_scores = score_pose_vina(receptor_pdbqt, ligand_pdbqt, temp_dir)
        
        # Calculate molecular properties
        mol_properties = calculate_molecular_properties(ligand)
        
        # Combine results
        results = {**vina_scores, **mol_properties}
        
        return results


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--complex', type=str, default='/Users/austin/MPro-docked/MPro_0387_Gen3L_21.pdb', required=False)
    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()
    scores = score_and_build(args.complex)
    print(scores)
