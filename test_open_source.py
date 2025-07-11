#!/usr/bin/env python3
"""
Test script to verify open-source MMGBSA package functionality
"""

import os
import sys
import tempfile
from pathlib import Path

def test_imports():
    """Test that all required packages can be imported"""
    print("Testing imports...")
    
    try:
        import numpy as np
        print("✓ NumPy imported successfully")
    except ImportError as e:
        print(f"✗ NumPy import failed: {e}")
        return False
    
    try:
        import mdtraj as md
        print("✓ MDTraj imported successfully")
    except ImportError as e:
        print(f"✗ MDTraj import failed: {e}")
        return False
    
    try:
        from rdkit import Chem
        print("✓ RDKit imported successfully")
    except ImportError as e:
        print(f"✗ RDKit import failed: {e}")
        return False
    
    try:
        import simtk.openmm as mm
        from simtk import unit
        from simtk.openmm import app
        print("✓ OpenMM imported successfully")
    except ImportError as e:
        print(f"✗ OpenMM import failed: {e}")
        return False
    
    try:
        from pymbar import timeseries
        print("✓ PyMBAR imported successfully")
    except ImportError as e:
        print(f"✗ PyMBAR import failed: {e}")
        return False
    
    try:
        from mmgpbsa.systemloader import SystemLoader
        from mmgpbsa.parameter_generator import AmberParameterGenerator
        from mmgpbsa.openmm_mmgbsa import OpenMMGBSACalculator
        print("✓ MMGBSA modules imported successfully")
    except ImportError as e:
        print(f"✗ MMGBSA modules import failed: {e}")
        return False
    
    return True

def test_rdkit_functionality():
    """Test RDKit functionality"""
    print("\nTesting RDKit functionality...")
    
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        # Create a simple molecule
        mol = Chem.MolFromSmiles("CCO")
        if mol is None:
            print("✗ Failed to create molecule from SMILES")
            return False
        
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates
        AllChem.EmbedMolecule(mol, randomSeed=42)
        
        print("✓ RDKit basic functionality works")
        return True
        
    except Exception as e:
        print(f"✗ RDKit functionality test failed: {e}")
        return False

def test_openmm_functionality():
    """Test OpenMM functionality"""
    print("\nTesting OpenMM functionality...")
    
    try:
        import simtk.openmm as mm
        from simtk import unit
        from simtk.openmm import app
        
        # Create a simple system
        topology = app.Topology()
        
        # Add a simple molecule (methane)
        chain = topology.addChain()
        residue = topology.addResidue("MOL", chain)
        
        # Add atoms
        c1 = topology.addAtom("C", app.Element.getBySymbol("C"), residue)
        h1 = topology.addAtom("H", app.Element.getBySymbol("H"), residue)
        h2 = topology.addAtom("H", app.Element.getBySymbol("H"), residue)
        h3 = topology.addAtom("H", app.Element.getBySymbol("H"), residue)
        h4 = topology.addAtom("H", app.Element.getBySymbol("H"), residue)
        
        # Add bonds
        topology.addBond(c1, h1)
        topology.addBond(c1, h2)
        topology.addBond(c1, h3)
        topology.addBond(c1, h4)
        
        # Create force field
        forcefield = app.ForceField('amber14-all.xml')
        
        # Create system
        system = forcefield.createSystem(
            topology,
            nonbondedMethod=app.CutoffNonPeriodic,
            nonbondedCutoff=1.0*unit.nanometer,
            constraints=app.HBonds,
            implicitSolvent=app.GBn2,
            removeCMMotion=True
        )
        
        print("✓ OpenMM basic functionality works")
        return True
        
    except Exception as e:
        print(f"✗ OpenMM functionality test failed: {e}")
        return False

def test_parameter_generator():
    """Test parameter generator"""
    print("\nTesting parameter generator...")
    
    try:
        from mmgpbsa.parameter_generator import AmberParameterGenerator
        from rdkit import Chem
        
        # Create a simple molecule
        mol = Chem.MolFromSmiles("CCO")
        mol = Chem.AddHs(mol)
        
        # Create temporary directory
        with tempfile.TemporaryDirectory() as temp_dir:
            param_gen = AmberParameterGenerator(temp_dir)
            
            # Generate parameters
            params = param_gen.generate_ligand_parameters(mol, "test")
            
            # Check that files were created
            required_files = ['mol2', 'pdb', 'charged_mol2', 'prmtop', 'inpcrd']
            for file_type in required_files:
                if file_type in params:
                    file_path = params[file_type]
                    if os.path.exists(file_path):
                        print(f"✓ Generated {file_type} file")
                    else:
                        print(f"✗ Failed to generate {file_type} file")
                        return False
                else:
                    print(f"✗ Missing {file_type} in parameters")
                    return False
        
        print("✓ Parameter generator works")
        return True
        
    except Exception as e:
        print(f"✗ Parameter generator test failed: {e}")
        return False

def test_external_dependencies():
    """Test external dependencies (Open Babel, Vina)"""
    print("\nTesting external dependencies...")
    
    # Test Open Babel
    try:
        import subprocess
        result = subprocess.run(['obabel', '--version'], 
                              capture_output=True, text=True)
        if result.returncode == 0:
            print("✓ Open Babel is available")
        else:
            print("⚠ Open Babel not available (optional)")
    except FileNotFoundError:
        print("⚠ Open Babel not found (optional)")
    
    # Test Vina
    try:
        result = subprocess.run(['vina', '--version'], 
                              capture_output=True, text=True)
        if result.returncode == 0:
            print("✓ Vina is available")
        else:
            print("⚠ Vina not available (optional)")
    except FileNotFoundError:
        print("⚠ Vina not found (optional)")
    
    return True

def main():
    """Run all tests"""
    print("Testing Open-Source MMGBSA Package")
    print("=" * 40)
    
    tests = [
        test_imports,
        test_rdkit_functionality,
        test_openmm_functionality,
        test_parameter_generator,
        test_external_dependencies
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        try:
            if test():
                passed += 1
        except Exception as e:
            print(f"✗ Test {test.__name__} failed with exception: {e}")
    
    print("\n" + "=" * 40)
    print(f"Tests passed: {passed}/{total}")
    
    if passed == total:
        print("🎉 All tests passed! The open-source MMGBSA package is ready to use.")
        return 0
    else:
        print("⚠ Some tests failed. Please check the installation.")
        return 1

if __name__ == "__main__":
    sys.exit(main()) 