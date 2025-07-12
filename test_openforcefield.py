#!/usr/bin/env python3
"""
Test script for OpenForceField-based MMGBSA calculation
"""

import sys
import os
sys.path.append('.')

from mmgpbsa.openforcefield_mmgbsa import run_openforcefield_mmgbsa

def main():
    """Test OpenForceField-based MMGBSA calculation"""
    
    complex_pdb = "../complex.pdb"
    trajectory_file = "../complex.xtc"
    
    print("Testing OpenForceField-based MMGBSA calculation...")
    print(f"Complex PDB: {complex_pdb}")
    print(f"Trajectory: {trajectory_file}")
    
    # Run with verbose=2 for maximum debug output
    result = run_openforcefield_mmgbsa(
        complex_pdb=complex_pdb,
        trajectory_file=trajectory_file,
        method='gbsa',
        igb=5,
        start_frame=0,
        end_frame=10,  # Just test first 10 frames
        stride=1,
        verbose=2
    )
    
    if result[0] is not None:
        print(f"Success! ΔG = {result[0]:.4f} ± {result[1]:.4f} kcal/mol")
    else:
        print("Calculation failed")

if __name__ == "__main__":
    main() 