import argparse
import logging
import os
import shutil
import warnings

from mmgpbsa.xtc_mmgbsa import run_xtc_mmgbsa


def get_args():
    parser = argparse.ArgumentParser(
        description='Open-source MMGBSA analysis package',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with complex PDB and trajectory
  python run_mmgpbsa.py --com complex.pdb --traj trajectory.xtc --odir results

  # Using GBSA method with specific parameters
  python run_mmgpbsa.py --com complex.pdb --traj trajectory.xtc --method gbsa --odir results

  # Using PBSA method (simplified implementation)
  python run_mmgpbsa.py --com complex.pdb --traj trajectory.xtc --method pbsa --odir results

  # Advanced usage with frame selection
  python run_mmgpbsa.py --com complex.pdb --traj trajectory.xtc --calcFrames 100 --mbar 5 --odir results
        """
    )
    
    # Input file options
    input_group = parser.add_argument_group('Input Files')
    input_group.add_argument('--com', type=str, default=None,
                            help='PDB file with ligand and protein in complex')
    input_group.add_argument('--traj', type=str, default=None,
                            help='Trajectory file (XTC, DCD, TRR format)')
    input_group.add_argument('--apo', type=str, default=None,
                            help='PDB file without ligand (optional, will be extracted from complex)')
    input_group.add_argument('--lig', type=str, default=None,
                            help='Ligand file (optional, will be extracted from complex)')
    
    # Analysis options
    analysis_group = parser.add_argument_group('Analysis Options')
    analysis_group.add_argument('--calcFrames', type=int, default=None,
                               help='Frame stride for analysis (default: 1, analyze all frames)')
    analysis_group.add_argument('--mbar', type=int, default=5,
                               help='GB model: 5=GBn, 7=GBn2, 8=OBC (default: 5)')
    analysis_group.add_argument('--method', type=str, choices=['gbsa', 'pbsa'], default='gbsa',
                               help='Use PBSA or GBSA (default: gbsa)')
    
    # Output options
    output_group = parser.add_argument_group('Output Options')
    output_group.add_argument('--odir', type=str, default=None,
                             help='Output directory for results (default: auto-generated)')
    output_group.add_argument('-v', '--verbose', type=int, choices=[0, 1, 2], default=1,
                             help='Verbosity level (default: 1)', dest='v')
    
    args = parser.parse_args()
    
    # Validate input arguments
    if args.com is None:
        parser.error("--com (complex PDB file) must be provided")
    if args.traj is None:
        parser.error("--traj (trajectory file) must be provided")
    
    return args

def setup_folder(args):
    if args.odir is None:
        # Use a default directory name based on complex file
        if args.com is not None:
            args.odir = f'{os.getcwd()}/{args.com.split("/")[-1].split(".")[0]}_mmgbsa'
        else:
            args.odir = f'{os.getcwd()}/mmgbsa_analysis'

    try:
        os.mkdir(args.odir)
    except OSError:
        if not os.path.isdir(args.odir):
            print("Creation of the directory %s failed" % args.odir)
            exit()
        else:
            shutil.rmtree(args.odir)
            os.mkdir(args.odir)
    else:
        if args.v >= 1:
            print("Successfully created the directory %s " % args.odir)

if __name__ == '__main__':
    args = get_args()

    setup_folder(args)

    # No AMBER dependency needed for XTC-based analysis

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.ERROR)
    logging.getLogger('openmm').setLevel(logging.ERROR)
    logging.getLogger('openforcefield').setLevel(logging.ERROR)
    logging.getLogger('openmmtools').setLevel(logging.ERROR)
    warnings.filterwarnings("ignore")

    # Run MMGBSA calculation using XTC trajectory
    print(f"Running MMGBSA calculation...")
    print(f"Complex PDB: {args.com}")
    print(f"Trajectory: {args.traj}")
    print(f"Method: {args.method}")
    print(f"Output directory: {args.odir}")
    
    # Calculate frame range if specified
    start_frame = 0
    end_frame = None
    stride = 1
    
    if args.calcFrames is not None:
        # If calcFrames is specified, use it as stride to subsample
        stride = max(1, args.calcFrames)
    
    deltag, std = run_xtc_mmgbsa(
        complex_pdb=args.com,
        trajectory_file=args.traj,
        method=args.method,
        igb=args.mbar,  # Use mbar parameter as igb value
        start_frame=start_frame,
        end_frame=end_frame,
        stride=stride,
        verbose=args.v
    )
    
    if deltag is not None:
        print(f"MMGBSA Result: {deltag:.4f} ± {std:.4f} kcal/mol")
    else:
        print("MMGBSA calculation failed")