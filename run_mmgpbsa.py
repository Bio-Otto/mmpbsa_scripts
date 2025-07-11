import argparse
import logging
import os
import shutil
import warnings

from mmgpbsa.amber_mmgpbsa import run_amber
from mmgpbsa.systemloader import SystemLoader


def get_args():
    parser = argparse.ArgumentParser(
        description='Open-source MMGBSA analysis package',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with complex file
  python run_mmgpbsa.py --com complex.pdb --odir results

  # Using separate protein and ligand files
  python run_mmgpbsa.py --apo protein.pdb --lig ligand.mol2 --odir results

  # Advanced usage with specific parameters
  python run_mmgpbsa.py --com complex.pdb --odir results --platform CUDA --ps 1000 --mbar 50
        """
    )
    
    # Input file options
    input_group = parser.add_argument_group('Input Files')
    input_group.add_argument('--com', type=str, default=None,
                            help='PDB file with ligand and protein in complex (single chain)')
    input_group.add_argument('--apo', type=str, default=None,
                            help='PDB file without ligand (single chain)')
    input_group.add_argument('--lig', type=str, default=None,
                            help='Ligand file (MOL2, SDF, PDB format)')
    
    # AMBER options (optional)
    amber_group = parser.add_argument_group('AMBER Options (Optional)')
    amber_group.add_argument('--amber_path', type=str, default=None,
                            help='Path to AMBER installation folder (e.g., /home/user/amber20)')
    
    # Simulation options
    sim_group = parser.add_argument_group('Simulation Options')
    sim_group.add_argument('--platform', type=str, choices=['CPU', 'CUDA', 'OpenCL'], default='CPU',
                          help='OpenMM platform to use (default: CPU)')
    sim_group.add_argument('--ps', type=float, default=10.0,
                          help='Picoseconds to run simulation (default: 10.0)')
    sim_group.add_argument('--equil_ps', type=float, default=10.0,
                          help='Picoseconds for equilibration (default: 10.0)')
    sim_group.add_argument('--calcFrames', type=int, default=None,
                          help='Number of frames for calculation (default: auto)')
    sim_group.add_argument('--mbar', type=int, default=30,
                          help='Use PyMBAR to subsample frames (default: 30)')
    sim_group.add_argument('--method', type=str, choices=['gbsa', 'pbsa'], default='gbsa',
                          help='Use PBSA or GBSA (default: gbsa)')
    
    # Output options
    output_group = parser.add_argument_group('Output Options')
    output_group.add_argument('--odir', type=str, default=None,
                             help='Output directory for results (default: auto-generated)')
    output_group.add_argument('-v', '--verbose', type=int, choices=[0, 1, 2], default=1,
                             help='Verbosity level (default: 1)')
    
    args = parser.parse_args()
    
    # Validate input arguments
    if args.com is None and (args.apo is None or args.lig is None):
        parser.error("Either --com (complex file) or both --apo and --lig (separate files) must be provided")
    
    return args

def setup_folder(args):
    if args.odir is None:
        # Use a default directory name if no input file is provided
        if args.com is not None:
            args.odir = f'{os.getcwd()}/{args.com.split("/")[-1].split(".")[0]}'
        elif args.apo is not None and args.lig is not None:
            args.odir = f'{os.getcwd()}/mmgbsa_analysis'
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

    if args.amber_path is None:
        try:
            args.amber_path = os.environ['AMBERHOME']
        except KeyError:
            args.amber_path = None
            if args.v >= 1:
                print("AMBERHOME not set. Using open-source alternatives.")

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.ERROR)
    logging.getLogger('openmm').setLevel(logging.ERROR)
    logging.getLogger('openforcefield').setLevel(logging.ERROR)
    logging.getLogger('openmmtools').setLevel(logging.ERROR)
    warnings.filterwarnings("ignore")

    if args.com is None:
        if args.apo is not None and args.lig is not None:
            print("if args.com is None, args.apo and args.lig is not None")
            exit()
        else:
            sl = SystemLoader(dirpath=args.odir,
                              verbose=args.v,
                              input_pdb=args.com,  # Fixed: changed args.pdb to args.com
                              ps=args.ps,
                              calcFrames=args.calcFrames,
                              equil_ps=args.equil_ps,
                              platform_name=args.platform,
                              mbar=args.mbar
                              )
    else:
        sl = SystemLoader(dirpath=args.odir,
                          verbose=args.v,
                          input_pdb=args.com,  # Fixed: changed args.pdb to args.com
                          ps=args.ps,
                          calcFrames=args.calcFrames,
                          equil_ps=args.equil_ps,
                          platform_name=args.platform,
                          mbar=args.mbar
                          )
    sl.prepare_simulation()

    deltag, std = sl.run_amber(args.method, args.amber_path)
    print(f"{deltag} ± {std} (kcal/mol)")