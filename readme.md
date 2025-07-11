# Open-Source MMGBSA Analysis Package

A fully open-source Python package for performing MMGBSA (Molecular Mechanics Generalized Born Surface Area) analysis for drug discovery and computational chemistry applications.

## Features

- **Fully Open Source**: No proprietary software dependencies
- **MMGBSA Calculations**: Binding free energy calculations using Generalized Born models
- **Molecular Docking**: Integration with Vina for molecular docking
- **Molecular Properties**: Comprehensive molecular property calculations
- **Trajectory Analysis**: MDTraj-based trajectory processing
- **Statistical Analysis**: PyMBAR integration for advanced statistical analysis

## Open Source Alternatives

This package replaces proprietary software with open-source alternatives:

| Original (Proprietary) | Open Source Alternative | Purpose |
|------------------------|------------------------|---------|
| OpenEye Toolkits | RDKit + Open Babel | Molecular file I/O, charge assignment |
| AMBER | OpenMM | Molecular dynamics, force fields |
| OpenEye Docking | Vina/Smina | Molecular docking |
| AMBER MMGBSA | OpenMM GB | MMGBSA calculations |

## Installation

### Prerequisites

1. **Python 3.7+**
2. **Open Babel** (for molecular format conversion)
3. **Vina** (for molecular docking, optional)

### Install Open Babel

**Ubuntu/Debian:**
```bash
sudo apt-get install openbabel
```

**macOS:**
```bash
brew install open-babel
```

**Windows:**
Download from [Open Babel website](http://openbabel.org/wiki/Get_Open_Babel)

### Install Vina (Optional)

**Ubuntu/Debian:**
```bash
sudo apt-get install autodock-vina
```

**macOS:**
```bash
brew install vina
```

### Install Python Package

```bash
# Clone the repository
git clone https://github.com/yourusername/mmgpbsa-open.git
cd mmgpbsa-open

# Install the package
pip install -e .

# Or install from requirements
pip install -r requirements.txt
```

## Quick Start

### Basic MMGBSA Analysis

```bash
# Run MMGBSA analysis on a protein-ligand complex
python run_mmgpbsa.py --com complex.pdb --odir results --platform CPU
```

### Component Scoring

```bash
# Score a protein-ligand complex
python component_scoring.py --complex complex.pdb
```

### Advanced Usage

```bash
# Run with specific parameters
python run_mmgpbsa.py \
    --com complex.pdb \
    --odir results \
    --platform CUDA \
    --ps 1000 \
    --mbar 50 \
    --method gbsa \
    -v 2
```

## Usage Examples

### MMGBSA Analysis

```python
from mmgpbsa.systemloader import SystemLoader

# Initialize system
sl = SystemLoader(
    dirpath="results",
    verbose=1,
    input_pdb="complex.pdb",
    ps=1000,
    calcFrames=50,
    platform_name="CPU"
)

# Prepare and run simulation
sl.prepare_simulation()

# Calculate MMGBSA
delta_g, std = sl.run_amber("gbsa", amber_path=None)
print(f"ΔG = {delta_g:.2f} ± {std:.2f} kcal/mol")
```

### Molecular Docking

```python
from component_scoring import score_and_build

# Score a complex
scores = score_and_build("complex.pdb")
print(scores)
```

## Configuration

### Force Fields

The package supports multiple force fields through OpenMM:
- AMBER14 (default)
- CHARMM36
- OPLS-AA

### GB Models

Available Generalized Born models:
- GBn (igb=5)
- GBn2 (igb=7)
- OBC2 (igb=8)

### Platforms

Supported OpenMM platforms:
- CPU (default)
- CUDA (NVIDIA GPUs)
- OpenCL (AMD GPUs)

## Output Files

The analysis generates several output files:

- `trajectory.dcd`: Molecular dynamics trajectory
- `com.pdb`: Complex structure
- `apo.pdb`: Protein structure
- `lig.pdb`: Ligand structure
- `results.log`: Analysis log file

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

### Development Setup

```bash
# Clone repository
git clone https://github.com/yourusername/mmgpbsa-open.git
cd mmgpbsa-open

# Install in development mode
pip install -e .[dev]

# Run tests
pytest

# Format code
black .

# Type checking
mypy mmgpbsa/
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this package in your research, please cite:

```bibtex
@software{mmgpbsa_open,
  title={Open-Source MMGBSA Analysis Package},
  author={Open Source MMGBSA Team},
  year={2024},
  url={https://github.com/yourusername/mmgpbsa-open}
}
```

## Acknowledgments

This package builds upon several excellent open-source projects:
- [OpenMM](http://openmm.org/) - Molecular dynamics
- [RDKit](https://www.rdkit.org/) - Molecular informatics
- [MDTraj](http://mdtraj.org/) - Trajectory analysis
- [PyMBAR](https://github.com/choderalab/pymbar) - Statistical analysis
- [Vina](http://vina.scripps.edu/) - Molecular docking

## Support

- **Documentation**: [https://mmgpbsa-open.readthedocs.io/](https://mmgpbsa-open.readthedocs.io/)
- **Issues**: [GitHub Issues](https://github.com/yourusername/mmgpbsa-open/issues)
- **Discussions**: [GitHub Discussions](https://github.com/yourusername/mmgpbsa-open/discussions)

