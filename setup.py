#!/usr/bin/env python3
"""
Setup script for Open-Source MMGBSA Analysis Package
"""

from setuptools import setup, find_packages
import os

# Read the README file
def read_readme():
    readme_path = os.path.join(os.path.dirname(__file__), 'readme.md')
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            return f.read()
    return "Open-Source MMGBSA Analysis Package"

# Read requirements
def read_requirements():
    requirements_path = os.path.join(os.path.dirname(__file__), 'requirements.txt')
    if os.path.exists(requirements_path):
        with open(requirements_path, 'r', encoding='utf-8') as f:
            return [line.strip() for line in f if line.strip() and not line.startswith('#')]
    return []

setup(
    name="mmgpbsa-open",
    version="1.0.0",
    author="Open Source MMGBSA Team",
    author_email="contact@example.com",
    description="Open-source MMGBSA analysis package",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/mmgpbsa-open",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    python_requires=">=3.7",
    install_requires=read_requirements(),
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.8",
            "mypy>=0.800",
        ],
        "docs": [
            "sphinx>=3.0",
            "sphinx-rtd-theme>=0.5",
        ],
    },
    entry_points={
        "console_scripts": [
            "mmgpbsa=run_mmgpbsa:main",
            "component-scoring=component_scoring:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
    keywords="molecular dynamics, mmgbsa, free energy, drug discovery, computational chemistry",
    project_urls={
        "Bug Reports": "https://github.com/yourusername/mmgpbsa-open/issues",
        "Source": "https://github.com/yourusername/mmgpbsa-open",
        "Documentation": "https://mmgpbsa-open.readthedocs.io/",
    },
) 