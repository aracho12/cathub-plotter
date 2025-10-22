#!/usr/bin/env python3
"""
Setup script for cathub-plotter package
"""

from setuptools import setup, find_packages
import os

# Read the README file
def read_readme():
    with open("README.md", "r", encoding="utf-8") as fh:
        return fh.read()

# Read requirements
def read_requirements():
    with open("requirements.txt", "r", encoding="utf-8") as fh:
        return [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="cathub-plotter",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A Python package for plotting free energy diagrams from catalysis-hub.org data",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/cathub-plotter",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=read_requirements(),
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov",
            "black",
            "flake8",
            "mypy",
        ],
        "jupyter": [
            "jupyter",
            "ipywidgets",
        ],
    },
    entry_points={
        "console_scripts": [
            "cathub-plotter=cathub_plotter.cli.main:main",
        ],
    },
    include_package_data=True,
    package_data={
        "cathub_plotter": [
            "data/*.csv",
            "data/*.xlsx",
            "examples/*.yaml",
            "examples/*.txt",
        ],
    },
    keywords="catalysis, free energy, DFT, reaction mechanisms, plotting",
    project_urls={
        "Bug Reports": "https://github.com/yourusername/cathub-plotter/issues",
        "Source": "https://github.com/yourusername/cathub-plotter",
        "Documentation": "https://github.com/yourusername/cathub-plotter#readme",
    },
)
