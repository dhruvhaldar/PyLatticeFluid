# PyLatticeFluid

A simple, 2D Computational Fluid Dynamics (CFD) solver based on the Lattice Boltzmann Method (LBM) in Python. 

## Features
- **D2Q9 Lattice Model:** Supports 2D flow simulations with 9 discrete velocity directions.
- **BGK Collision Operator:** Single-relaxation-time (SRT) approximation for particle collisions.
- **Boundary Conditions:**
  - On-node Bounce-Back for stationary solid obstacles.
  - Zou/He velocity boundary condition for moving walls (e.g., lid-driven cavity).
- **Vectorized Implementation:** Utilizes NumPy for efficient tensor operations.
- **Validation & Testing:** Includes a suite of Pytest unit tests for mass conservation and boundary handling, and a validation script for the classic lid-driven cavity problem.

## File Architecture
- `lbm.py`: Core `LBMSolver` class containing initialization, collision, streaming, and boundary condition logic.
- `lbm_test.py`: A validation script simulating a lid-driven cavity flow using `LBMSolver`.
- `test_lbm.py`: Pytest suite for verifying physical properties like mass conservation.
- `specifications.md`: This file, detailing the project overview.

## Requirements
- Python 3.8+
- NumPy
- Matplotlib
- Pytest

## Quick Start
1. Install dependencies: `pip install numpy matplotlib pytest`
2. Run tests: `pytest test_lbm.py`
3. Run validation simulation: `python3 lbm_test.py`
