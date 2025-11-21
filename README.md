# Transieata_calculator

Overview:
This repository automates the complex and often tedious workflow of Quantum Transport Calculations using Siesta and TranSIESTA (NEGF + DFT). 

It was designed specifically to address the stability issues often encountered in 1D/2D/3D systems in transport calculations (e.g., Density Matrix divergence or "explosions"). It handles the entire pipeline—from initial convergence testing to the final I-V characteristics—with zero manual file copying. 

The Problem:
1.  Divergence: Transport calculations often fail with `dDmax` values exploding (e.g., > 40.0) instead of converging.
2.  Manual Labor:  Users waste hours manually testing Mesh/K-points and copying `.TSHS` and `.DM` files between directories.
3.  Unit Mismatches: Inconsistent K-point blocks or Energy units often crash runs.

The Solution:
This script uses a final optimized propagation logic. It optimizes parameters sequentially and automatically injects "High Stability " settings (Safe Mixing and Complex Contour Integration) into the device calculation to ensure convergence.

Key Features

1. Sequential Parameter Optimization
The script performs convergence tests in a strict order, using the optimized result of step $N$ as the input for step $N+1$:
Mesh Cutoff: Optimizes the grid fineness.
K-Points: Optimizes Brillouin zone sampling (supports 1D, 2D, and 3D grids).
Lattice Constant: Optimizes the electrode Z-vector (critical for transport).

2. Device "Anti-Explosion" Stabilization
To prevent the density matrix from diverging, the script automatically injects robust parameters into the `device.fdf`:
Safe Mixing: Forces `SCF.Mixer.Method Pulay`, a low weight (`0.01`), and a high history (`50` steps) to dampen oscillations.
Complex Contour Integration: Injects `Eq` and `NonEq` contour circles (C-Left/Right, T-Left/Right) essential for convergence in open systems.
Extended Iterations: Automatically sets `MaxSCFIterations` to 5000 to allow slow-converging systems to settle.

3. Smart Voltage Sweeps
Chained Restarts: Uses the converged Density Matrix from $0.0V$ as the initial guess for $0.1V$, and so on.
Auto-Plotting: Automatically generates a Transmission Spectrum (`.png`) after every single calculation using `sisl`and tbtrans. so you can monitor results in real-time.

Installation & Prerequisites

1. System Requirements
Siesta / TranSIESTA: Must be installed and accessible (e.g., via `mpirun`).
Python 
Dependencies: numpy, matplotlib, sisl

2. Files and Directories:

electrode.fdf (electrode structure fdf file with initial parameters, Note include here K points, MeshCutoff, Z-vector for the initial guess)
device.fdf (just need to change the lattice coordinates in the fdf and other parameters if needed, in this case z vector is changed)
Au.psf (pseudo potential file for Au).
calc.py (python file). 


## What could be the user input once you start the calculation, an example:

Symmetrical Electrodes? y 

MPI Processors: 24

Convergence Threshold: 0.01 eV

Mesh Sequence: 150 200

K-Point Sequence: 50 75 100 

Lattice Sequence: 2.80 2.85

K-Point Dimension: 1D

Fixed Transverse K-points (kx​ ,ky): 1 and 1

Voltage List: 0.01 0.05


