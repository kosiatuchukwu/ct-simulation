# Interactive CT Simulation Engine

## What it does
A polychromatic CT simulation engine in MATLAB for teaching 
image quality trade-offs.

## Requirements
- MATLAB R2021a or later with Image Processing Toolbox
- Python 3.x with SpekPy (for spectrum generation only)

## How to run
1. Add spectrum CSV files using export_spectra.py
2. Run ct_simulate(struct()) from the MATLAB command window
3. For the GUI: run ct_app.m

## Files
- ct_simulate.m — core engine
- ct_app.m — graphical interface
- ct_sweep.m — parameter sweep used in Chapter 4
