# AMS 595 - Project 2: Coastline of the Mandelbrot Set

This repository contains the MATLAB implementation for approximating the circumference of the Mandelbrot set's upper boundary, exploring the fractal nature of coastlines as described in Mandelbrot's famous paper "How Long Is the Coast of Britain?"

## Project Overview

The project implements a numerical method to:
- Compute divergence iterations for complex points in the Mandelbrot set
- Locate the set boundary using bisection search along vertical lines
- Fit a polynomial to sampled boundary points
- Calculate the curve length through numerical integration

## Files

- `ams595_hw2.m` - Main MATLAB script containing all functions and execution code
- `AMS595_HW2.pdf` - Complete project report with analysis and results

## Implementation

The single file `ams595_hw2.m` contains:

### Main Script
- Samples the Mandelbrot boundary at multiple x-points
- Performs polynomial fitting (order 15) to boundary data
- Computes the integrated length of the fitted curve
- Generates visualization plots

## Key Results

- **Estimated Boundary Length**: Approximately 2.88 units
- **Polynomial Fit**: 15th order polynomial
- **Sampling Resolution**: 1001 points across x ∈ [-2, 1]
- **Valid Boundary Points**: 714 samples used for fitting

## Usage

### MATLAB
```matlab
% Run the main script
ams595_hw2
