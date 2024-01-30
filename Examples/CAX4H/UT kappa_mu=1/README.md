# Example for CAX4H elements: UT kappa_mu=1

This folder contains an example showcasing the use of the [subroutine](/UMAT_KLP_RK5_hybrid.f) for CAX4H 4-node bilinear quadrilateral axisymmetric element (without twist), hybrid with constant pressure.

The example corresponds to a specimen made out of a highly compressible viscoelastic elastomer ($\kappa=\mu$) subjected to uniaxial tension. 

## Usage

The example can be run from the command line as follows:
```
abaqus job=patch_CAX4H.inp user=UMAT_KLP_RK5_hybrid.f
```

A Python script `comparison.py` is included to post-process the finite-element solution and compare it to the exact solution included in the `exact_solution.csv` file. The post-processing script is run in command line as follows:
```
abaqus cae -script=comparison.py
```
The script generates and saves .png files of the stress-strain plot and the plot of the strain energies (per deformed unit volume) in the equilibrium (psiE) and non-equilibrium (psiNE) branches:

![Stress vs strain](S22%20vs%20LE22.png)

![Energy densities](psis%20vs%20t.png)