# Examples for CPEG4H elements

This folder contains examples showcasing the use of the [subroutine](/UMAT_KLP_RK5_hybrid.f) for CPEG4H 4-node bilinear quadrilateral generalized plane strain element, hybrid with constant pressure.

All the examples correspond to a specimen made out of a viscoelastic elastomer subjected to homogeneous deformations. 

Examples for a highly compressible viscoelastic elastomer ($\kappa=\mu$) and for a nearly-incompressible viscoelastic elastomer ($\kappa=10000\mu$) are included.

## Folder naming convention

Folders containing the examples are named according to the applied loading (EBT: equi-biaxial tension, PT: planar tension, UT: uniaxial tension) and the compressibility of the material (kappa_mu=1: highly compressible, kappa_mu=10000, nearly-incompressible)

## Usage

The examples can be run from the command line (in each of the example folders) as follows:
```
abaqus job=patch_CPEG4H.inp user=UMAT_KLP_RK5_hybrid.f
```

A Python script `comparison.py` is included to post-process the finite-element solution and compare it to the exact solution included in the `exact_solution.csv` file. The post-processing script is run in command line (where the `patch_CPEG4H.odb` file is located) as follows:
```
abaqus cae -script=comparison.py
```
The script generates and saves .png files of the stress-strain plot and the plot of the strain energies (per deformed unit volume) in the equilibrium (psiE) and non-equilibrium (psiNE) branches. 