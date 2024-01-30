# Examples for UMAT_Lefevre_Sozio_Lopez-Pamies

## Description

This folder contains examples showcasing the use of the [subroutine](/UMAT_KLP_RK5_hybrid.f) for:
- a specimen made out of a viscoelastic elastomer subjected to homogeneous deformations for different types of 2D, axisymmetric, and 3D Abaqus finite elements:
  - [CPE4H](Examples/CPE4H/): 4-node bilinear quadrilateral plane strain element, hybrid with constant pressure
  - [CPEG4H](Examples/CPEG4H): 4-node bilinear quadrilateral generalized plane strain element, hybrid with constant pressure
  - [CAX4H](Examples/CAX4H): 4-node bilinear quadrilateral axisymmetric element (without twist), hybrid with constant pressure
  - [C3D8H](Examples/C3D8H): 8-node linear brick element, hybrid with constant pressure
  
  Examples for a highly compressible viscoelastic elastomer ($\kappa=\mu$) and for a nearly-incompressible viscoelastic elastomer ($\kappa=10000\mu$) are included.

- a single-edge notch tension ([SENT](Examples/SENT)) specimen made out of a viscoelastic elastomer