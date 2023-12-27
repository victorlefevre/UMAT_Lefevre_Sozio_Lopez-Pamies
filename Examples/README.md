# Examples for UMAT_Lefevre_Sozio_Lopez-Pamies

## Description

This folder contains examples showcasing the use of the [subroutine](/UMAT_KLP_RK5_hybrid.f) for different types of 2D, axisymmetric, and 3D Abaqus finite elements:
- [CPE4H](CPE4H): 4-node bilinear quadrilateral plane strain element, hybrid with constant pressure
- [CPEG4H](CPEG4H): 4-node bilinear quadrilateral generalized plane strain element, hybrid with constant pressure
- [CAX4H](CAX4H): 4-node bilinear quadrilateral axisymmetric element (without twist), hybrid with constant pressure
- [C3D8H](C3D8H): 8-node linear brick element, hybrid with constant pressure

All the examples correspond to a specimen made out of a viscoelastic elastomer subjected to homogeneous deformations. Examples for a highly compressible viscoelastic elastomer ($\kappa=\mu$) and for a nearly-incompressible viscoelastic elastomer ($\kappa=10000\mu$) are included.