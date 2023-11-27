!**********************************************************************
! Legal notice: UMAT_KLP_RK5_hybrid.f 
!
! Copyright (C) 2023 Victor Lefèvre
!                    Fabio Sozio (fsozio@illinois.edu)
!                    Oscar Lopez-Pamies (pamies@illinois.edu)
!
! This Abaqus UMAT subroutine implements a family of internal-variable-
! based constitutive models for the viscoelastic response of elastomers
! of any compressibility — including fully incompressible elastomers — 
! undergoing finite deformations, that introduced by Kumar and 
! Lopez-Pamies [1]. The models can be chosen to account for a wide 
! range of non-Gaussian elasticities, as well as for a wide range of 
! nonlinear viscosities. 
!
! The present subroutine is implemented for the choice of two-term 
! Lopez-Pamies non-Gaussian elasticities [2] and deformation-enhanced 
! shear-thinning nonlinear viscosity proposed in [1]. It relies on a 
! hybrid formulation [3] to deal with nearly or fully incompressible 
! elastomers.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see https://www.gnu.org/licenses/
!
!**********************************************************************
! Usage:
!
! The subroutine is to be used as an USER MATERIAL with 15 material 
! CONSTANTS and a TOTAL HYBRID FORMULATION, e.g., 
! *USER MATERIAL,CONSTANTS=15,HYBRID FORMULATION=TOTAL
! in the input (.inp) file. Due to the hybrid nature of the underlying
! formulation, hybrid finite elements must be used in the simulation.
!
! 1/ The 15 materials properties required by the model to be input to 
! the subroutine via the PROPS array are listed in order below:
!
!  AMU1    = PROPS(1)  ! PARAMETER #1 OF THE EQUILIBRIUM ELASTICITY
!  AALPHA1 = PROPS(2)  ! EXPONENT #1 OF THE EQUILIBRIUM ELASTICITY
!  AMU2    = PROPS(3)  ! PARAMETER #2 OF THE EQUILIBRIUM ELASTICITY
!  AALPHA2 = PROPS(4)  ! EXPONENT #2 OF THE EQUILIBRIUM ELASTICITY
!  AKAPPAE = PROPS(5)  ! BULK MODULUS OF THE EQUILIBRIUM ELASTICITY 
!  AM1     = PROPS(6)  ! PARAMETER #1 OF THE NON-EQUILIBRIUM ELASTICITY
!  AA1     = PROPS(7)  ! EXPONENT #1 OF THE NON-EQUILIBRIUM ELASTICITY
!  AM2     = PROPS(8)  ! PARAMETER #2 OF THE NON-EQUILIBRIUM ELASTICITY
!  AA2     = PROPS(9)  ! EXPONENT #2 OF THE NON-EQUILIBRIUM ELASTICITY
!  AETA0   = PROPS(10) ! PARAMETER #1 OF THE DEVIATORIC VISCOSITY
!  AETAI   = PROPS(11) ! PARAMETER #2 OF THE DEVIATORIC VISCOSITY
!  ABETA1  = PROPS(12) ! EXPONENT #1 OF THE DEVIATORIC VISCOSITY
!  ABETA2  = PROPS(13) ! EXPONENT #2 OF THE DEVIATORIC VISCOSITY
!  AK1     = PROPS(14) ! FACTOR #1 OF THE DEVIATORIC VISCOSITY
!  AK2     = PROPS(15) ! FACTOR #2 OF THE DEVIATORIC VISCOSITY
!
! The material parameters AMU1, AMU2 characterizing the equilibrium     
! elasticity are non-negative real numbers 
! (AMU1 >= 0, AMU2 >= 0) with strictly positive sum (AMU1 + AMU2 > 0). 
! The two exponents AALPHA1, AALPHA2 are non-zero real numbers 
! (AALPHA1 ≠ 0, AALPHA2 ≠ 0) with MAX{AALPHA1, AALPHA2} > 0 leading to  
! a strongly elliptic strain energy (see eq. (22) in [2]). This is left 
! to the user to check.
!
! The material parameter AKAPPAE charcterizing the compressibility of
! the equilibrium elasticity is a strictly positive real number
! (AKAPPAE > 0).
!
! The material parameters AM1, AM2 characterizing the non-equilibrium     
! elasticity are non-negative real numbers 
! (AM1 >= 0, AM2 >= 0) with strictly positive sum (AM1 + AM2 > 0). 
! The two exponents AA1, AA2 are non-zero real numbers 
! (AA1 ≠ 0, AA2 ≠ 0) with MAX{AA1, AA2} > 0 leading to a 
! strongly elliptic strain energy (see eq. (22) in [2]). This is left 
! to the user to check.
!
! The material parameters AETA0, AETAI, ABETA1, ABETA2, AK1, AK2
! characterizing the viscous dissipation are all non-negative real  
! numbers (AETA0 >= 0, AETAI >= 0, ABETA1 >= 0, ABETA2 >= 0, AK1 >= 0,
! AK2 >= 0).
!
! 2/ The subroutine uses solution-dependent state variables for the
! independent components of the symmetric tensorial internal variable 
! Cvi. The space for the state variables must be allocated in the input 
! file using the *DEPVAR keyword:
! *DEPVAR
! NSTATV
! where NSTATV is the number of state variables required. NSTATV depends 
! on the type of finite elements used:
! - NSTATV = 4 for elements with only 1 non-zero shear component
!   (CPE3H, CPE4H, CPE6H, CPE8H, CAX3H, CAX4H, CAX6H, CAX8H, CPEG3H, 
!    CPEG4H, CPEG6H, CPEG8H, ...)
! - NSTATV = 6 for elements with all 3 non-zero shear components 
!   (C3D4H, C3D8H, C3D10H, C3D20H, CGAX3H, CGAX4H, CGAX6H, CGAX8H, ...)
! Note that an incorrect value of NSTATV will trigger an error in the 
! subroutine and the error message in the .msg file will indicate the 
! expected value for NSTATV.
! 
! 3/ The values of state variables must be initialized according to the 
! initial condition Cvi=Id using the *INITIAL CONDITIONS keyword as
! follows:
! - NSTATV = 4:
! *INITIAL CONDITIONS, TYPE=SOLUTION
! ALLELEMENTSET,1.0,1.0,1.0,0.0
! - NSTATV = 6:
! *INITIAL CONDITIONS, TYPE=SOLUTION
! ALLELEMENTSET,1.0,1.0,1.0,0.0,0.0,0.0
! where ALLELEMENTSET is an element set of all finite elements 
! associated with the UMAT material.
!
! 4/ The subroutine can create at the user's request a user output 
! variable UVARM1 for the determinant of the deformation gradient. The 
! space for the user output variable must be allocated in the input 
! file using the *USER OUTPUT VARIABLES keyword:
! *USER OUTPUT VARIABLES
! 1
! The user must also include UVARM in the field output request.
!**********************************************************************
! Additional information:
!
! This subroutine doesn't create predefined field variables. 
!
! Please consult the ABAQUS Documentation for additional references 
! regarding the use of the UMAT subroutine for (near-)incompressible
! material models.
!
!**********************************************************************
! References:
!
! [1] Kumar, A., Lopez-Pamies, O., 2016. On the two-potential constitu-
!     tive modeling of rubber viscoelastic materials, Comptes Rendus
!     Mecanique 344, 102–-112.
! [2] Lopez-Pamies, O., 2010. A new I1-based hyperelastic model for 
!     rubber elastic materials. C. R. Mec. 338, 3--11.
! [3] Lefèvre, V., Sozio, F., Lopez-Pamies, O., 2023. Abaqus implemen-
!     tation of a large family of finite viscoelasticity models.
!     Submitted.
!
!**********************************************************************