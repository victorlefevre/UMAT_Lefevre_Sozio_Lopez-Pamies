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
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS),
     3 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
C ---------------------------------------------------------------------
C
C    IMPORTANT ABAQUS VARIABLES
C
C    DFGRD0 - DEFORMATION GRADIENT AT THE BEGINNING OF THE 
C             INCREMENT
C    DFGRD1 - DEFORMATION GRADIENT AT THE END OF THE INCREMENT
C    DDSDDE - JACOBIAN MATRIX
C    STRESS - CAUCHY STRESS AT THE BEGINNING OF THE INCREMENT.
C             MUST BE UPDATED TO CAUCHY STRESS AT THE END OF 
C             THE INCREMENT
C    STATEV - STATE VARIABLES
C             INDEPENDENT COMPONENTS OF CVI, THE INVERSE OF 
C             THE TENSORIAL INTERNAL VARIABLE CV AT THE BEGINNING 
C             OF THE INCREMENT. 
C             MUST BE UPDATED TO CVI AT THE END OF THE INCREMENT. 
C             COMPONENTS STORED AS CVI_11, CVI_22, CVI_33, CVI_12
C             WHEN NSTATV=4, AND AS CVI_11, CVI_22, CVI_33, CVI_12, 
C             CVI_13, CVI_23 WHEN NSTATV=6.
C ---------------------------------------------------------------------
C
C    IMPORTANT LOCAL VARIABLES:
C
C    KINEMATICS:
C       BBARTENS  - DEVIATORIC LEFT CAUCHY-GREEN TENSOR STORED AS 
C                   NTENS VECTOR 
C       BEBARTENS - DEVIATORIC SYMMETRIC SECOND ORDER TENSOR DEFINED
C                   AS FBAR.CVI.FBAR^T AND STORED AS NTENS VECTOR
C       CVI       - INVERSE OF THE TENSORIAL INTERNAL VARIABLE CV
C       CVIDOT    - 1ST ORDER TIME DERIVATIVE OF CVI
C       DETF      - DETERMINANT OF F AT THE END OF THE INCREMENT 
c                   DETF=DET(DFGRD1)
C
C    MATERIAL COEFFICIENTS
C       AMU1    = PROPS(1)  
C       AALPHA1 = PROPS(2)  
C       AMU2    = PROPS(3)  
C       AALPHA2 = PROPS(4) 
C       AKAPPAE = PROPS(5)   
C       AM1     = PROPS(6)  
C       AA1     = PROPS(7)  
C       AM2     = PROPS(8)  
C       AA2     = PROPS(9)  
C       AETA0   = PROPS(10)  
C       AETAI   = PROPS(11)  
C       ABETA1  = PROPS(12)  
C       ABETA2  = PROPS(13)  
C       AK1     = PROPS(14)  
C       AK2     = PROPS(15)  
C
C    RUNGE-KUTTA INTEGRATION SCHEME:
C       NRK     - NUMBER OF STAGES FOR THE RUNGE-KUTTA INTEGRATION
C                 SCHEME
C       RKAIJS  - A_IJ COEFFICIENTS OF THE RUNGE-KUTTA INTEGRATION 
C                 SCHEME STORED AS NRKXNRK ARRAY
C       RKBIS   - B_I COEFFICIENTS OF THE RUNGE-KUTTA INTEGRATION 
C                 SCHEME STORED AS NRK ARRAY
C       RKCIS   - C_I COEFFICIENTS OF THE RUNGE-KUTTA INTEGRATION 
C                 SCHEME STORED AS NRK ARRAY
C       RKCVI   - RUNGE-KUTTA INTERMEDIATE EVALUATION OF CVI STORED AS 
C                 NTENS VECTOR
C       RKCBARTENS - RUNGE-KUTTA INTERMEDIATE EVALUATION OF THE  
C                    DEVIATORIC RIGHT CAUCHY-GREEN TENSOR STORED AS 
C                    NTENS VECTOR
C       RKKIJS  - RUNGE-KUTTA INTERMEDIATE EVALUATIONS OF THE RHS OF 
C                 THE EVOLUTION EQUATION (I.E., CVIDOT) (K^(1)_IJ,
C                 K^(2)_IJ,...,K^(NRK)_IJ) STORED AS NTENSxNRK ARRAY    
C    
C ---------------------------------------------------------------------
C
      DIMENSION RKCBARTENS(NTENS), RKCVI(NTENS), CVIDOT(NTENS),
     1 BBARTENS(NTENS), BEBARTENS(NTENS), DEVBBAR(NTENS), 
     2 DEVBEBAR(NTENS), CVI0(NTENS), CVITEMP(NTENS)

      DOUBLE PRECISION, ALLOCATABLE :: RKKS(:,:), RKAIJS(:,:), RKBIS(:), 
     1 RKCIS(:)
C
      NRK=6
C            
      ALLOCATE(RKKS(NTENS,NRK), RKAIJS(NRK,NRK), RKBIS(NRK), RKCIS(NRK))
C
      IF ((KSTEP.EQ.1).AND.(KINC.EQ.1)) THEN
        CALL KCHECKS(STATEV, NSTATV, NTENS, PROPS, NPROPS)
      END IF
C ---------------------------------------------------------------------
C
C INTEGRATE EVOLUTION EQUATION FOR CVI
C 
C
C FIFTH ORDER RUNGE-KUTTA SCHEME COEFFICIENTS
C     
      CALL KINITIA(RKAIJS,NRK*NRK)
      RKAIJS(2,1)=0.5
      RKAIJS(3,1)=3.0/16.0
      RKAIJS(3,2)=1.0/16.0
      RKAIJS(4,3)=0.5
      RKAIJS(5,2)=-3.0/16.0
      RKAIJS(5,3)=6.0/16.0
      RKAIJS(5,4)=9.0/16.0
      RKAIJS(6,1)=1.0/7.0
      RKAIJS(6,2)=4.0/7.0
      RKAIJS(6,3)=6.0/7.0
      RKAIJS(6,4)=-12.0/7.0
      RKAIJS(6,5)=8.0/7.0
C      
      RKBIS(1)=7.0/90.0
      RKBIS(2)=0.0/90.0
      RKBIS(3)=32.0/90.0
      RKBIS(4)=12.0/90.0
      RKBIS(5)=32.0/90.0
      RKBIS(6)=7.0/90.0
C      
      RKCIS(1)=0.0
      RKCIS(2)=0.5
      RKCIS(3)=0.25 
      RKCIS(4)=0.5 
      RKCIS(5)=0.75 
      RKCIS(6)=1.0  
C
C RUNGE-KUTTA INTEGRATION STEPS
C      
      CALL KINITIA(RKKS,NTENS*NRK)
      DO KRKSTEP=1,NRK   
C       COMPUTE RKCBARTENS
        RKC=RKCIS(KRKSTEP) 
        CALL KRKCBARTENS(RKCBARTENS,DFGRD0,DFGRD1,RKC,NSHR,NTENS)
C       COMPUTE RKCVI
        DO I = 1,NTENS
          RKCVI(I) = STATEV(I)
          IF (KRKSTEP.GT.1) THEN
            DO J=1,KRKSTEP-1   
              RKCVI(I) = RKCVI(I)+DTIME*RKAIJS(KRKSTEP,J)*RKKS(I,J)
            END DO
          END IF          
        END DO   
C       COMPUTE CVIDOT
        CALL KCVIDOT(CVIDOT,RKCBARTENS,RKCVI,PROPS,NPROPS,NDI,
     1        NSHR,NTENS)
C       STORE RESULT IN COLUMN #KRKSTEP OF RKKS ARRAY
        CALL KINSERTCOLUMN(CVIDOT,RKKS,KRKSTEP,NTENS,NRK)
      END DO
C
C COMPUTE RK TIME INTEGRATION
C      
      DO I = 1,NTENS
        CVI0(I) = STATEV(I)
        CVITEMP(I) = STATEV(I)
        DO KRKSTEP = 1, NRK 
          CVITEMP(I) = CVITEMP(I) + 
     1                 DTIME*RKBIS(KRKSTEP)*RKKS(I,KRKSTEP)
        END DO
      END DO
C
C CORRECTOR STEP FOR CONSTRAINT DET(CVI)=1 AND UPDATE INTERNAL VARIABLE  
C FOR CVI AT END OF INCREMENT 
C   
      CALL KRKCORRECTOR(STATEV,CVI0,CVITEMP,NSHR,NTENS)
C ---------------------------------------------------------------------
C
C COMPUTE STRESS AND JACOBIAN DDSDDE AT END OF INCREMENT
C 
C
C COMPUTE DETF AT END OF INCREMENT
C
      CALL KDET(DETF,DFGRD1,NSHR)
      AJHAT=STRESS(NTENS+1)
C
C COMPUTE BBARTENS, BEBARTENS AT END OF INCREMENT
C
      CALL KBBAR(BBARTENS, DFGRD1, DETF, NSHR, NTENS)
      CALL KBEBAR(BEBARTENS, DFGRD1, STATEV, DETF, NSHR, NTENS)
C
C COMPUTE DEVBBAR, DEVBEBAR AT END OF INCREMENT
C
      CALL KTRACEVECTOR(AI1B, BBARTENS, NDI, NTENS)
      CALL KBDEVA(DEVBBAR, BBARTENS, AI1B, NSHR, NTENS) 
      CALL KTRACEVECTOR(AI1BE, BEBARTENS, NDI, NTENS) 
      CALL KBDEVA(DEVBEBAR, BEBARTENS, AI1BE, NSHR, NTENS)      
C
C MATERIAL PARAMETERS
C   
      AMU1    = PROPS(1)  
      AALPHA1 = PROPS(2)  
      AMU2    = PROPS(3)  
      AALPHA2 = PROPS(4) 
      AKAPPAE = PROPS(5)      
      AM1     = PROPS(6)  
      AA1     = PROPS(7)  
      AM2     = PROPS(8)  
      AA2     = PROPS(9)  
C
C 1ST and 2ND DERIVATIVES OF EQUILIBRIUM ELASTICITY
C                                 
      DUEDI1B = (3**(1-AALPHA1)*AI1B**(AALPHA1-1)*AMU1)/2. +
     1          (3**(1-AALPHA2)*AI1B**(AALPHA2-1)*AMU2)/2.
      D2UEDI1B2 = (3**(1-AALPHA1)*(AALPHA1-1)*AI1B**(AALPHA1-2)*AMU1)/2.+
     1            (3**(1-AALPHA2)*(AALPHA2-1)*AI1B**(AALPHA2-2)*AMU2)/2.
      DUEDJHAT = (AJHAT-1.0)*AKAPPAE
      D2UEDJHAT2 = AKAPPAE
C
C 1ST and 2ND DERIVATIVES OF NON-EQUILIBRIUM ELASTICITY
C	        
      DUNEDI1BE = (3**(1-AA1)*AI1BE**(AA1-1)*AM1)/2. + 
     1            (3**(1-AA2)*AI1BE**(AA2-1)*AM2)/2.
      D2UNEDI1BE2 = (3**(1-AA1)*(AA1-1)*AI1BE**(AA1-2)*AM1)/2. +
     1              (3**(1-AA2)*(AA2-1)*AI1BE**(AA2-2)*AM2)/2.
C     
C SIMPLIFYING NOTATIONS
C       
      BB1 = DUEDI1B/DETF
      BB2 = 2.0*BB1
      BB3 = BB2/3.0
      BB4 = 2.0*BB3
C      
      BB1E = DUNEDI1BE/DETF
      BB2E = 2.0*BB1E
      BB3E = BB2E/3.0
      BB4E = 2.0*BB3E  
C      
      CC1  = 4.0/DETF*D2UEDI1B2
      CC1E = 4.0/DETF*D2UNEDI1BE2
C
      DD1 = AI1B/3.0
      DD1E = AI1BE/3.0
C	  	  
      PHAT = DUEDJHAT      
      AKHAT = DETF*D2UEDJHAT2
C
C COMPUTE STRESS
C      
      DO K1=1,NDI
        STRESS(K1)=BB2*DEVBBAR(K1)+BB2E*DEVBEBAR(K1)+
     1             PHAT
      END DO
      DO K1=NDI+1,NDI+NSHR
        STRESS(K1)=BB2*BBARTENS(K1)+BB2E*BEBARTENS(K1)
      END DO
      STRESS(NTENS+2)=AKHAT
      STRESS(NTENS+3)=0.0 
C
C COMPUTE JACOBIAN DDSDDE
C      
      DDSDDE(1, 1)=CC1*DEVBBAR(1)**2.0+
     1             BB4*(BBARTENS(1)+DD1)+
     2             CC1E*DEVBEBAR(1)**2.0+
     3             BB4E*(BEBARTENS(1)+DD1E)+AKHAT
      DDSDDE(2, 2)=CC1*DEVBBAR(2)**2.0+
     1             BB4*(BBARTENS(2)+DD1)+
     2             CC1E*DEVBEBAR(2)**2.0+
     3             BB4E*(BEBARTENS(2)+DD1E)+AKHAT
      DDSDDE(3, 3)=CC1*DEVBBAR(3)**2.0+
     1             BB4*(BBARTENS(3)+DD1)+
     2             CC1E*DEVBEBAR(3)**2.0+
     3             BB4E*(BEBARTENS(3)+DD1E)+AKHAT
C     
      DDSDDE(1, 2)=CC1*DEVBBAR(1)*DEVBBAR(2)-
     1             BB4*(DEVBBAR(1)+BBARTENS(2))+
     2             CC1E*DEVBEBAR(1)*DEVBEBAR(2)-
     3             BB4E*(DEVBEBAR(1)+BEBARTENS(2))+AKHAT
      DDSDDE(1, 3)=CC1*DEVBBAR(1)*DEVBBAR(3)-
     1             BB4*(DEVBBAR(1)+BBARTENS(3))+
     2             CC1E*DEVBEBAR(1)*DEVBEBAR(3)-
     3             BB4E*(DEVBEBAR(1)+BEBARTENS(3))+AKHAT
      DDSDDE(2, 3)=CC1*DEVBBAR(2)*DEVBBAR(3)-
     1             BB4*(DEVBBAR(2)+BBARTENS(3))+
     2             CC1E*DEVBEBAR(2)*DEVBEBAR(3)-
     3             BB4E*(DEVBEBAR(2)+BEBARTENS(3))+AKHAT
C      
      DDSDDE(1, 4)=(CC1*DEVBBAR(1)+BB3)*BBARTENS(4)+
     1             (CC1E*DEVBEBAR(1)+BB3E)*BEBARTENS(4) 
      DDSDDE(2, 4)=(CC1*DEVBBAR(2)+BB3)*BBARTENS(4)+
     1             (CC1E*DEVBEBAR(2)+BB3E)*BEBARTENS(4) 
      DDSDDE(3, 4)=(CC1*DEVBBAR(3)-BB4)*BBARTENS(4)+
     1             (CC1E*DEVBEBAR(3)-BB4E)*BEBARTENS(4) 
C     
      DDSDDE(4, 4)=CC1*BBARTENS(4)**2.0+
     1             BB1*(BBARTENS(1)+BBARTENS(2))+
     2             CC1E*BEBARTENS(4)**2.0+
     3             BB1E*(BBARTENS(1)+BBARTENS(2))
C     
      IF(NSHR.EQ.3) THEN
        DDSDDE(1, 5)=(CC1*DEVBBAR(1)+BB3)*BBARTENS(5)+
     1             (CC1E*DEVBEBAR(1)+BB3E)*BEBARTENS(5) 
        DDSDDE(2, 5)=(CC1*DEVBBAR(2)-BB4)*BBARTENS(5)+
     1             (CC1E*DEVBEBAR(2)-BB4E)*BEBARTENS(5) 
        DDSDDE(3, 5)=(CC1*DEVBBAR(3)+BB3)*BBARTENS(5)+
     1             (CC1E*DEVBEBAR(3)+BB3E)*BEBARTENS(5) 
C     
        DDSDDE(1, 6)=(CC1*DEVBBAR(1)-BB4)*BBARTENS(6)+
     1             (CC1E*DEVBEBAR(1)-BB4E)*BEBARTENS(6) 
        DDSDDE(2, 6)=(CC1*DEVBBAR(2)+BB3)*BBARTENS(6)+
     1             (CC1E*DEVBEBAR(2)+BB3E)*BEBARTENS(6) 
        DDSDDE(3, 6)=(CC1*DEVBBAR(3)+BB3)*BBARTENS(6)+
     1             (CC1E*DEVBEBAR(3)+BB3E)*BEBARTENS(6) 
C        
        DDSDDE(5, 5)=CC1*BBARTENS(5)**2.0+
     1             BB1*(BBARTENS(1)+BBARTENS(3))+
     2             CC1E*BEBARTENS(5)**2.0+
     3             BB1E*(BBARTENS(1)+BBARTENS(3))
        DDSDDE(6, 6)=CC1*BBARTENS(6)**2.0+
     1             BB1*(BBARTENS(2)+BBARTENS(3))+
     2             CC1E*BEBARTENS(6)**2.0+
     3             BB1E*(BBARTENS(2)+BBARTENS(3))
C        
        DDSDDE(4,5)=CC1*BBARTENS(4)*BBARTENS(5)+
     1              BB1*BBARTENS(6)+
     2              CC1E*BEBARTENS(4)*BEBARTENS(5)+
     3              BB1E*BEBARTENS(6)
        DDSDDE(4,6)=CC1*BBARTENS(4)*BBARTENS(6)+
     1              BB1*BBARTENS(5)+
     2              CC1E*BEBARTENS(4)*BEBARTENS(6)+
     3              BB1E*BEBARTENS(5)
        DDSDDE(5,6)=CC1*BBARTENS(5)*BBARTENS(6)+
     1              BB1*BBARTENS(4)+
     2              CC1E*BEBARTENS(5)*BEBARTENS(6)+
     3              BB1E*BEBARTENS(4)
      END IF
      DO K1=1, NTENS
        DO K2=1, K1-1
          DDSDDE(K1, K2)=DDSDDE(K2, K1)
        END DO
      END DO
C	  
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE KCVIDOT(CVIDOT, CTENS, CVI, PROPS, NPROPS, NDI, 
     1 NSHR, NTENS)
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION CVIDOT(NTENS), CTENS(NTENS), CVI(NTENS), PROPS(NPROPS),
     1 CVICCVI(NTENS) 
C     
C    RETURNS CVIDOT, 1ST ORDER TIME DERIVATIVE OF CVI
C
C    IMPORTANT LOCAL VARIABLES:
C       AETAK     - DEVIATORIC VISCOSITY
C       CVI       - INVERSE OF TENSORIAL INTERNAL VARIABLE CV
C       CVICCVI   - SYMMETRIC SECOND ORDER TENSOR DEFINED AS CVI.C.CVI 
C                   AND STORED AS NTENS VECTOR
C       CVIDOT    - 1ST ORDER TIME DERIVATIVE OF CVI
C       CTENS     - DEVIATORIC RIGHT CAUCHY-GREEN TENSOR STORED AS 
C                   NTENS VECTOR     
C
C COMPUTE CVICCVI AND INVARIANTS I1E, I2E, I1V
C     
      CALL KCVICCVI(CVICCVI,CTENS,CVI,NSHR,NTENS)
      CALL KAI1E(AI1E,CTENS,CVI,NSHR,NTENS)
      CALL KAI2E(AI2E,CVICCVI,CTENS,AI1E,NSHR,NTENS)
      CALL KAI1V(AI1V,CVI,NSHR,NTENS)     
C
C MATERIAL PARAMETERS
C         
      AM1  = PROPS(6)  
      AA1  = PROPS(7)  
      AM2  = PROPS(8)  
      AA2  = PROPS(9)       
      AETA0  = PROPS(10)  
      AETAI  = PROPS(11)  
      ABETA1 = PROPS(12)  
      ABETA2 = PROPS(13)  
      AK1    = PROPS(14)  
      AK2    = PROPS(15)    
C
C 1ST DERIVATIVE OF NON-EQUILIBRIUM BRANCH AND CAUCHY STRESS INVARIANT
C 
      DUNEDI1E = (3**(1 - AA1)*AI1E**(-1 + AA1)*AM1)/2. + 
     1           (3**(1 - AA2)*AI1E**(-1 + AA2)*AM2)/2.        
      AJ2NE=(AI1E**2./3.-AI2E)*(2.0*DUNEDI1E)**2.
C
C COMPUTE AETAK
C 
      AETAK=AETAI+(AETA0-AETAI+AK1*(AI1V**ABETA1-3.**ABETA1))/
     1      (1.+(AK2*AJ2NE)**ABETA2) 
C
C COMPUTE CVIDOT
C
      TAUI = 2.0*DUNEDI1E/AETAK           
      DO I= 1,NTENS
         CVIDOT(I)=TAUI*(AI1E/3.*CVI(I)-CVICCVI(I))
      END DO
C                  
      RETURN
      END      
C
C***********************************************************************
C
      SUBROUTINE KRKCORRECTOR(YKPONE,YK,YKPONEHAT,NSHR,NTENS)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION YKPONE(NTENS),YK(NTENS),YKPONEHAT(NTENS)
C
      CALL KDETVECTOR(DETYKPONEHAT,YKPONEHAT,NSHR,NTENS)
      AFAC=(1.0/DETYKPONEHAT)**(1.0/3.0)
C      
      DO I=1,NTENS
        YKPONE(I)=AFAC*YKPONEHAT(I)
      ENDDO
C      
      RETURN
      END      
C
C***********************************************************************
C
      SUBROUTINE KBBAR(BBAR, F, DETF, NSHR, NTENS)
C  
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION F(3,3), FBAR(3,3), BBAR(NTENS)
C
C     RETURNS THE DEVIATORIC LEFT CAUCHY-GREEN TENSOR AS A NTENS 
C     VECTOR
C
C    JACOBIAN AND DISTORTION TENSOR
C      
      AA=DETF**(-1.0/3.0)
      DO I=1, 3
        DO J=1, 3
          FBAR(I, J)=AA*F(I, J)
        END DO
      END DO
C
C    COMPUTE LEFT CAUCHY-GREEN TENSOR
C
      CALL KAATTENS(BBAR, FBAR, NSHR, NTENS)  
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE KCBAR(CBAR, F, DETF, NSHR, NTENS)
C  
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION F(3,3), FBAR(3,3), CBAR(NTENS)
C
C     RETURNS THE DEVIATORIC RIGHT CAUCHY-GREEN TENSOR AS A NTENS 
C     VECTOR      
C
C    JACOBIAN AND DISTORTION TENSOR
C      
      AA=DETF**(-1.0/3.0)
      DO I=1, 3
        DO J=1, 3
          FBAR(I, J)=AA*F(I, J)
        END DO
      END DO
C
C    COMPUTE RIGHT CAUCHY-GREEN TENSOR
C
      CALL KATATENS(CBAR, FBAR, NSHR, NTENS)  

      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE KBEBAR(BEBARTENS, F, CVI, DETF, NSHR, NTENS)
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION BEBARTENS(NTENS),CVI(NTENS),F(3,3)
C     
C    RETURNS BEBARTENS, SYMMETRIC SECOND ORDER TENSOR DEFINED
C    AS FBAR.CVI.FBAR^T AND STORED AS NTENS VECTOR
C      
      AFAC=(1.0/DETF)**(2.0/3.0)
      IF (NSHR.EQ.3) THEN
       BEBARTENS(1)=(F(1,3)*(CVI(5)*F(1,1)+CVI(6)*F(1,2)+CVI(3)*F(1,3))+ 
     1         F(1,1)*(CVI(1)*F(1,1) + CVI(4)*F(1,2) + CVI(5)*F(1,3)) + 
     2     F(1,2)*(CVI(4)*F(1,1) + CVI(2)*F(1,2) + CVI(6)*F(1,3)))*AFAC
       BEBARTENS(2)=(F(2,3)*(CVI(5)*F(2,1)+CVI(6)*F(2,2)+CVI(3)*F(2,3))+ 
     1         F(2,1)*(CVI(1)*F(2,1) + CVI(4)*F(2,2) + CVI(5)*F(2,3)) + 
     2     F(2,2)*(CVI(4)*F(2,1) + CVI(2)*F(2,2) + CVI(6)*F(2,3)))*AFAC
       BEBARTENS(3)=(F(3,3)*(CVI(5)*F(3,1)+CVI(6)*F(3,2)+CVI(3)*F(3,3))+ 
     1         F(3,1)*(CVI(1)*F(3,1) + CVI(4)*F(3,2) + CVI(5)*F(3,3)) + 
     2     F(3,2)*(CVI(4)*F(3,1) + CVI(2)*F(3,2) + CVI(6)*F(3,3)))*AFAC
       BEBARTENS(4)=((CVI(1)*F(1,1)+CVI(4)*F(1,2)+CVI(5)*F(1,3))*F(2,1)+ 
     1         (CVI(4)*F(1,1) + CVI(2)*F(1,2) + CVI(6)*F(1,3))*F(2,2)+ 
     2     (CVI(5)*F(1,1) + CVI(6)*F(1,2) + CVI(3)*F(1,3))*F(2,3))*AFAC
       BEBARTENS(5)=((CVI(1)*F(1,1)+CVI(4)*F(1,2)+CVI(5)*F(1,3))*F(3,1)+ 
     1           (CVI(4)*F(1,1) + CVI(2)*F(1,2) + CVI(6)*F(1,3))*F(3,2)+ 
     2     (CVI(5)*F(1,1) + CVI(6)*F(1,2) + CVI(3)*F(1,3))*F(3,3))*AFAC
       BEBARTENS(6)=((CVI(1)*F(2,1)+CVI(4)*F(2,2)+CVI(5)*F(2,3))*F(3,1)+ 
     1           (CVI(4)*F(2,1) + CVI(2)*F(2,2) + CVI(6)*F(2,3))*F(3,2)+ 
     2      (CVI(5)*F(2,1) + CVI(6)*F(2,2) + CVI(3)*F(2,3))*F(3,3))*AFAC
      ELSE
        BEBARTENS(1)=(F(1,2)*(CVI(4)*F(1,1) + CVI(2)*F(1,2)) + 
     1                F(1,1)*(CVI(1)*F(1,1) + CVI(4)*F(1,2)))*AFAC
        BEBARTENS(2)=(F(2,2)*(CVI(4)*F(2,1) + CVI(2)*F(2,2)) + 
     1                F(2,1)*(CVI(1)*F(2,1) + CVI(4)*F(2,2)))*AFAC
        BEBARTENS(3)=(CVI(3)*F(3,3)**2.0)*AFAC
        BEBARTENS(4)=((CVI(1)*F(1,1) + CVI(4)*F(1,2))*F(2,1) +
     1                (CVI(4)*F(1,1) + CVI(2)*F(1,2))*F(2,2))*AFAC       
      END IF
C      
      RETURN
      END 
C
C***********************************************************************
C
      SUBROUTINE KRKCBARTENS(RKCBARTENS,DFGRD0,DFGRD1,RKC,NSHR,NTENS)
C  
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION DFGRD0(3,3), DFGRD1(3,3), RKCBARTENS(NTENS), DFGRINT(3,3)
C
C     RETURNS THE DEVIATORIC RIGHT CAUCHY-GREEN TENSOR USED IN THE R-K
C     INTEGRATION SCHEME AS A NTENS VECTOR
C
      DO I= 1,3
        DO J= 1,3
          DFGRINT(I,J)=DFGRD0(I,J)+RKC*(DFGRD1(I,J)-DFGRD0(I,J))
        END DO
      END DO   
C            
      CALL KDET(DETDFGRINT, DFGRINT, NSHR)      
C
      CALL KCBAR(RKCBARTENS, DFGRINT, DETDFGRINT, NSHR, NTENS)
      RETURN
      END  
C
C***********************************************************************
C
      SUBROUTINE KCVICCVI(CVICCVI,CTENS,CVI,NSHR,NTENS)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION CVICCVI(NTENS),CTENS(NTENS),CVI(NTENS)      
C
C     RETURNS THE SYMMETRIC SECOND ORDER TENSOR CVI.C.CVI STORED AS  
C     NTENS VECTOR
C
      IF(NSHR.EQ.3) THEN
        CVICCVI(1)=0.0+
     1    CVI(1)*(CTENS(1)*CVI(1) + CTENS(4)*CVI(4) + CTENS(5)*CVI(5)) + 
     2    CVI(4)*(CVI(1)*CTENS(4) + CVI(4)*CTENS(2) + CVI(5)*CTENS(6)) + 
     3    CVI(5)*(CVI(1)*CTENS(5) + CVI(4)*CTENS(6) + CVI(5)*CTENS(3))
        CVICCVI(2)=0.0+
     1    CVI(4)*(CTENS(1)*CVI(4) + CTENS(4)*CVI(2) + CTENS(5)*CVI(6)) + 
     2    CVI(2)*(CTENS(4)*CVI(4) + CTENS(2)*CVI(2) + CTENS(6)*CVI(6)) + 
     3    CVI(6)*(CVI(4)*CTENS(5) + CVI(2)*CTENS(6) + CVI(6)*CTENS(3))
        CVICCVI(3)=0.0+
     1    CVI(5)*(CTENS(1)*CVI(5) + CTENS(4)*CVI(6) + CTENS(5)*CVI(3)) + 
     2    CVI(6)*(CTENS(4)*CVI(5) + CTENS(2)*CVI(6) + CTENS(6)*CVI(3)) + 
     3    CVI(3)*(CTENS(5)*CVI(5) + CTENS(6)*CVI(6) + CTENS(3)*CVI(3))
        CVICCVI(4)=0.0+
     1    CVI(4)*(CTENS(1)*CVI(1) + CTENS(4)*CVI(4) + CTENS(5)*CVI(5)) + 
     2    CVI(2)*(CVI(1)*CTENS(4) + CVI(4)*CTENS(2) + CVI(5)*CTENS(6)) + 
     3    CVI(6)*(CVI(1)*CTENS(5) + CVI(4)*CTENS(6) + CVI(5)*CTENS(3))
        CVICCVI(5)=0.0+
     1    CVI(5)*(CTENS(1)*CVI(1) + CTENS(4)*CVI(4) + CTENS(5)*CVI(5)) + 
     2    (CVI(1)*CTENS(4) + CVI(4)*CTENS(2) + CVI(5)*CTENS(6))*CVI(6) + 
     3    (CVI(1)*CTENS(5) + CVI(4)*CTENS(6) + CVI(5)*CTENS(3))*CVI(3)
        CVICCVI(6)=0.0+
     1    CVI(5)*(CTENS(1)*CVI(4) + CTENS(4)*CVI(2) + CTENS(5)*CVI(6)) + 
     2    CVI(6)*(CTENS(4)*CVI(4) + CTENS(2)*CVI(2) + CTENS(6)*CVI(6)) + 
     3    (CVI(4)*CTENS(5) + CVI(2)*CTENS(6) + CVI(6)*CTENS(3))*CVI(3)        
      ELSE
        CVICCVI(1)=0.0+
     1    CVI(1)*(CTENS(1)*CVI(1) + CTENS(4)*CVI(4)) + 
     2    CVI(4)*(CVI(1)*CTENS(4) + CVI(4)*CTENS(2))
        CVICCVI(2)=0.0+
     1    CVI(4)*(CTENS(1)*CVI(4) + CTENS(4)*CVI(2)) + 
     2    CVI(2)*(CTENS(4)*CVI(4) + CTENS(2)*CVI(2))
        CVICCVI(3)=CTENS(3)*CVI(3)**2.0
        CVICCVI(4)=0.0+
     1    CVI(4)*(CTENS(1)*CVI(1) + CTENS(4)*CVI(4)) + 
     2    CVI(2)*(CVI(1)*CTENS(4) + CVI(4)*CTENS(2))
      END IF
C
      RETURN
      END
C     
C***********************************************************************
C
      SUBROUTINE KAI1E(AI1E,CTENS,CVI,NSHR,NTENS)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION CTENS(NTENS),CVI(NTENS)
C
C     RETURNS THE FIRST INVARIANT OF C.CVI
C
      IF(NSHR.EQ.3) THEN
        AI1E=CTENS(1)*CVI(1) + 2.0*CTENS(4)* CVI(4) + 
     1   CTENS(2)*CVI(2) + 2.0*CTENS(5)*CVI(5) + 
     2   CTENS(3)*CVI(3) + 2.0*CTENS(6)*CVI(6)     
      ELSE
        AI1E=CTENS(1)*CVI(1) + 2.0*CTENS(4)* CVI(4) + 
     1   CTENS(2)*CVI(2) + CTENS(3)*CVI(3)   
      END IF
C
      RETURN
      END
C     
C***********************************************************************
C
      SUBROUTINE KAI2E(AI2E,CVICCVI,CTENS,AI1E,NSHR,NTENS)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION CVICCVI(NTENS),CTENS(NTENS)
C
C     RETURNS THE SECOND INVARIANT OF C.CVI
C
      CALL KAI1E(ATRCVICCVIC,CTENS,CVICCVI,NSHR,NTENS)
      AI2E = (AI1E**2.0-ATRCVICCVIC)/2.0
C
      RETURN
      END
C     
C***********************************************************************
C
      SUBROUTINE KAI1V(AI1V,CVI,NSHR,NTENS)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION CVI(NTENS)
C
C     RETURNS THE FIRST INVARIANT OF CV IN TERMS OF CVI
C
      IF(NSHR.EQ.3) THEN      
        CALL KDETVECTOR(DETCVI,CVI,NSHR,NTENS)
        AI1V = (CVI(1)*CVI(2)+CVI(2)*CVI(3)+CVI(3)*CVI(1)-
     1          CVI(4)**2.-CVI(5)**2.-CVI(6)**2.)/DETCVI
      ELSE
        CALL KDETVECTOR(DETCVI,CVI,NSHR,NTENS)
        AI1V = (CVI(1)*CVI(2)+CVI(2)*CVI(3)+CVI(3)*CVI(1)-
     1          CVI(4)**2.)/DETCVI
      END IF
C
      RETURN
      END   
C
C***********************************************************************
C
      SUBROUTINE KCHECKS(STATEV, NSTATV, NTENS, PROPS, NPROPS)
C  
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#INCLUDE <SMAASPUSERSUBROUTINES.HDR>  
      DIMENSION STATEV(NSTATV), PROPS(NPROPS)
C
C     PERFORMS CHECKS ON THE INPUT TO THE SUBROUTINE
C
C
C     STDB_ABQERR AND GET_THREAD_ID INITIALIZATION
C
      DIMENSION INTV(2),REALV(6)
      CHARACTER*8 CHARV(1)
      CHARACTER*100 STRING1, STRING2, STRING3, STRING4, STRING5
      CHARACTER*500 STRING
C      
      INTEGER MYTHREADID  
C      
      INTV(1)=0
      INTV(2)=0
      REALV(1)=0.
      REALV(2)=0.
      REALV(3)=0.
      REALV(4)=0.
      REALV(5)=0.
      REALV(6)=0.
      CHARV(1)=''
C
      MYTHREADID = GET_THREAD_ID()
C
C     INPUT CHECKS
C   
      IF (MYTHREADID.EQ.0) THEN
        IF (NSTATV.NE.NTENS) THEN  
          INTV(1)=NSTATV
          INTV(2)=NTENS
          STRING1='RECEIVED NSTATV = %I BUT NSTATV SHOULD BE EQUAL TO'
          STRING2='%I. UPDATE *DEPVAR KEYWORD IN THE INPUT FILE.'
          STRING = TRIM(STRING1) // ' ' // TRIM(STRING2)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        ELSE IF (NSTATV.EQ.4) THEN
          IF ((STATEV(1).NE.1.0).OR.(STATEV(2).NE.1.0).OR.
     1       (STATEV(3).NE.1.0).OR.(STATEV(4).NE.0.0)) THEN    
            REALV(1)=STATEV(1)  
            REALV(2)=STATEV(2)
            REALV(3)=STATEV(3)  
            REALV(4)=STATEV(4)  
            STRING1='RECEIVED (%R, %R, %R, %R) AS INITIAL CONDITIONS'
            STRING2='FOR THE INTERNAL VARIABLE STATEV. THE INITIAL'
            STRING3='CONDITIONS MUST BE (1.0,1.0,1.0,0.0). UPDATE THE'
            STRING4='*INITIAL CONDITIONS KEYWORD IN THE INPUT FILE.'
            STRING = TRIM(STRING1) // ' ' // TRIM(STRING2) // ' ' // 
     1               TRIM(STRING3) // ' ' // TRIM(STRING4)
            CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
          END IF
        ELSE IF (NSTATV.EQ.6) THEN
          IF ((STATEV(1).NE.1.0).OR.(STATEV(2).NE.1.0).OR.
     1       (STATEV(3).NE.1.0).OR.(STATEV(4).NE.0.0).OR.
     2       (STATEV(5).NE.0.0).OR.(STATEV(6).NE.0.0)) THEN    
            REALV(1)=STATEV(1)  
            REALV(2)=STATEV(2)
            REALV(3)=STATEV(3)  
            REALV(4)=STATEV(4)  
            REALV(5)=STATEV(5)  
            REALV(6)=STATEV(6)  
            STRING1='RECEIVED (%R, %R, %R, %R, %R, %R) AS INITIAL'
            STRING2='CONDITIONS FOR THE INTERNAL VARIABLE STATEV.'
            STRING3='THE INITIAL CONDITIONS MUST BE (1.0,1.0,1.0,'
            STRING4='0.0,0.0,0.0). UPDATE THE *INITIAL CONDITIONS'
            STRING5='KEYWORD IN THE INPUT FILE.'
            STRING = TRIM(STRING1) // ' ' // TRIM(STRING2) // ' ' // 
     1               TRIM(STRING3) // ' ' // TRIM(STRING4) // ' ' //
     2               TRIM(STRING5)
            CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
          END IF
        END IF
C
        IF (NPROPS.NE.15) THEN    
          INTV(1)=NPROPS  
          STRING1='RECEIVED %I MATERIAL PROPERTIES. THE SUBROUTINE'
          STRING2='REQUIRES 15 MATERIAL PROPERTIES.'
          STRING = TRIM(STRING1) // ' ' //  TRIM(STRING2)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        END IF
      END IF     
C
C     MATERIAL PARAMETERS
C  
      AMU1    = PROPS(1) 
      AALPHA1 = PROPS(2) 
      AMU2    = PROPS(3) 
      AALPHA2 = PROPS(4) 
      AKAPPAE = PROPS(5)  
      AM1     = PROPS(6) 
      AA1     = PROPS(7) 
      AM2     = PROPS(8) 
      AA2     = PROPS(9) 
      AETA0   = PROPS(10)
      AETAI   = PROPS(11)
      ABETA1  = PROPS(12)
      ABETA2  = PROPS(13)
      AK1     = PROPS(14)
      AK2     = PROPS(15)  
C
C     PARTIAL CHECKS OF MATERIAL PARAMETERS
C  
      IF (MYTHREADID.EQ.0) THEN
        IF ((AMU1.LT.0.).OR.(AMU2.LT.0.).OR.(AMU1+AMU2.LE.0.)) THEN
          REALV(1)=AMU1      
          REALV(2)=AMU2      
          REALV(3)=AMU1+AMU2      
          STRING1='RECEIVED AMU1 = %R AND AMU2 = %R, AMU1 + AMU2 = %R.'
          STRING2='THE PARAMETERS AMU1 AND AMU2 MUST BE NON-NEGATIVE'
          STRING3='AND AMU1 + AMU2 MUST BE GREATER THAN ZERO.'
          STRING = TRIM(STRING1) // ' ' // TRIM(STRING2) // ' ' // 
     1             TRIM(STRING3)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        END IF           
C        
        IF ((AALPHA1.EQ.0.).OR.(AALPHA2.EQ.0.).OR.
     1   (MAX(AALPHA1,AALPHA2).LE.0.)) THEN
          REALV(1)=AALPHA1      
          REALV(2)=AALPHA2      
          STRING1='RECEIVED AALPHA1 = %R AND AALPHA2 = %R.'
          STRING2='THE EXPONENTS AALPHA1 AND AALPHA2 MUST BE NON-ZERO'
          STRING3='AND MAX(AALPHA1,AALPHA2) MUST BE GREATER THAN ZERO.'
          STRING = TRIM(STRING1) // ' ' // TRIM(STRING2) // ' ' // 
     1             TRIM(STRING3)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        END IF 
C         
        IF (AKAPPAE.LT.0.) THEN
          REALV(1)=AKAPPAE     
          STRING1='RECEIVED AKAPPAE = %R.  THE PARAMETER AKAPPAE'
          STRING2='MUST BE GREATER THAN ZERO.'
          STRING = TRIM(STRING1) // ' ' // TRIM(STRING2)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        END IF           
C          
        IF ((AM1.LT.0.).OR.(AM2.LT.0.).OR.(AM1+AM2.LE.0.)) THEN
          REALV(1)=AM1      
          REALV(2)=AM2      
          REALV(3)=AM1+AM2      
          STRING1='RECEIVED AM1 = %R AND AM2 = %R, AM1 + AM2 = %R.'
          STRING2='THE PARAMETERS AM1 AND AM2 MUST BE NON-NEGATIVE'
          STRING3='AND AM1 + AM2 MUST BE GREATER THAN ZERO.'
          STRING = TRIM(STRING1) // ' ' // TRIM(STRING2) // ' ' // 
     1             TRIM(STRING3)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        END IF           
C        
        IF ((AA1.EQ.0.).OR.(AA2.EQ.0.).OR.(MAX(AA1,AA2).LE.0.)) THEN
          REALV(1)=AA1      
          REALV(2)=AA2      
          STRING1='RECEIVED AA1 = %R AND AA2 = %R.'
          STRING2='THE EXPONENTS AA1 AND AA2 MUST BE NON-ZERO'
          STRING3='AND MAX(AA1,AA2) MUST BE GREATER THAN ZERO.'
          STRING = TRIM(STRING1) // ' ' // TRIM(STRING2) // ' ' // 
     1             TRIM(STRING3)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        END IF  
C
        IF ((AETA0.LT.0.).OR.(AETAI.LT.0.)) THEN
          REALV(1)=AETA0      
          REALV(2)=AETAI        
          STRING1='RECEIVED AETA0 = %R AND AETAI = %R. THE '
          STRING2='PARAMETERS AETA0 AND AETAI MUST BE NON-NEGATIVE.'
          STRING = TRIM(STRING1) // ' ' // TRIM(STRING2)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        END IF  
C
        IF ((ABETA1.LT.0.).OR.(ABETA2.LT.0.)) THEN
          REALV(1)=ABETA1      
          REALV(2)=ABETA2        
          STRING1='RECEIVED ABETA1 = %R AND ABETA2 = %R. THE'
          STRING2='PARAMETERS ABETA1 AND ABETA2 MUST BE NON-NEGATIVE.'
          STRING = TRIM(STRING1) // ' ' // TRIM(STRING2)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        END IF  
C
        IF ((AK1.LT.0.).OR.(AK2.LT.0.)) THEN
          REALV(1)=AK1      
          REALV(2)=AK2        
          STRING1='RECEIVED AK1 = %R AND AK2 = %R. THE'
          STRING2='PARAMETERS AK1 AND AK2 MUST BE NON-NEGATIVE.'
          STRING = TRIM(STRING1) // ' ' // TRIM(STRING2)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        END IF    
      END IF           
      RETURN
      END                
C  
C
C***********************************************************************
C
      SUBROUTINE KBDEVA(DEVA, A, TRA, NSHR, NTENS)
C  
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION DEVA(NTENS), A(NTENS)
C
C     RETURNS THE TRACELESS SYMMETRIC TENSOR A-TR(A)/3.ID AS A NTENS 
C     VECTOR
C
      DEVA(1)=A(1)-TRA/3.0
      DEVA(2)=A(2)-TRA/3.0
      DEVA(3)=A(3)-TRA/3.0
      DEVA(4)=A(4)
      IF(NSHR.EQ.3) THEN
        DEVA(5)=A(5)
        DEVA(6)=A(6)
      END IF    
      
      RETURN
      END  
C
C***********************************************************************
C
      SUBROUTINE KAATTENS(AAT, A, NSHR, NTENS)
C  
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION A(3,3), AAT(NTENS)
C
C     RETURNS THE SYMMETRIC TENSOR A.AT (SIMILAR TO B, THE LEFT CAUCHY-GREEN TENSOR)
c     AS A NTENS VECTOR
C
      AAT(1)=A(1, 1)**2.0+A(1, 2)**2.0+A(1, 3)**2.0
      AAT(2)=A(2, 1)**2.0+A(2, 2)**2.0+A(2, 3)**2.0
      AAT(3)=A(3, 3)**2.0+A(3, 1)**2.0+A(3, 2)**2.0
      AAT(4)=A(1, 1)*A(2, 1)+A(1, 2)*A(2, 2)
     1       +A(1, 3)*A(2, 3)
      IF(NSHR.EQ.3) THEN
        AAT(5)=A(1, 1)*A(3, 1)+A(1, 2)*A(3, 2)
     1         +A(1, 3)*A(3, 3)
        AAT(6)=A(2, 1)*A(3, 1)+A(2, 2)*A(3, 2)
     1         +A(2, 3)*A(3, 3)
      END IF         
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE KATATENS(ATA, A, NSHR, NTENS)
C  
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION A(3,3), ATA(NTENS)
C
C     RETURNS THE SYMMETRIC TENSOR AT.A (SIMILAR TO C, THE RIGHT CAUCHY-GREEN TENSOR)
C     AS A NTENS VECTOR
C
      ATA(1)=A(1, 1)**2.0+A(2, 1)**2.0+A(3, 1)**2.0
      ATA(2)=A(1, 2)**2.0+A(2, 2)**2.0+A(3, 2)**2.0
      ATA(3)=A(1, 3)**2.0+A(2, 3)**2.0+A(3, 3)**2.0
      ATA(4)=A(1, 1)*A(1, 2)+A(2, 1)*A(2, 2)
     1       +A(3, 1)*A(3, 2)
      IF(NSHR.EQ.3) THEN
        ATA(5)=A(1, 1)*A(1, 3)+A(2, 1)*A(2, 3)
     1         +A(3, 1)*A(3, 3)
        ATA(6)=A(1, 2)*A(1, 3)+A(2, 2)*A(2, 3)
     1         +A(3, 2)*A(3, 3)
      END IF         
      RETURN
      END  
C
C***********************************************************************
C
      SUBROUTINE KINSERTCOLUMN(ACOL, A, NCOL, NASZ1, NASZ2)
C            
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION ACOL(NASZ1), A(NASZ1,NASZ2)            
C
C     COPIES COMPONENTS OF VECTOR ACOL TO COLUMN #NCOL OF ARRAY A
C 
      DO I= 1,NASZ1
         A(I,NCOL)=ACOL(I)
      END DO
C      
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE KDET(DETA,A,NSHR)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(3,3)
C
C     RETURNS THE DETERMINANT OF 3X3 MATRIX
C
      IF (NSHR.EQ.3) THEN
         DETA = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2))-
     1          A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1))+
     2          A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
      ELSE
         DETA = A(3,3)*(A(1,1)*A(2,2) - A(1,2)*A(2,1))      
      END IF

      RETURN
      END      
C
C***********************************************************************
C
      SUBROUTINE KTRACEVECTOR(TRA, A, NDI, NTENS)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION A(NTENS)           
C
C     RETURNS THE TRACE TRA OF SYMMETRIC TENSOR A STORED AS NTENS VECTOR
C 
      TRA=0.0
      DO I=1,NDI
        TRA=TRA+A(I)
      ENDDO
C
      RETURN
      END   
C
C***********************************************************************
C
      SUBROUTINE KDETVECTOR(DETA,A,NSHR,NTENS)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION A(NTENS)         
C
C     RETURNS THE DETERMINANT DETA OF SYMMETRIC TENSOR A STORED AS NTENS 
C     VECTOR
C 
      IF(NSHR.EQ.3) THEN
        DETA=A(1)*A(2)*A(3)+2.0*A(4)*A(5)*A(6)
     1      -A(3)*A(4)**2.0-A(2)*A(5)**2.0-A(1)*A(6)**2.0
      ELSE
        DETA=A(1)*A(2)*A(3)-A(3)*A(4)**2.0
      END IF
      IF (DABS(DETA).LT.1.D-30) THEN
         WRITE(*,*) 'SINGULAR DET FOUND - ** PROGRAM STOP **'
         STOP
      ENDIF
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE KINITIA(A,N)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(N)
C
C     RETURNS ZERO VECTOR OF DIMENSION N
C
      DO I=1,N
        A(I)=0.0
      END DO
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
C
C     USER OUTPUT VARIABLE FOR THE DETERMINANT OF THE DEFORMATION 
C     GRADIENT
C      
      CALL GETVRM('DGP',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     1 MATLAYO,LACCFLA)
      DETF=ARRAY(1)*ARRAY(2)*ARRAY(3)
      UVAR(1) = DETF
C
      RETURN
      END