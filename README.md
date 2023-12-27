# UMAT_Lefevre_Sozio_Lopez-Pamies

## Description
This Abaqus UMAT [subroutine](UMAT_KLP_RK5_hybrid.f) implements a family of internal-variable-based constitutive models for the viscoelastic response of elastomers of any compressibility — including fully incompressible elastomers — undergoing finite deformations, that introduced by Kumar and  Lopez-Pamies [1]. The models can be chosen to account for a wide range of non-Gaussian elasticities, as well as for a wide range of nonlinear viscosities. 

The present [subroutine](UMAT_KLP_RK5_hybrid.f) is implemented for the choice of two-term Lopez-Pamies non-Gaussian elasticities [2] and deformation-enhanced shear-thinning nonlinear viscosity proposed in [1]. It relies on a hybrid formulation [3] to deal with nearly or fully incompressible elastomers.

## Usage
The [subroutine](UMAT_KLP_RK5_hybrid.f) is to be used as a `USER MATERIAL` with 15 material `CONSTANTS` and a `TOTAL HYBRID FORMULATION`, e.g.,
```
*USER MATERIAL,CONSTANTS=15,HYBRID FORMULATION=TOTAL
```
in the input (.inp) file. Due to the hybrid nature of the underlying formulation, hybrid finite elements must be used in the simulation.

### Material properties
The 15 material properties of the model are the primary inputs of the [subroutine](UMAT_KLP_RK5_hybrid.f) and are listed in order below:

- `AMU1`: parameter #1 of the equilibrium elasticity
- `AALPHA1`: exponent #1 of the equilibrium elasticity
- `AMU2`: parameter #2 of the equilibrium elasticity
- `AALPHA2`: exponent #2 of the equilibrium elasticity
- `AKAPPAE`: bulk modulus of the equilibrium elasticity 
- `AA1`: exponent #1 of the non-equilibrium elasticity
- `AM2`: parameter #2 of the non-equilibrium elasticity
- `AA2`: exponent #2 of the non-equilibrium elasticity
- `AM1`: parameter #1 of the non-equilibrium elasticity
- `AETA0`: parameter #1 of the deviatoric viscosity
- `AETAI`: parameter #2 of the deviatoric viscosity
- `ABETA1`: exponent #1 of the deviatoric viscosity
- `ABETA2`: exponent #2 of the deviatoric viscosity
- `AK1`: factor #1 of the deviatoric viscosity
- `AK2`: factor #2 of the deviatoric viscosity

They are entered in the input (.inp) file as follows:
```
*USER MATERIAL,CONSTANTS=15,HYBRID FORMULATION=TOTAL
AMU1, AALPHA1, AMU2, AALPHA2, AKAPPAE, AM1, AA1, AM2,
AA2, AETA0, AETAI, ABETA1, ABETA2, AK1, AK2
```
  
The material parameters `AMU1`, `AMU2` characterizing the equilibrium elasticity are non-negative real numbers (`AMU1` ≥ 0, `AMU2` ≥ 0) with strictly positive sum (`AMU1` + `AMU2` > 0). 

The two exponents `AALPHA1`, `AALPHA2` are non-zero real numbers (`AALPHA1` ≠ 0, `AALPHA2` ≠ 0) with MAX{`AALPHA1`, `AALPHA2`} > 0 leading to  a strongly elliptic strain energy (see eq. (22) in [2]). This is left to the user to check.

The material parameter `AKAPPAE` characterizing the compressibility of the equilibrium elasticity is a strictly positive real number (`AKAPPAE` > 0).
The material parameters `AM1`, `AM2` characterizing the non-equilibrium elasticity are non-negative real numbers (`AM1` ≥ 0, `AM2` ≥ 0) with strictly positive sum (`AM1` + `AM2` > 0).

The two exponents `AA1`, `AA2` are non-zero real numbers (`AA1` ≠ 0, `AA2` ≠ 0) with MAX{`AA1`, `AA2`} > 0 leading to a strongly elliptic strain energy (see eq. (22) in [2]). This is left to the user to check.

The material parameters `AETA0`, `AETAI`, `ABETA1`, `ABETA2`, `AK1`, `AK2` characterizing the viscous dissipation are all non-negative real numbers (`AETA0` ≥ 0, `AETAI` ≥ 0, `ABETA1` ≥ 0, `ABETA2` ≥ 0, `AK1` ≥ 0, `AK2` ≥ 0).

### Solution-dependent state variables
- The [subroutine](UMAT_KLP_RK5_hybrid.f) uses solution-dependent state variables for the independent components of the symmetric tensorial internal variable ${\bold C}^{v-1}$. The space for the state variables must be allocated in the input file using the `*DEPVAR` keyword:
    ```
    *DEPVAR
    NSTATV
    ```
    where `NSTATV` is the number of state variables required. `NSTATV` depends on the type of finite elements used:
    - `NSTATV` = 4 for elements with only 1 non-zero shear component
    (CPE3H, CPE4H, CPE6H, CPE8H, CAX3H, CAX4H, CAX6H, CAX8H, CPEG3H, 
    CPEG4H, CPEG6H, CPEG8H, ...)
    - `NSTATV` = 6 for elements with all 3 non-zero shear components 
    (C3D4H, C3D8H, C3D10H, C3D20H, CGAX3H, CGAX4H, CGAX6H, CGAX8H, ...)

    Note that an incorrect value of `NSTATV` will trigger an error in the subroutine and the error message in the .msg file will indicate the expected value for `NSTATV`.

- The values of state variables must be initialized according to the initial condition ${\bold C}^{v-1}={\bold I}$ using the `*INITIAL CONDITIONS` keyword as follows:
    - `NSTATV` = 4: 
        ```
        *INITIAL CONDITIONS, TYPE=SOLUTION
        ALLELEMENTSET,1.0,1.0,1.0,0.0    
        ```
    - `NSTATV` = 6:
        ```
        *INITIAL CONDITIONS, TYPE=SOLUTION
        ALLELEMENTSET,1.0,1.0,1.0,0.0,0.0,0.0
        ```
    where `ALLELEMENTSET` is an element set of all finite elements     associated with the UMAT material.

### Energy outputs

- The [subroutine](UMAT_KLP_RK5_hybrid.f) uses the variable `SSE` to output $\psi^E = \Psi^E / J$, the strain energy density per deformed unit volume in the equilibrium branch. 

    The corresponding contour plot can be obtained by including the output variable identifier `SENER` in the element output request keyword `*ELEMENT OUTPUT`. 

    The corresponding total energy (i.e., the volume internal of `SSE`) can be obtained by including the output variable identifier `ALLSE` in the history output request keyword `*ENERGY OUTPUT`.

- The [subroutine](UMAT_KLP_RK5_hybrid.f) uses the variable `SPD` to output $\psi^{NE} = \Psi^{NE} / J$, the strain energy density per deformed unit volume in the non-equilibrium branch. Note that this use of the `SPD` variable is different from the typical meaning of the `SPD` variable in UMAT subroutines, which is the plastic dissipation density.

    The corresponding contour plot can be obtained by including the output variable identifier `PENER` in the element output request keyword `*ELEMENT OUTPUT`. 

    The corresponding total energy (i.e., the volume internal of `SPD`) can be obtained by including the output variable identifier `ALLPD` in the history output request keyword `*ENERGY OUTPUT`. 

### Contour plots
- Contour plots of the components of the internal variable ${\bold C}^{v-1}$ can be obtained by including the output variable identifier `SDV` in the element output request keyword `*ELEMENT OUTPUT`. 

- A UVARM subroutine is also included to create two user output variable, `UVARM1` = $I_1$ = $\text{tr}({\bold F}^T{\bold F})$ and `UVARM2` = $J$ = $\text{det(F)}$. The space for these user output variables must be allocated in the input file using the `*USER OUTPUT VARIABLES` keyword:
    ```
    *USER OUTPUT VARIABLES
    2
    ```

    The corresponding contour plots can be obtained by including the output variable identifier `UVARM` in the element output request keyword `*ELEMENT OUTPUT`. 

### Additional information
- This [subroutine](UMAT_KLP_RK5_hybrid.f) doesn't create predefined field variables. 
- Please consult the Abaqus Documentation for additional references regarding the use of the UMAT subroutine for (near-)incompressible material models.

## Examples
The [Examples](Examples) folder contains examples of how to use the present [subroutine](UMAT_KLP_RK5_hybrid.f).

## References
The references below can be found in the [References](References) folder.

1. Kumar, A., Lopez-Pamies, O., 2016. [On the two-potential constitutive modeling of rubber viscoelastic materials](http://pamies.cee.illinois.edu/Publications_files/CRM_2016.pdf). *Comptes Rendus Mécanique* 344, 102–112. 
2. Lopez-Pamies, O., 2010. [A new I1-based hyperelastic model for rubber elastic materials](http://pamies.cee.illinois.edu/Publications_files/CMAME2019.pdf). *Comptes Rendus Mécanique* 338, 3–11.
3. Lefèvre, V., Sozio, F., Lopez-Pamies, O., 2023. [Abaqus implementation of a large family of finite viscoelasticity models](https://arxiv.org/pdf/2311.13751.pdf). Submitted.