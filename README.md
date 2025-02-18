# User manual
---
## 1. Introduction
### 1.1 Overview
This program implements high-fidelity steady-state and transient Neutron/Thermal-hydraulics coupling between RMC and ANSYS Fluent. 
- RMC is a 3D neutron transport code for reactors based on the Monte Carlo method.
- ANSYS Fluent is a CFD-based thermal-hydraulics code.

Theoretically, this program can be quickly applied to any type of reactors with various geometries and materials.
### 1.2 Prerequisites
* Operating System Requirements: CentOS 7, Ubuntu 22.04, or higher versions of Linux distributions. This version cannot run directly on Windows systems，where modifications and adaptations are required.
* RMC 3.5.0 or higher, available upon resonable requests from the **'REAL'** group of the Department of Engineering Physics, Tsinghua University at http://reallab.ep.tsinghua.edu.cn.
* ANSYS Fluent 2021 R1 **(recommended)** or higher.
* hdf5, the version must be consistent with the RMC compilation dependencies and ANSYS Fluent’s hdf5 architecture, **v1.10.5 is recommended**.
* mpich, v3.2.1 or higher is recommended.
* cmake, v3.26.0 or higher is recommended.
### 1.3 Program Structure
* The `\lib_h5rw` folder contains the dynamic link libraries compiling external hdf5 read/write functions.
* The `\libudf` folder contains the UDF compilation path and library files for Fluent, which are `.dll` and `.lib` on **Windows**, and `.so` on **Linux**.
* `compile.sh` is the script to compile the above dynamic link libraries on Linux.
* `init.py` is the coupling initialization script.
* `iter.py` is the coupling iteration control script.
* `runFluent.jou` is the Fluent control script.
* **coupling.dat** is the **Coupling Input Card**.
  #### 1.3.1 Coupling Input Card Description
  * Line 1: **Maximum Picard iteration step** for steady-state or transient single time step.
  * Line 2: **Total power** for the steady-state or **current** transient time step, in MW.
  * Line 3: Photon or other direct heat generation contribution, usually set to 0.0.
  * Lines 4-6: Fluent geometry model's offset relative to RMC geometry, in cm. It is recommended that the geometry be aligned in the modeling stage, in which case these values should be 0.0.
  * Lines 7-10: Default values for fuel, moderator, reflector, and coolant temperatures in steady-state, followed by coolant density, all in SI units.
  * `#Project Name` is the name of the RMC input card and Fluent case file (without the extension).
  * `#Coupling mode` defines the calculation mode, 0 for steady-state and 1 for transient. 
  
    * *`Large time step transient`, i.e., burn-up calculation, is currently under development.*
  * `#time step` defines the time step settings, **important**: In steady-state, only the 5<sup>th</sup> parameter matters; for transient, all parameters are used:
    * The 1<sup>st</sup> value is the **transient neutronic-thermal coupling** time step (outer iteration time step), currently only fixed time steps are supported.
    * The 2<sup>nd</sup> value is the total duration for the **transient neutronic-thermal coupling**.
    * The 3<sup>rd</sup> value is the sub-time step for updating the neutron physics field's point-kinetics parameters, also known as the neutron inner iteration time step.
    * The 4<sup>th</sup> value is the time step for the **thermo-hydraulic field** transient equations, also known as the thermal-hydraulics inner iteration time step.
    * The 5<sup>th</sup> value is the maximum iteration for each time step of the thermal-hydraulics field steady-state or transient equations.
  * `#DIM` to `#OUTDIMR` define the divisions for the meshes that store power distribution, and fuel, coolant, moderator, and reflector thermal parameter distributions in the Cartesian coordinate directions (xyz). Currently, only uniform structured meshes are supported.
  * `#P_bndry` defines the boundaries for the mesh that stores the power distribution, **in cm**.
  * `#F_bndry` to `#R_bndry` define the boundaries for meshes that store fuel, coolant, moderator, and reflector thermal parameter distributions, **in m**.
  * `#Multilevel Flag` controls the multiscale model switch with 0 or 1, commonly used for models containing TRISO particles like HTGRs.
  * `#Omega` is the coupling relaxation factor, typically set to 0.5. **Currently, power is always relaxed**, to avoid relaxation, set this value to 1.
  * `#SOR flag for th_calculation` defines whether relaxation is applied to thermal-hydraulics data. **Current research shows little benefit, recommended to disable**.
  * `#Residual factor` is the upper α quantile for power, typically set to 3.0; for complex models, this value can be slightly relaxed.
  * `#tmax discrepancy tolerance` defines the maximum allowable relative temperature error.
  * `#Parallel processes` defines the number of parallel **processes** for both RMC and Fluent. These should not exceed the number of **cores** on your machine, **NOT** threads.
  * `#Output all iter flag` determines whether to output all iterations unconditionally, while the convergence will still be output.
  * `#Coupling order` defines the coupling sequence, **only valid for transient**: 0 means start with RMC, 1 means start with Fluent, default is 0.
---

## 2. Usage Guidelines

### 2.1 Steady-State Workflow

1. **Single-physics solution**. Prepare the RMC input card and a **Fluent case without UDF**.
    * **Neutron Physics**: Prepare the `.rmc` format input card.
      * Neutron generation information **must be set manually**, all Tally **must be set manually**, and MeshTally1 should have its mesh information on a separate line. Other Tally must **NOT** share the same line strategy as MeshTally1.
      * Use the Plot function to verify the geometry of the neutron physics model, paying attention to modifying the **geometry repetition strategy** for the grid elements that need to consider coupling effects.
      * Manually specify fixed grid element temperatures and check if the critical calculation works correctly.
      * Modify the thermohydraulic parameter transmission grid settings in the `coupling.dat` input card. Open a command line in the working directory and run:
        ```
        python init.py
        ```
        This generates an test thermohydraulic parameter distribution file in `.h5` format . Modify the input card according to RMC usage specifications, check if the external thermohydraulic parameter file can be read correctly, and perform a critical calculation. This step checks whether the model meets coupling conditions and whether memory overflow occurs. You can monitor memory usage with:
        ```
        watch free -g
        ```
      * After passing the tests, rename the `.rmc` input card. The name (without the extension) must match the project name in `coupling.dat`.
    * **Thermohydraulic Field**: Prepare a `.cas.h5` format case file.
      * Prepare the geometry and mesh model files.
      * Import the Fluent solver, do **NOT** load any UDFs, use **constant heat source**, and set boundary conditions to check if calculations can proceed. This verifies model correctness.
      * ***Optional***: If material properties are defined with UDFs, compile and load them to verify the property UDFs.
      * Perform grid independence study and turbulent wall function *y*^+^ verification. 
        * For the Standard *k*-ε model, it should satisfy 30 < *y*^+^ < 200
        * For the Standard *k*-ω model, it should satisfy *y*^+^ < 5.
      * After passing all tests, save the case. **It is recommended to save this case (without UDFs) separately for later use**.
2. **VERY IMPORTANT:** open the `.cas.h5` in GUI mode and record the **zone ID** of regions corresponding to different materials that need to consider coupling effects. Modify the corresponding lines in `udf.c`.
3. Open the `coupling.dat` input card, **check and modify the coupling calculation settings line by line**.
4. Open `iter.py` and **only modify** the file names.
5. Currently, the steady-state coupling working directory should contain all the interface files for the project, along with the `.rmc` and `.cas.h5` files. Running:
```
sh compile.sh
```
will generate the initial thermohydraulic parameter files such as `info_fuel.h5` and the dynamic link library `libh5rw.so`.
* **Manually load** the `libudf`, **save as a new file**, and **exit**. **The new file name must match the one in the coupling input card**.
6. Run:
```
python iter.py
```
---
### 2.2 Transient Workflow

1. **Perform steady-state coupling calculation first** (preferably in a separate directory to avoid confusion). Successful steady-state coupling means the model geometry, file reading/writing, memory usage, and UDF are all correct. For transient calculations, the **steady-state coupling *k*~eff~ value must converge within 1-2 standard deviations from 1**.
2. Copy the interface files from the steady-state coupling to the transient coupling working directory.
3. **Single-physics solution**. Prepare the RMC and Fluent input files.
   * **Neutron Physics**: Prepare `.rmc`, `.innerproduct`, and `.State.h5` format input cards.
     * When the *k*~eff~ condition is met, modify the steady-state `.rmc` input card to change the critical calculation to steady-state space-time dynamics (QUASISTATIC_S calculation), using the converged thermohydraulic parameter file. **Modify outer iteration time step and neutron number**, generating `.innerproduct` and `.State.h5` files.
     * Copy the QUASISTATIC_S input card to the transient working directory, change it to QUASISTATIC_D, and add inner iteration steps (this value can be arbitrarily set, the script will modify it based on `coupling.dat`), **but DO NOT change** neutron generation information or outer iteration settings. Also copy the `.innerproduct` and `.State.h5` files to the transient coupling working directory. Since `.State.h5` will be overwritten during transient calculations, **rename the `.State.h5` file**, for example, `Singlepin_init.State.h5`.
     * Open `coupling.dat` and modify the inner iteration time step count for QUASISTATIC_D to maintain a consistent inner iteration step size. For example, if the outer iteration time step is 0.5s or 1s, set the inner iteration count to 50 or 100, so the inner iteration time step is fixed as 0.01s.
  
        *`The selection of coupling outer iteration and neutron inner iteration time steps is still under research`*
   * **Thermohydraulic Field**: Prepare a `.cas.h5` format case file and a `.dat.h5` format initial condition file.
     * In the GUI, import the **steady-state case without any UDF** into the Fluent solver and load the converged steady-state CFD result `.dat.h5`. Use **constant heat source and switch to transient mode** to check if calculations can proceed and verify model correctness.
     * Perform Courant number validation for the inner iteration time step. According to the ANSYS Fluent user manual, the maximum Courant number should be between 20-40, and each inner iteration time step should converge within 10 iterations. Ideally, the Courant number should be close to 1.
     * After passing all tests, **compile and load UDFs, and save**. Since transient calculations are based on steady-state calculations, do not modify any data transmission grid settings, so no UDF modification is needed. Simply compile and load.
     * Copy the converged `.dat.h5` format CFD result file from steady-state coupling to the transient coupling working directory **to be used as the initial condition**. Since `.dat.h5` will be overwritten during transient calculations, **rename the `.dat.h5` file**, for example, `pin_init.dat.h5`.
4. Copy the converged power and thermohydraulic data transfer files from steady-state coupling, such as `MeshTally1.h5`, `info_fuel.h5`, to the transient coupling working directory.
5. Open the `coupling.dat` input card, **check and modify the coupling calculation settings line by line**.
6. Open `iter.py` and **only modify** the file names for archiving.
7. Currently, the transient coupling working directory should contain all interface files for the project, `.rmc` and `.cas.h5` files, and the eight files including initial conditions: `.innerproduct`, `_init.State.h5`, `MeshTally1.h5`, `_init.dat.h5`, `info_fuel.h5`, and others. Run:
```
python iter.py
```
