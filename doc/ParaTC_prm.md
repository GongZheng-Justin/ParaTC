# input parameter for `ParaTC` 

&emsp;Input file for `ParaTC` can be divided into 8 part. Every part contains one **namelist**, i.e. `BasicParam`, `MeshSection`, `MeshOptions`, `SpectraOptions`, `Spectra2DOpt`, `IO_Options`, `SaveScalarOption`, and `ScalarFlowOptions`.

&emsp;**NOTE:**
* The parameter sequence in a **namelist** is not important.
* A **namelist** is started with `&NamelistName` (e.g. `&BasicParam`), and ended with `/End`.
* You can also add blank lines or comments (after sign `!`) in the input file, just like what you do in a free format Fortran code.

&emsp;Consider the following input file:

```
!===================
&BasicParam
!===================

  RestartFlag = T         ! restart or not

  ! Flow type ( 1=Channel, 2=Half channel )
  FlowType  = 2
  ubulk= 1.4352
  IsUxConst = F
  IsUseCRF  = F          ! use Converting Reference Frame or not
  uCRF= 0.0              ! Velocity of the Converting Reference Frame

  ! Mesh options
  xlx = 12.56637061      ! domain length in x-dir
  yly = 1.0              ! domain length in y-dir 
  zlz = 4.188790205      ! domain length in z-dir
  nxp =  193             ! grid point number in x-dir
  nyp =   97             ! grid point number in y-dir
  nzp =  129             ! grid point number in z-dir

  ! Physical properties
  xnu = 4.66591665597541E-04   ! kinematic viscosity
  gravity =0.0 0.0 0.0   ! Gravity or other constant body forces (if any)
  FluidDensity = 1000    ! fluid density

  ! Time stepping
  dtMax= 0.015           ! Maxium time step
  iCFL = 1               ! Use CFL condition to change time step dynamically( 1: yes, 2:no ).
  CFLc = 1.5             ! CFL parameter
  ifirst= 1              ! First iteration
  ilast = 40000          ! Last iteration

  ! Numerical scheme options
  ischeme = 3                 ! (1=AB2, 2=RK2, 3=RK3)
  FFTW_plan_type = 1          ! (1=FFTW_MEASURE,  2=FFTW_ESTIMATE)

  ! Boundary conditions
  !                    x-,  x+,  y-,  y+,  z-,  z+
  ! And the velocity Bc values will ONLY be used corresponding to the no-slip Bc
  ! While at the same time, you are allow to assign a  transpiration Bc for uy in y- bottom(uyBcValue(3) dosen't work).
  uxBcValue = 0.0, 0.0,  0.0,  0.0,  0.0,  0.0
  uyBcValue = 0.0, 0.0,  0.0,  0.0,  0.0,  0.0    
  uzBcValue = 0.0, 0.0,  0.0,  0.0,  0.0,  0.0 

  ! I/O, Statistics
  ivstats  = 20                  ! time step interval for statistics calculation
  saveStat = 5000                ! Output Statistics file frequency
  SaveVisu = 5000                ! Output visualizing file frequency
  BackupFreq =10000              ! Output Restarting file frequency
  RunName  = "TurbidityCurrent"  ! Run name
  Res_Dir  = "./CFD/Results/"    ! Result directory
  RestartDir = "./CFD/Restart/"  ! Restart directory
  Cmd_LFile_Freq  = 5            ! Report frequency in the terminal
  LF_file_lvl     = 5            ! Logfile report level
  LF_cmdw_lvl     = 3            ! Terminal report level

  ! Decomp2d options
  p_row = 2 !7
  p_col = 4 !7

  ! limited velocity and div
  vel_limit = 4.0
  div_limit = 0.2

/End of NAMELIST "&BasicParam"

!=================
&MeshSection
!=================

  nSection      =   2       ! yly will be diveded into "nSection" part
  
/End of NAMELIST "MeshSection"

!=================
&MeshOptions
!=================

  SectionLength =   1     1 ! e.g. If nSection=2, and SectionLength=[1,3], yly is further divided into 1/4*yly and 3/4*yly
  nycSection    =  48    48 ! sum(nycSection)=nyc
  StretType     =   2     2 ! 0:Uniform; 1:Tangent hyperbolic function; 2:Sine/cosine function; 3:Proportional sequence
  StretOption   =   0     1 ! 0:bottom;  1:top
  SectioncStret =  1.0  1.0 ! Stretching parameter. if StretType=0, this parameter doesn't work.
  
/End of NAMELIST "MeshOptions"

!=================
&SpectraOptions
!=================

  clcSpectra1D =  T
  ivSpec   = 20
  jForLCS  = 10

  clcSpectra2D =  F
  nySpec2D = 9
  
/End of NAMELIST "SpectraOptions"

!=================
&Spectra2DOpt
!=================

  iySpec2D = 5 10 15 20 25 30 35 40 45

/End of NAMELIST "Spectra2DOpt"

!=================
&IO_Options
!=================
  
  iskip   = 1
  jskip   = 1
  kskip   = 1
  XDMF_SET_TYPE=1  ! 0: cell, 1:Node

  save_ux    = T
  save_uy    = T
  save_uz    = T
  save_wx    = F
  save_wy    = F
  save_wz    = F
  save_wMag  = F
  save_pr    = F
  save_Q_vor = F
  save_lamda2= F

  WriteHistOld = T
  ReadHistOld  = T

/End of NAMELIST "&IO_Options"

!=================
&SaveScalarOption
!=================

  save_scalar    = T

/End of NAMELIST "&SaveScalarOption"

!=================
&ScalarFlowOptions
!=================

  FallingVel= 1.679729996151150E-03   ! Particle Settling Velocity.
  SchmidtNumber= 1.0            ! Schmidt number for gravity flow (Or called Prandtl number sometimes). Sc= nu/K, where K is the diffusivity coefficient.
  GravityEff= 16.1865           ! Effective Gravity magnitude. g_eff=Rg=(rho_s-rho_f)/rho_f*g for gravity flow.
  GravityDir= 8.71557427476582E-02 -9.961946980917460E-01  0.0 ! Unit gravity vector.                                  

  ! Scalar boundary conditions
  ! -1: Dirichlet Bc. C     = F
  ! -2: Neumann Bc.   dC/dy = S
  ! -3: Robin Bc.     dC/dy + F*C = S
  ! y-,  y+
  ScalarBcOption  =  -3,    -3
  ScalarBcValues  =   3.58630091313028 3.58630091313028   0.0,  0.0 !(F_b,F_t,S_b,S_t), 4.66219118706937

  Scalar_InitValue= 5.0E-3 ! Initially, a uniform scalar field is set to the value of Scalar_InitValue.
  
/End of NAMELIST "ScalarFlowOptions"
```
&emsp;As you can see, most parameters are followed by some comments to illustrate its function.

## BasicParam
&emsp;**BasicParam** specifies several basic parameters, e.g. input/output options, physical properties.

* `RestartFlag`: logical type. Restart or not. If RestartFlag=.true., the simulation will start from a prestored restarting file.
* `FlowType`: integer type. Specify the flow type, 1=Channel, 2=Half channel.
* `ubulk`: real type. Mean bulk velocity in x-direction. This parameter is only used for wall-bounded turbulent channel flows (channel or half channel) to set the initial mean streamwise velocity.
* `IsUxConst`:  logical type. Whether ux is constant or not. This parameter is only used for wall-bounded turbulent channel flows (channel or half channel). If IsUxConst=.true., the turbulent channel flow is driven by a uniform pressure gradient, which varies in time to maintain a constant mass flow rate in streamwise (x) direction. If IsUxConst=.false., the uniform pressure gradient will be also constant in time, which is set by parameter `gravity`.
* `IsUseCRF`: logical type. Whether use Converting Reference Frame or not. If so, the computations are performed in a moving reference frame for which the bulk velocity (net streamwise mass flux) is zero.
* `uCRF`: Velocity of the Converting Reference Frame. If `IsUseCRF=F`, this parameter does not work.
* `xlx`: real type. Domain length in x-direction.
* `yly`: real type. Domain length in y-direction.
* `zlz`: real type. Domain length in z-direction.
* `nxp`: integer type. Grid point number in x-direction. **Note:** If nxp=129, there will be 129 grid points in x-direction, and the x-domain will be divided into 128 parts.
* `nyp`: integer type. Grid point number in y-direction.
* `nzp`: integer type. Grid point number in z-direction.
* `xnu`: real type. Kinematic viscosity.
* `gravity`: real vector containing 3 components. Gravity or other constant body forces (if any), which can be used to drive the flow.
* `FluidDensity`: fluid density
* `dtMax`: real type. Maxium allowable time step. If iCFL=1, `dt` will adjust according to the allowable CFL number. `dt=min{dtMax, CFLc/max(|u|/dx+|v|/dy+|w|/dz)}`. If iCFL=2, `dt=dtMax`. 
* `iCFL`: integer type. Whether use CFL condition to change time step dynamically or not. 1:yes, 2:no.
* `CFLc`: real type. Allowable CFL parameter. If iCFL=2, this parameter will not work.
* `ifirst`: integer type. First iteration.
* `ilast`:  integer type. Last iteration.
* `ischeme`: integer type. Specify the time integral scheme. 1=AB2, 2=RK2, 3=RK3.
* `FFTW_plan_type`: integer type. 1=FFTW_MEASURE, 2=FFTW_ESTIMATE. **Note:** In my practice, using *FFTW_MEASURE* is faster than *FFTW_ESTIMATE*, while *FFTW_MEASURE* might lead to different *FFTW_plan*, even for the same simulation case. And further, different *FFTW_plan* will result in slight different DFT values (The last several digits of the decimal point might be different). If you want to debug the code, please use *FFTW_ESTIMATE*, which will lead to the same DFT values for a specific case.
* `uxBcValue`: real vector containing 6 components to specify the u-velocity Bc values, and these values will ONLY be used corresponding to the no-slip Bc.
* `uyBcValue`: real vector containing 6 components to specify the v-velocity Bc values.
* `uzBcValue`: real vector containing 6 components to specify the w-velocity Bc values.
* `ivstats`: integer type. Time step interval for statistics calculation. The statistic value will be calculated every `ivstats` time step.
* `SaveStat`: integer type. Output Statistics file frequency. The statistic value will be written to the folder `Res_Dir` every `SaveStat` time step.
* `SaveVisu`: integer type. Output visualizing file frequency. The visualizing data will be written to the folder `Res_Dir` every `SaveVisu` time step.
* `BackupFreq`: integer type. Output restarting file frequency. The restarting file will be be written to the folder `RestartDir` every `BackupFreq` time step, and the previous restarting file will be deleted in order to save memory room.
* `RunName`: character type. Run name string.
* `Res_Dir`: character type. The folder to store result data.
* `RestartDir`: character type. The folder to store restart data.
* `Cmd_LFile_Freq`: integer type. Reporting frequency. Reporting information  will be printed to the terminal and written to the logfile every `Cmd_LFile_Freq` time step.
* `LF_file_lvl`: integer type. Logfile report level. From 1 to 5. `5` means every message will be reported into logfile, smaller `LF_file_lvl` corresponds to less reported message.
* `LF_cmdw_lvl`: integer type. Terminal report level. From 1 to 5.
* `p_row`: integer type. Number of processors in row pencil.
* `p_col`: integer type. Number of processors in column pencil.
* `vel_limit`: real type. Maximum allowable velocity. If `min{|u|,|v|,|w|}>vel_limit`, the program will abort.
* `div_limit`: real type. Maximum allowable divergence. If `max{abs(du/dx +dv/dy + dw/dz)}>div_limit`, the program will abort.

## MeshSection
* `nSection`: integer type. yly will be diveded into "nSection" part

## MeshOptions
&emsp;**MeshOptions** specifies several parameters for y-mesh.

* `SectionLength`: real vector type containing `nSection` components. relative length for every section. If nSection=2, and SectionLength=1,3, yly is further divided into 1/4*yly and 3/4*yly.
* `nycSection`: integer vector type containing `nSection` components. y mesh number for every section. sum(nycSection)=nyc.
* `StretType`: integer vector type containing `nSection` components. Whether y-mesh is stretched or not. 0:Uniform; 1:Tangent hyperbolic function; 2:Sine/cosine function; 3:Proportional sequence.
* `StretOption`: integer vector type containing `nSection` components. Which side is streatched. 0:bottom;  1:top.
* `SectioncStret`: integer vector type containing `nSection` components. Stretching parameter. if StretType=0, this parameter doesn't work.

## SpectraOptions
&emsp;**SpectraOptions** designates options for energy spectra calculation.

* `clcSpectra1D`: logical type. Calculate 1D energy spectra or not.
*  `ivSpec`: integer type. Time step interval for energy spectra calculation. The energy spectr will be calculated every `ivSpec` time step.
* `jForLCS`: integer type. Reference y-index for Linear coherence spectrum.
* `clcSpectra2D`: logical type. Calculate 2D energy spectra or not.
* `nySpec2D`: integer type. Number of y-locations for 2D energy spectra.

## Spectra2DOpt
* `iySpec2D`: integer vector type containing `nySpec2D` components. index of y-locations for 2D energy spectra.

## IO_Options
&emsp;**IO_Options** sets input/output options.

* `iskip`: integer type. Flow data will be written into output visualizing file every `iskip` grid mesh in x-dir.
* `jskip`: integer type. Flow data will be written into output visualizing file every `jskip` grid mesh in y-dir.
* `kskip`: integer type. Flow data will be written into output visualizing file every `kskip` grid mesh in z-dir.
* `XDMF_SET_TYPE`: integer type. The output XDMF_SET_TYPE, 0: cell, 1:Node. 
* `save_ux`: logical type. Save u-velocity or not.
* `save_uy`: logical type. Save v-velocity or not.
* `save_uz`: logical type. Save w-velocity or not.
* `save_wx`: logical type. Save x-vorticity or not.
* `save_wy`: logical type. Save y-vorticity or not.
* `save_wz`: logical type. Save z-vorticity or not.
* `save_wMag`: logical type. Save vorticity magnitude or not.
* `save_pr`: logical type, Save pressure or not.
* `save_Q_vor`: logical type. Save Q vortex criterion value or not.
* `save_lamda2`: logical type. Save lamda2 vortex criterion value or not.
* `WriteHistOld`: logical type. Write the old convective data (more exactly, the terms which are treated explicitly) to restart value or not.
* `ReadHistOld`: logical type. Read the old convective data (more exactly, the terms which are treated explicitly) to restart value or not. This parameter, coupled with the previous one `WriteHistOld`, can help us to make the program restarting more flexible. Sometimes, we want to restart a `AB2` case from `RK3` case. `AB2` integral approach need the old convective data, while `RK3` method not. So, if we want to restart a `AB2` case from `RK3` case, we can set `ReadHistOld=F` for `AB2` case. If we want to restart a `AB2` case from a previous `AB2` case, we can set `ReadHistOld=T` for the new `AB2` case, and `WriteHistOld=T` for the previous case. If we want to restart a `RK3` case from `AB2` case, we can set `WriteHistOld=F` for the previous `AB2` case.

## SaveScalarOption
* `save_scalar`: logical type. Save save_scalar field or not.

## ScalarFlowOptions
&emsp;**ScalarFlowOptions** sets B.C. options for scalar field.

`FallingVel`: real type. Particle Settling Velocity.  
`SchmidtNumber`: real type. Schmidt number. Schmidt number for gravity flow (Or called Prandtl number sometimes). Sc= nu/K, where K is the diffusivity coefficient. 
`GravityEff`: real type. Effective Gravity magnitude. g_eff=Rg=(rho_s-rho_f)/rho_f* g for gravity flow. 
`GravityDir`: real vector type containing 3 components. Unit gravity vector pointing at gravity direction. 
`ScalarBcOption`: integer vector type containing 2 components. Scalar boundary conditions for bottom and top walls, respectively.  -1: Dirichlet Bc. C = F; -2: Neumann Bc. dC/dy = S; -3: Robin Bc. dC/dy + F* C = S 
`ScalarBcValues`: real vector type containing 4 components. Assign F_b,F_t,S_b,S_t. 
`Scalar_InitValue`: real type. Initially, a uniform scalar field is set to the value of Scalar_InitValue. 
