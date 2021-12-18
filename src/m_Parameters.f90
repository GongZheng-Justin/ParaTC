module m_Parameters
  use m_TypeDef
  use m_LogInfo
  implicit none  
  private
  
  ! Flow type option
  integer,parameter,public:: FT_CH=1  ! Channel
  integer,parameter,public:: FT_HC=2  ! Half channel
  integer, public,save:: FlowType
  real(RK),public,save:: ubulk
  logical,public,save :: IsUxConst
  logical,public,save:: IsUseCRF
  real(RK),public::PrGradAver,uCRF
  
  ! mesh options
  real(RK),public,save::xlx, yly, zlz  ! domain length in three directions
  integer,public,save::nxp,nyp,nzp     ! grid point number in three directions
  integer,public,save::nxc,nyc,nzc     ! grid center number in three directions. nxc = nxp-1 
  
  integer,parameter,public:: x_pencil=1
  integer,parameter,public:: y_pencil=2
  integer,parameter,public:: z_pencil=3
  integer,parameter,public:: xm_dir=1  ! x- direction
  integer,parameter,public:: xp_dir=2  ! x+ direction
  integer,parameter,public:: ym_dir=3  ! y- direction
  integer,parameter,public:: yp_dir=4  ! y+ direction
  integer,parameter,public:: zm_dir=5  ! z- direction
  integer,parameter,public:: zp_dir=6  ! z+ direction  
  
  ! Physical properties
  real(RK),save,public:: xnu                    ! Kinematic viscosity
  real(RK),save,public:: FluidDensity           ! Fluid density 
  real(RK),dimension(3),save,public:: gravity   ! Gravity or  other constant body forces (if any)
  
  ! Time stepping scheme and Projection method options
  integer,parameter,public:: FI_AB2 =1          !  AB2 for convective term
  integer,parameter,public:: FI_RK2 =2          !  RK2 for convective term
  integer,parameter,public:: FI_RK3 =3          !  RK3 for convective term
  real(RK),save,public:: dt         ! current time step
  real(RK),save,public:: dtMax      ! Maxium time step
  real(RK),save,public:: SimTime    ! Real simulation time
  integer,save,public :: iCFL       ! Use CFL condition to change time step dynamically( 1: yes, 2:no ).
  real(RK),save,public:: CFLc       ! CFL parameter
  integer,save,public::  itime      ! current time step
  integer,save,public::  ifirst     ! First iteration
  integer,save,public::  ilast      ! Last iteration 
  integer,save,public::  iadvance
  integer,save,public::  ischeme
  real(RK),dimension(3),save,public :: pmGammaConst
  real(RK),dimension(3),save,public :: pmThetaConst
  real(RK),dimension(3),save,public :: pmAlphaConst
  real(RK),save,public::pmGamma
  real(RK),save,public::pmTheta
  real(RK),save,public::pmAlpha
  real(RK),save,public::pmAlphaC
  real(RK),save,public::pmBeta
  real(RK),save,public::pmBetaT
#ifdef ScalarFlow
  real(RK),save,public::pmBetaSc
#endif

  ! FFTW_option
  integer,save,public:: FFTW_plan_type
  
  ! Boundary conditions
  real(RK),public,dimension(6):: uxBcValue
  real(RK),public,dimension(6):: uyBcValue
  real(RK),public,dimension(6):: uzBcValue

  ! I/O, Statistics
  integer,public,save:: ivstats            ! time step interval for statistics calculation 
  integer,public,save:: SaveVisu           ! Output visulizing file frequency
  integer,public,save:: BackupFreq         ! Output Restarting file frequency
  integer,public,save:: SaveStat           ! Output Statistics file frequency

  logical,public,save::       RestartFlag  ! restart or not
  character(64),public,save:: RunName      ! Run name
  character(64),public,save:: Res_Dir      ! Result directory
  character(64),public,save:: RestartDir   ! Restart directory
  integer,public,save:: Cmd_LFile_Freq= 1  ! report frequency in the terminal 
  integer,public,save:: LF_file_lvl   = 5  ! logfile report level      
  integer,public,save:: LF_cmdw_lvl   = 3  ! terminal report level

  ! Decomp2d options
  integer,public,save:: p_row,p_col

  ! limited velocity and div
  real(RK),public,save:: vel_limit
  real(RK),public,save:: div_limit

#ifdef ScalarFlow
  !=== Scalar Options ===!
  ! Scalar boundary conditions
  ! -1: Dirichlet Bc. C     = F
  ! -2: Neumann Bc.   dC/dy = S
  ! -3: Robin Bc.     dC/dy + S*C = F
  ! y-,  y+
  real(RK),public,save:: SDiffCoe ! Scalar Diffusivity Coefficient
  real(RK),public,save:: Scalar_InitValue
  real(RK),public,dimension(3):: GravityEffVec,FallingVelVec
  integer, public,dimension(2):: ScalarBcOption
  real(RK),public,dimension(4):: ScalarBcValues !(F_b,F_t,S_b,S_t)
#endif
  
  public:: ReadAndInitParameters,DumpReadedParam,PMcoeUpdate
contains
    
  !******************************************************************
  ! InitParameters
  !****************************************************************** 
  subroutine ReadAndInitParameters(chFile)
    implicit none 
    character(*),intent(in)::chFile
    
    ! locals
    integer:: nUnitFile,myistat
    NAMELIST/BasicParam/FlowType,ubulk,xlx,yly,zlz,nxp,nyp,nzp,xnu,dtMax,iCFL,CFLc,ifirst,ilast,  &
                        ischeme,FFTW_plan_type,gravity,uxBcValue,uyBcValue,uzBcValue,ivstats,     &
                        BackupFreq,SaveStat,SaveVisu,IsUxConst,IsUseCRF,uCRF,RestartFlag,         &
                        RunName,Res_Dir,RestartDir,Cmd_LFile_Freq,LF_file_lvl,LF_cmdw_lvl,p_row,  &
                        p_col,vel_limit,div_limit,FluidDensity
#ifdef ScalarFlow
    real(RK)::RealTemp
    real(RK)::FallingVel    !Particle Settling Velocity.
    real(RK)::GravityEff    !Effective Gravity magnitude. g_eff=Rg=(rho_s-rho_f)/rho_f*g for gravity flow.
    real(RK)::SchmidtNumber !Schmidt number for gravity flow (Or Prandtl number somewhere). Sc= nu/K, where K is the diffusivity coefficient
    real(RK),dimension(3)::GravityDir ! Unit gravity vector
    NAMELIST/ScalarFlowOptions/SchmidtNumber,Scalar_InitValue,GravityEff,GravityDir,FallingVel,ScalarBcOption,ScalarBcValues
#endif

    nUnitFile = GetFileUnit()
    open(unit=nUnitFile, file=chFile, status='old',form='formatted',IOSTAT=myistat )
    if(myistat/=0) then
      print*,"Cannot open file: "//trim(chFile); STOP
    endif
    read(nUnitFile, nml=BasicParam)
#ifdef ScalarFlow
    rewind(nUnitFile)
    read(nUnitFile, nml=ScalarFlowOptions)
#endif
    close(nUnitFile,IOSTAT=myistat)  

    nxc= nxp-1
    nyc= nyp-1
    nzc= nzp-1
    if(mod(nxc,2)/=0   .or.  mod(nzc,2)/=0) then
      print*,'mod(nxc,2)/=0 || mod(nzc,2)/=0 WRONG';STOP
    endif
    if(FlowType==FT_CH .and. mod(nyc,2)/=0) then
      print*,'mod(nyc,2)/=0 For channel geomotry. WRONG';STOP
    endif
    if(IsUseCRF) then
      uxBcValue= uxBcValue -uCRF
    else
      uCRF= zero
    endif
  
    if(ischeme==FI_AB2) then
      iadvance = 1
      pmGammaConst = (/ three/two, zero, zero/)
      pmThetaConst = (/     -half, zero, zero/)
    elseif(ischeme==FI_RK2) then
      iadvance = 2
      pmGammaConst = (/ half,  one, zero/)
      pmThetaConst = (/ zero,-half, zero/)
    elseif(ischeme==FI_RK3) then
      iadvance = 3
      pmGammaConst = (/ eight/fifteen,       five/twelve,    three/four/)
      pmThetaConst = (/          zero,  -seventeen/sixty,  -five/twelve/)
    else
      print*,"Time scheme WRONG!!! ischeme=",ischeme; STOP      
    endif
    pmAlphaConst = pmGammaConst + pmThetaConst

#ifdef ScalarFlow
    RealTemp=sqrt(GravityDir(1)*GravityDir(1)+GravityDir(2)*GravityDir(2)+GravityDir(3)*GravityDir(3))
    if(RealTemp>1.0E-10_RK) then
      GravityDir=GravityDir/RealTemp
    else
      GravityDir=[zero, -one, zero]
    endif
    SDiffCoe=xnu/SchmidtNumber
    GravityEffVec=GravityDir*GravityEff
    FallingVelVec=GravityDir*FallingVel
    FallingVelVec(1)=FallingVelVec(1)+uCRF
#endif
  end subroutine ReadAndInitParameters

  !******************************************************************
  ! DumpReadedParam
  !****************************************************************** 
  subroutine DumpReadedParam()
    implicit none

    ! locals
    NAMELIST/BasicParam/FlowType,ubulk,xlx,yly,zlz,nxp,nyp,nzp,xnu,dtMax,iCFL,CFLc,ifirst,ilast,  &
                        ischeme,FFTW_plan_type,gravity,uxBcValue,uyBcValue,uzBcValue,ivstats,     &
                        BackupFreq,SaveStat,SaveVisu,IsUxConst,IsUseCRF,uCRF,RestartFlag,         &
                        RunName,Res_Dir,RestartDir,Cmd_LFile_Freq,LF_file_lvl,LF_cmdw_lvl,p_row,  &
                        p_col,vel_limit,div_limit,FluidDensity
#ifdef ScalarFlow
    NAMELIST/ScalarFlowOptions/SDiffCoe,Scalar_InitValue,GravityEffVec,FallingVelVec,ScalarBcOption,ScalarBcValues
#endif

    write(MainLog%nUnit, nml=BasicParam)
#ifdef ScalarFlow
    write(MainLog%nUnit, nml=ScalarFlowOptions)
#endif    
  end subroutine DumpReadedParam

  !******************************************************************
  ! PMcoeUpdate
  !******************************************************************
  subroutine PMcoeUpdate(ns)
    implicit none
    integer,intent(in)::ns
      
    pmGamma = pmGammaConst(ns) * dt
    pmTheta = pmThetaConst(ns) * dt
    pmAlpha = pmAlphaConst(ns) * dt
    pmBeta  = half * pmAlphaConst(ns) *xnu * dt
#ifdef ScalarFlow
    pmBetaSc= half * pmAlphaConst(ns) *SDiffCoe * dt
#endif
    SimTime=SimTime+dt*real(pmAlphaConst(ns),RK)

    pmAlphaC= pmAlphaConst(ns)              ! for pressure gradient purpose only
    if(ns==1) PrGradAver=zero
  end subroutine PMcoeUpdate
  
end module m_Parameters
