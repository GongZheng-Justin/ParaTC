module m_Variables
  use m_TypeDef
  use m_LogInfo
  use m_Parameters
  use m_decomp2d
  implicit none
  private
  
  ! define all major arrays here 
  real(RK), save,public,allocatable, dimension(:,:,:) :: ux
  real(RK), save,public,allocatable, dimension(:,:,:) :: uy
  real(RK), save,public,allocatable, dimension(:,:,:) :: uz
  real(RK), save,public,allocatable, dimension(:,:,:) :: HistxOld
  real(RK), save,public,allocatable, dimension(:,:,:) :: HistyOld
  real(RK), save,public,allocatable, dimension(:,:,:) :: HistzOld
  real(RK), save,public,allocatable, dimension(:,:,:) :: pressure
  
  real(RK), save,public,allocatable, dimension(:,:,:) :: RealArr1
  real(RK), save,public,allocatable, dimension(:,:,:) :: RealArr2
  real(RK), save,public,allocatable, dimension(:,:,:) :: Realhalo

  real(RK), save,public,allocatable, dimension(:,:) :: uy_ym, duy_ym

  type(MatBound),save,public:: mb1        ! matrix bound type 1
  type(HaloInfo),save,public:: hi1        ! halo info type 1
  type(HaloInfo),save,public:: hi_pr      ! halo info type for pressure
  type(HaloInfo),save,public:: hi_uxPrSrc ! halo info type for ux2
  type(HaloInfo),save,public:: hi_uzPrSrc ! halo info type for uz2

#ifdef ScalarFlow
  real(RK),save,public,allocatable,dimension(:,:,:):: scalar
  real(RK),save,public,allocatable,dimension(:,:,:):: HistCOld
  real(RK),save,public,allocatable,dimension(:,:,:):: RealArrC
#endif
  public:: AllocateVariables
contains  

  !******************************************************************
  ! InverseTridiagonal
  !****************************************************************** 
  subroutine AllocateVariables()
    implicit none
      
    ! locals
    integer::iErr01,iErr02,iErr03,iErr04,iErr05,iErr06,iErr07,iErrSum

    mb1%pencil = y_pencil  
    mb1%xme=2;  mb1%xpe=2
    mb1%yme=1;  mb1%ype=1
    mb1%zme=2;  mb1%zpe=2

    hi1%pencil = y_pencil
    hi1%xmh=2;  hi1%xph=2
    hi1%ymh=0;  hi1%yph=0
    hi1%zmh=2;  hi1%zph=2

    hi_pr%pencil = y_pencil
    hi_pr%xmh=2;  hi_pr%xph=1
    hi_pr%ymh=0;  hi_pr%yph=0
    hi_pr%zmh=2;  hi_pr%zph=1

    hi_uxPrSrc%pencil = y_pencil
    hi_uxPrSrc%xmh=1;  hi_uxPrSrc%xph=2
    hi_uxPrSrc%ymh=0;  hi_uxPrSrc%yph=0
    hi_uxPrSrc%zmh=0;  hi_uxPrSrc%zph=0

    hi_uzPrSrc%pencil = y_pencil
    hi_uzPrSrc%xmh=0;  hi_uzPrSrc%xph=0
    hi_uzPrSrc%ymh=0;  hi_uzPrSrc%yph=0
    hi_uzPrSrc%zmh=1;  hi_uzPrSrc%zph=2
       
    !-------------------------------------------------
    ! Arrays with ghost cells
    !-------------------------------------------------
    call myallocate(ux, mb1, opt_global=.true.)
    call myallocate(uy, mb1, opt_global=.true.)
    call myallocate(uz, mb1, opt_global=.true.)
    call myallocate(pressure, mb1, opt_global=.true.)
    call myallocate(RealHalo, mb1, opt_global=.true.)
#ifdef ScalarFlow
    call myallocate(scalar,   mb1, opt_global=.true.)
    call myallocate(RealArrC, mb1, opt_global=.true.)
    allocate(HistCOld(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)))
    scalar=Scalar_InitValue; RealArrC=zero; HistCOld=zero
#endif
    !-------------------------------------------------
    ! Arrays without ghost cells
    !-------------------------------------------------
    allocate(HistxOld(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)), Stat = iErr01)
    allocate(HistyOld(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)), Stat = iErr02)
    allocate(HistzOld(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)), Stat = iErr03)
    allocate(RealArr1(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)), Stat = iErr04)
    allocate(RealArr2(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)), Stat = iErr05)
    allocate(uy_ym(y1start(1):y1end(1),  y1start(3):y1end(3)), Stat = iErr06 )
    allocate(duy_ym(y1start(1):y1end(1), y1start(3):y1end(3)), Stat = iErr07 )
      
    iErrSum=abs(iErr01)+abs(iErr02)+abs(iErr03)+abs(iErr04)+abs(iErr05)+abs(iErr06)+abs(iErr07)
    if(iErrSum /= 0) call MainLog%CheckForError(ErrT_Abort,"AllocateVariables","Allocation failed")    
    
    ux=zero;         uy=zero;        uz=zero
    HistxOld=zero;   HistyOld=zero;  HistzOld=zero
    pressure=zero;   RealArr1=zero;  RealArr2=zero;  RealHalo=zero;
    uy_ym=zero;      duy_ym=zero
      
  end subroutine AllocateVariables
    
end module m_Variables
