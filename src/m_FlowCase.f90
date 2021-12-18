module m_FlowCase
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_Parameters
  use m_MeshAndMetries
  use m_decomp2d
  use m_Variables,only:mb1
  use m_Tools,only:CalcUxAver
  use iso_c_binding
  implicit none
  private
  include "fftw3.f03"

  ! statistics variabls
  integer,save:: nfstime
  real(RK),save:: PrGradsum
  real(RK),save,allocatable,dimension(:,:):: SumStat

  ! SpectraOptions
  integer,save:: ivSpec,jForLCS,nySpec2D
  logical,save:: clcSpectra1D,clcSpectra2D
  integer,allocatable,dimension(:)::iySpec2D
  
  ! Spectra variables
  type(C_PTR),save::fft_plan_x,fft_plan_z,fft_plan_z2
  integer,save:: nxh,nxhp,nzh,nzhp,nSpectime  
  real(RK),save,allocatable,dimension(:,:,:)::EnergySpecX,EnergySpecZ
  real(RK),save,allocatable,dimension(:,:,:,:)::EnergySpec2D
  
  public:: InitVelocity, Update_uy_ym, InitStatVar, clcStat
contains
!#define SAVE_SINGLE_Spec2D

#if defined(ScalarFlow)
#define NCHASTAT 40
#define NEnergySpec1D 9
#define NEnergySpec2D 6

#else 
#define NCHASTAT 35
#define NEnergySpec1D 8
#define NEnergySpec2D 5

#endif
  !******************************************************************
  ! InitVelocity
  !******************************************************************
  subroutine InitVelocity(ux,uy,uz,Deviation)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)),intent(inout)::Deviation
  
    ! locals
    integer :: ii,code,i,j,k,m1,m2
    real(RK):: retau_guass,utau_guass,height,rem,wx,wz,xlxPlus,zlzPlus
    real(RK):: xplus,yplus,zplus,yct,ybar,xp,zp,ratiot,uzmean(nyc),uzmeanR(nyc)

    height=yly
    if(FlowType==FT_CH)height=half*yly

    ux=zero; uy=zero; uz=zero
    rem = ubulk * height / xnu
    retau_guass = 0.1538_RK*rem**0.887741_RK
    utau_guass   = retau_guass*xnu/height
    if(nrank==0)print*,'************** retau_gauss=',retau_guass
    if(nrank==0)print*,'************** utau_gauss= ',utau_guass

    call system_clock(count=code); !code=0
    call random_seed(size = ii)
    call random_seed(put = code+63946*(/ (i - 1, i = 1, ii) /))
    call random_number(Deviation)
    Deviation= 0.2_RK* Deviation + 0.9_RK ! [0.8, 1.2]

    !modulation of the random noise + initial velocity profile
    uzmean=zero
    wx=twopi/500.0_RK; wz=twopi/200.0_RK
    xlxPlus=xlx*utau_guass/xnu;   zlzPlus=zlz*utau_guass/xnu;
    m1=floor(xlxPlus*wx/twopi)+1; wx=real(m1,RK)*twopi/xlxPlus
    m2=floor(zlzPlus*wz/twopi)+1; wz=real(m2,RK)*twopi/zlzPlus
    do j=y1start(2),y1end(2)
      yct = height-abs(height-yc(j)) 
      ybar= yct/height; yplus=utau_guass*yct/xnu
      do k=y1start(3),y1end(3)
        zp   =real(k-1,kind=RK)*dz+dz*half
        zplus=utau_guass*zp/xnu
        do i=y1start(1),y1end(1)
          xp   =real(i-1,kind=RK)*dx+dx*half
          xplus=utau_guass*xp/xnu
          !ux(i,j,k) = 0.0052_RK*ubulk*yplus*exp(-yplus*yplus/1800.0_RK)*cos(wz*zplus)*Deviation(i,j,k) ! original expression
          !uz(i,j,k) = 0.0050_RK*ubulk*yplus*exp(-yplus*yplus/1800.0_RK)*sin(wx*xplus)*Deviation(i,j,k) ! original expression
          ux(i,j,k) = ubulk*ybar*exp(-4.5_RK*ybar*ybar)*cos(wz*zplus)*Deviation(i,j,k)
          uz(i,j,k) = ubulk*ybar*exp(-4.5_RK*ybar*ybar)*sin(wx*xplus)*Deviation(i,j,k)
          ux(i,j,k) = ux(i,j,k)+ three*ubulk*(ybar-half*ybar*ybar)
          uzmean(j) = uzmean(j)+ uz(i,j,k)
        enddo
      enddo
      uzmean(j)=uzmean(j)/real(nxc*nzc,RK)
    enddo
    ratiot=ubulk/CalcUxAver(ux)
    call MPI_Bcast(ratiot,1,real_type,0,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(uzmean,uzmeanR,nyc,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    ux= ux*ratiot -uCRF
    do j=y1start(2),y1end(2)
      uz(:,j,:)=uz(:,j,:)-uzmeanR(j)
    enddo   
  end subroutine InitVelocity

  !******************************************************************
  ! Update_uy_ym
  !******************************************************************   
  subroutine Update_uy_ym(uy_ym, duy_ym)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(3):y1end(3)),intent(out):: uy_ym
    real(RK),dimension(y1start(1):y1end(1),y1start(3):y1end(3)),intent(inout):: duy_ym
  
    duy_ym = uy_ym
     
    ! update uy_ym here
    uy_ym = zero
    duy_ym= uy_ym - duy_ym
  end subroutine Update_uy_ym

  !******************************************************************
  ! InitStatVar
  !******************************************************************
  subroutine InitStatVar(ChannelPrm)
    implicit none 
    character(*),intent(in)::ChannelPrm

    ! locals
    character(len=50)::filename
    integer::ierror,nUnit,plan_type,jc
#ifdef OverWriteFFT
    real(RK),dimension(:,:,:),allocatable::Arr
#else
    real(RK),dimension(:,:,:),allocatable::Arr1,Arr2
#endif
    type(fftw_iodim),dimension(1)::iodim,iodim_howmany
    NAMELIST/SpectraOptions/clcSpectra1D,ivSpec,jForLCS,clcSpectra2D,nySpec2D
    NAMELIST/Spectra2DOpt/iySpec2D
    
    nUnit = GetFileUnit() 
    open(unit=nUnit, file=ChannelPrm, status='old',form='formatted',IOSTAT=ierror )
    if(ierror/=0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar", "Cannot open file: "//trim(ChannelPrm))
    read(nUnit, nml=SpectraOptions)
    close(nUnit,IOSTAT=ierror)
    
    if(nrank==0) then
      write(MainLog%nUnit, nml=SpectraOptions)
      if(mod(saveStat,ivstats)/=0 )    call MainLog%CheckForError(ErrT_Abort,"InitStatVar","ivstats wrong !!!")
      if(clcSpectra1D .and. (mod(saveStat,ivSpec)/=0 .or. mod(ivSpec,ivstats)/=0 )) then
        call MainLog%CheckForError(ErrT_Abort,"InitCAStatistics","ivSpec wrong !!!")
      endif
      if(IsUxConst)then
        nUnit=GetFileUnit()
        write(filename,'(A,I10.10)')trim(Res_dir)//"PrGrad",ilast
        open(nUnit,file=filename,status='replace',form='formatted',IOSTAT=ierror)
        if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitCAStatistics","Cannot open file: "//trim(filename))
        close(nUnit,IOSTAT=ierror)
      endif
    endif
    allocate(SumStat(NCHASTAT,nyp),Stat=ierror)
    if(ierror /= 0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar: ","Allocation failed")  
    nfstime=0;  nSpectime=0; SumStat=zero; PrGradsum=zero

    if(FFTW_plan_type == 1) then
      plan_type=FFTW_MEASURE
    else
      plan_type=FFTW_ESTIMATE
    endif    
    nxh=nxc/2; nxhp=nxh+1
    nzh=nzc/2; nzhp=nzh+1
    IF(clcSpectra1D) THEN
      iodim(1)%n  = x1size(1)
      iodim(1)%is = 1
      iodim(1)%os = 1
      iodim_howmany(1)%n  = x1size(2)*x1size(3)
      iodim_howmany(1)%is = x1size(1)
      iodim_howmany(1)%os = x1size(1)
      allocate(EnergySpecX(nxhp,x1size(2),NEnergySpec1D));EnergySpecX=zero
#ifdef OverWriteFFT
      allocate(Arr(x1size(1),x1size(2),x1size(3)))
      fft_plan_x = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,Arr,Arr,[FFTW_R2HC],plan_type)
      deallocate(Arr)
#else
      allocate(Arr1(x1size(1),x1size(2),x1size(3)),Arr2(x1size(1),x1size(2),x1size(3)))
      fft_plan_x = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,Arr1,Arr2,[FFTW_R2HC],plan_type)
      deallocate(Arr1,Arr2)
#endif

      iodim(1)%n  = z1size(3)
      iodim(1)%is = z1size(1)*z1size(2)
      iodim(1)%os = z1size(1)*z1size(2)
      iodim_howmany(1)%n  = z1size(1)*z1size(2)
      iodim_howmany(1)%is = 1
      iodim_howmany(1)%os = 1    
      allocate(EnergySpecZ(nzhp,z1size(2),NEnergySpec1D));EnergySpecZ=zero
#ifdef OverWriteFFT
      allocate(Arr(z1size(1),z1size(2),z1size(3)))
      fft_plan_z = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,Arr,Arr,[FFTW_R2HC],plan_type)
      deallocate(Arr)
#else
      allocate(Arr1(z1size(1),z1size(2),z1size(3)),Arr2(z1size(1),z1size(2),z1size(3)))
      fft_plan_z = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,Arr1,Arr2,[FFTW_R2HC],plan_type)
      deallocate(Arr1,Arr2)
#endif
    ENDIF
    if(nySpec2D<1) clcSpectra2D=.false.
    IF(clcSpectra2D) THEN
      allocate(iySpec2D(nySpec2D))
      nUnit = GetFileUnit() 
      open(unit=nUnit, file=ChannelPrm, status='old',form='formatted',IOSTAT=ierror)
      read(nUnit, nml=Spectra2DOpt)
      close(nUnit,IOSTAT=ierror)
      if(nySpec2D*NEnergySpec2D>nyc .and. nrank==0) then
        call MainLog%CheckForError(ErrT_Pass,"InitStatVar","So big nySpec2D")        
      endif 
      allocate(EnergySpec2D(y2size(1),nySpec2D,y2size(3),NEnergySpec2D))
      if(FlowType==FT_CH) then
        do jc=1,nySpec2D
          if((iySpec2D(jc)<2 .or. iySpec2D(jc)>nyc/2) .and. nrank==0) then
            call MainLog%CheckForError(ErrT_Abort,"InitStatVar","iySpec2D WRONG 1")
          endif
        enddo
      else
        do jc=1,nySpec2D
          if((iySpec2D(jc)<2 .or. iySpec2D(jc)>nyc)   .and. nrank==0) then
            call MainLog%CheckForError(ErrT_Abort,"InitStatVar","iySpec2D WRONG 2")
          endif
        enddo   
      endif
      EnergySpec2D=zero

      iodim(1)%n  = z2size(3)
      iodim(1)%is = z2size(1)*z2size(2)
      iodim(1)%os = z2size(1)*z2size(2)
      iodim_howmany(1)%n  = z2size(1)*z2size(2)
      iodim_howmany(1)%is = 1
      iodim_howmany(1)%os = 1
#ifdef OverWriteFFT
      allocate(Arr(z2size(1),z2size(2),z2size(3)))
      fft_plan_z2 = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,Arr,Arr,[FFTW_R2HC],plan_type)
      deallocate(Arr)
#else
      allocate(Arr1(z2size(1),z2size(2),z2size(3)),Arr2(z2size(1),z2size(2),z2size(3)))
      fft_plan_z2 = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,Arr1,Arr2,[FFTW_R2HC],plan_type)
      deallocate(Arr1,Arr2)
#endif
    ENDIF
  end subroutine InitStatVar

  !******************************************************************
  ! clcStat
  !******************************************************************
#ifdef ScalarFlow
  subroutine clcStat(ux,uy,uz,pressure,scalar,ArrTemp1,ArrTemp2)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure,scalar
#else
  subroutine clcStat(ux,uy,uz,pressure,ArrTemp1,ArrTemp2)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
#endif
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(inout)::ArrTemp1,ArrTemp2
   
    ! locals
    character(len=50)::filename
    integer(kind=MPI_OFFSET_KIND)::disp,disp_inc
    real(RK),dimension(:,:,:),allocatable::arrYplane
    real(RK),allocatable,dimension(:,:,:)::EnergySpecXR,EnergySpecZR
    real(RK),dimension(:,:,:),allocatable::arrx1,arrx2,arrz1,arrz2
#ifndef OverWriteFFT
    real(RK),dimension(:,:,:),allocatable::ArrFFT
#endif
    integer::ic,jc,kc,im,jm,km,ip,jp,kp,myistat,ierror,is,ks,iu,ku,it,jt,kt,nrankX,nrankZ,nUnit
    real(RK)::uxloc,uyloc,uzloc,prloc,uxCell,uyCell,uzCell,prloc2,inxz,infstime,cac,caj,cacU,dudyU,uy_xc,uy_xp
    real(RK)::dudxx,dudyy,dudzz,dvdxx,dvdyy,dvdzz,dwdxx,dwdyy,dwdzz,SumStatR(NCHASTAT,nyp),SumVec(NCHASTAT),rdxt,rdzt
    real(RK)::dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,vor_x,vor_y,vor_z,InterpY1,InterpY2,dudzC,dvdzC,dudxC,dvdxU
#ifdef ScalarFlow
    real(RK)::scloc
#endif

    rdxt=rdx/12.0_RK
    rdzt=rdz/12.0_RK
    inxz=one/(real(nxc,RK)*real(nzc,RK))
    DO jc=1,nyc
      jp=jc+1; jm=jc-1;
      InterpY1= YinterpCoe(jm); InterpY2=one-InterpY1
      cac=rdyc(jc);cacU=rdyc(jp); caj=rdyp(jc); SumVec=zero
      do kc=y1start(3),y1end(3)
        ks=kc-2;km=kc-1;kp=kc+1;ku=kc+2
        do ic=y1start(1),y1end(1)
          is=ic-2;im=ic-1;ip=ic+1;iu=ic+2

          uxloc= ux(ic,jc,kc)+uCRF
          uyloc= uy(ic,jc,kc)
          uzloc= uz(ic,jc,kc)
          prloc= pressure(ic,jc,kc)
          uxCell= InterpCoe1*ux(im,jc,kc) +InterpCoe2*ux(ic,jc,kc) +InterpCoe3*ux(ip,jc,kc) +InterpCoe4*ux(iu,jc,kc)+uCRF
          uyCell= (uyloc+ uy(ic,jp,kc))*half
          uy_xc = InterpCoe1*uy(is,jc,kc) +InterpCoe2*uy(im,jc,kc) +InterpCoe3*uy(ic,jc,kc) +InterpCoe4*uy(ip,jc,kc)
          uy_xp = InterpCoe1*uy(is,jp,kc) +InterpCoe2*uy(im,jp,kc) +InterpCoe3*uy(ic,jp,kc) +InterpCoe4*uy(ip,jp,kc)
          uzCell= InterpCoe1*uz(ic,jc,km) +InterpCoe2*uz(ic,jc,kc) +InterpCoe3*uz(ic,jc,kp) +InterpCoe4*uz(ic,jc,ku) 
          prloc2= InterpY1*(InterpCoe1*pressure(is,jm,kc) +InterpCoe2*pressure(im,jm,kc) + &
                            InterpCoe3*pressure(ic,jm,kc) +InterpCoe4*pressure(ip,jm,kc))+ &
                  InterpY2*(InterpCoe1*pressure(is,jc,kc) +InterpCoe2*pressure(im,jc,kc) + &
                            InterpCoe3*pressure(ic,jc,kc) +InterpCoe4*pressure(ip,jc,kc))  !2021-12-14
 
          dudx= dxCoe1*ux(im,jc,kc) +dxCoe2*ux(ic,jc,kc) +dxCoe3*ux(ip,jc,kc) +dxCoe4*ux(iu,jc,kc)
          dudy= (ux(ic,jc,kc)-ux(ic,jm,kc))*cac
          dudz= dzCoe1*ux(ic,jc,ks) +dzCoe2*ux(ic,jc,km) +dzCoe3*ux(ic,jc,kc) +dzCoe4*ux(ic,jc,kp)
          dvdx= dxCoe1*uy(is,jc,kc) +dxCoe2*uy(im,jc,kc) +dxCoe3*uy(ic,jc,kc) +dxCoe4*uy(ip,jc,kc)
          dvdy= (uy(ic,jp,kc)-uyloc)*caj
          dvdz= dzCoe1*uy(ic,jc,ks) +dzCoe2*uy(ic,jc,km) +dzCoe3*uy(ic,jc,kc) +dzCoe4*uy(ic,jc,kp)
          dwdx= dxCoe1*uz(is,jc,kc) +dxCoe2*uz(im,jc,kc) +dxCoe3*uz(ic,jc,kc) +dxCoe4*uz(ip,jc,kc)
          dwdy= (uzloc-uz(ic,jm,kc))*cac
          dwdz= dzCoe1*uz(ic,jc,km) +dzCoe2*uz(ic,jc,kc) +dzCoe3*uz(ic,jc,kp) +dzCoe4*uz(ic,jc,ku)

          dudxx= dxxCoe1*ux(is,jc,kc) +dxxCoe2*ux(im,jc,kc) +dxxCoe3*ux(ic,jc,kc) +dxxCoe4*ux(ip,jc,kc) +dxxCoe5*ux(iu,jc,kc)
          dudyy= ap2c(jc)*ux(ic,jp,kc)+ac2c(jc)*ux(ic,jc,kc)+am2c(jc)*ux(ic,jm,kc)
          dudzz= dzzCoe1*ux(ic,jc,ks) +dzzCoe2*ux(ic,jc,km) +dzzCoe3*ux(ic,jc,kc) +dzzCoe4*ux(ic,jc,kp) +dzzCoe5*ux(ic,jc,ku)
          dvdxx= dxxCoe1*uy(is,jc,kc) +dxxCoe2*uy(im,jc,kc) +dxxCoe3*uy(ic,jc,kc) +dxxCoe4*uy(ip,jc,kc) +dxxCoe5*uy(iu,jc,kc)
          dvdyy= ap2p(jc)*uy(ic,jp,kc)+ac2p(jc)*uy(ic,jc,kc)+am2p(jc)*uy(ic,jm,kc)
          dvdzz= dzzCoe1*uy(ic,jc,ks) +dzzCoe2*uy(ic,jc,km) +dzzCoe3*uy(ic,jc,kc) +dzzCoe4*uy(ic,jc,kp) +dzzCoe5*uy(ic,jc,ku)
          dwdxx= dxxCoe1*uz(is,jc,kc) +dxxCoe2*uz(im,jc,kc) +dxxCoe3*uz(ic,jc,kc) +dxxCoe4*uz(ip,jc,kc) +dxxCoe5*uz(iu,jc,kc)
          dwdyy= ap2c(jc)*uz(ic,jp,kc)+ac2c(jc)*uz(ic,jc,kc)+am2c(jc)*uz(ic,jm,kc)
          dwdzz= dzzCoe1*uz(ic,jc,ks) +dzzCoe2*uz(ic,jc,km) +dzzCoe3*uz(ic,jc,kc) +dzzCoe4*uz(ic,jc,kp) +dzzCoe5*uz(ic,jc,ku)
          vor_x= dwdy-dvdz
          vor_y= dudz-dwdx
          vor_z= dvdx-dudy

          dudxC= (-ux(iu,jc,kc)+8.0_RK*ux(ip,jc,kc)-8.0_RK*ux(im,jc,kc)+ux(is,jc,kc))*rdxt
          dvdxU= dxCoe1*uy(is,jp,kc) +dxCoe2*uy(im,jp,kc) +dxCoe3*uy(ic,jp,kc) +dxCoe4*uy(ip,jp,kc)
          dudyU= (ux(ic,jp,kc)-ux(ic,jc,kc))*cacU
          dudzC= (InterpY1*(-ux(ic,jm,ku)+8.0_RK*ux(ic,jm,kp)-8.0_RK*ux(ic,jm,km)+ux(ic,jm,ks) ) +&
                  InterpY2*(-ux(ic,jc,ku)+8.0_RK*ux(ic,jc,kp)-8.0_RK*ux(ic,jc,km)+ux(ic,jc,ks) ))*rdzt !2021-12-14
          dvdzC=((-uy(is,jc,ku)+8.0_RK*uy(is,jc,kp)-8.0_RK*uy(is,jc,km)+uy(is,jc,ks))*InterpCoe1+ &
                 (-uy(im,jc,ku)+8.0_RK*uy(im,jc,kp)-8.0_RK*uy(im,jc,km)+uy(im,jc,ks))*InterpCoe2+ &
                 (-uy(ic,jc,ku)+8.0_RK*uy(ic,jc,kp)-8.0_RK*uy(ic,jc,km)+uy(ic,jc,ks))*InterpCoe3+ &
                 (-uy(ip,jc,ku)+8.0_RK*uy(ip,jc,kp)-8.0_RK*uy(ip,jc,km)+uy(ip,jc,ks))*InterpCoe4)*rdzt !2021-12-14

          SumVec( 1)=SumVec( 1)+ uxloc
          SumVec( 2)=SumVec( 2)+ uyloc                     ! yp
          SumVec( 3)=SumVec( 3)+ uzloc
          SumVec( 4)=SumVec( 4)+ prloc
          SumVec( 5)=SumVec( 5)+ uxloc*uxloc
          SumVec( 6)=SumVec( 6)+ uyloc*uyloc               ! yp  
          SumVec( 7)=SumVec( 7)+ uzloc*uzloc
          SumVec( 8)=SumVec( 8)+ prloc*prloc
          SumVec( 9)=SumVec( 9)+ uxCell*uyCell
          SumVec(10)=SumVec(10)+ uyCell*uzCell
          SumVec(11)=SumVec(11)+ uxCell*uzCell
          SumVec(12)=SumVec(12)+ uxCell*prloc
          SumVec(13)=SumVec(13)+ uyCell*prloc
          SumVec(14)=SumVec(14)+ uzCell*prloc
          SumVec(15)=SumVec(15)+ uxCell*uxCell*uyCell 
          SumVec(16)=SumVec(16)+ uyloc *uyloc *uyloc       ! yp
          SumVec(17)=SumVec(17)+ uzCell*uzCell*uyCell
          SumVec(18)=SumVec(18)+ uxCell*uyCell*uyCell
          SumVec(19)=SumVec(19)+ uxloc*uxloc*uxloc
          SumVec(20)=SumVec(20)+ uzloc*uzloc*uzloc
          SumVec(21)=SumVec(21)+ uxloc*uxloc*uxloc*uxloc
          SumVec(22)=SumVec(22)+ uyloc*uyloc*uyloc*uyloc   ! yp
          SumVec(23)=SumVec(23)+ uzloc*uzloc*uzloc*uzloc
          SumVec(24)=SumVec(24)+ prloc*dudx
          SumVec(25)=SumVec(25)+ prloc*dvdy
          SumVec(26)=SumVec(26)+ prloc*dwdz
          SumVec(27)=SumVec(27)+ prloc2*(dudy+dvdx)        ! yp !
          SumVec(28)=SumVec(28)+ uxloc*(dudxx+dudyy+dudzz)
          SumVec(29)=SumVec(29)+ uyloc*(dvdxx+dvdyy+dvdzz) ! yp !
          SumVec(30)=SumVec(30)+ uzloc*(dwdxx+dwdyy+dwdzz)
          SumVec(31)=SumVec(31)+ (dudxC*(dvdxU+dvdx)+ (dudyU+dudy)*(uy_xp-uy_xc)*caj)*half !2021-12-14
          SumVec(32)=SumVec(32)+ dudzC*dvdzC          ! yp !
          SumVec(33)=SumVec(33)+ vor_x*vor_x          ! yp !
          SumVec(34)=SumVec(34)+ vor_y*vor_y
          SumVec(35)=SumVec(35)+ vor_z*vor_z          ! yp !
#ifdef ScalarFlow
          scloc=scalar(ic,jc,kc)
          SumVec(36)=SumVec(36)+ scloc
          SumVec(37)=SumVec(37)+ scloc*scloc
          SumVec(38)=SumVec(38)+ scloc*uxCell
          SumVec(39)=SumVec(39)+ scloc*uyCell
          SumVec(40)=SumVec(40)+ scloc*uzCell
#endif
        enddo
      enddo
      do kc=1,NCHASTAT
        SumStat(kc,jc)=SumStat(kc,jc)+ SumVec(kc)*inxz
      enddo
    ENDDO

    ! nyp only
    jc=nyp; jm=jc-1; cac=rdyc(jc); SumVec=zero
    InterpY1= YinterpCoe(jm); InterpY2=one-InterpY1
    do kc=y1start(3),y1end(3)
      km=kc-1;kp=kc+1
      do ic=y1start(1),y1end(1)
        is=ic-2;im=ic-1;ip=ic+1
        prloc2= InterpY1*(InterpCoe1*pressure(is,jm,kc) +InterpCoe2*pressure(im,jm,kc) +InterpCoe3*pressure(ic,jm,kc) +InterpCoe4*pressure(ip,jm,kc))+ &
                InterpY2*(InterpCoe1*pressure(is,jc,kc) +InterpCoe2*pressure(im,jc,kc) +InterpCoe3*pressure(ic,jc,kc) +InterpCoe4*pressure(ip,jc,kc)) !2021-12-14
        dudy= (ux(ic,jc,kc)-ux(ic,jm,kc))*cac
        dwdy= (uz(ic,jc,kc)-uz(ic,jm,kc))*cac
        vor_x=  dwdy
        vor_z= -dudy
        SumVec(27)=SumVec(27)+ prloc2*(dudy+zero)   ! yp !
        SumVec(33)=SumVec(33)+ vor_x*vor_x          ! yp !
        SumVec(35)=SumVec(35)+ vor_z*vor_z          ! yp !
      enddo
    enddo
    do kc=1,NCHASTAT
      SumStat(kc,jc)=SumStat(kc,jc)+ SumVec(kc)*inxz
    enddo

    ! shear stress and pressure gradient
    PrGradsum   = PrGradsum+ PrGradAver
    if(nrank==0 .and. IsUxConst) then
      nUnit=GetFileUnit()
      write(filename,'(A,I10.10)')trim(Res_dir)//"PrGrad",ilast
      open(nUnit,file=filename,status='old',position='append',form='formatted',IOSTAT=myistat)
      if(myistat/=0) then
        call MainLog%CheckForError(ErrT_Pass,"clcStat_CH","Cannot open file: "//trim(filename))
      else
        write(nUnit,'(I7,2ES24.15)')itime,SimTime,PrGradAver
      endif
      close(nUnit,IOSTAT=myistat)
    endif
    nfstime= nfstime + 1

    ! Calculate 1-D Energy Spectra
    IF(mod(itime,ivSpec)==0) THEN
    
      ! Determine ux in cell center
      DO kc=y1start(3),y1end(3)
        kt=kc-y1start(3)+1
        do jc=y1start(2),y1end(2)
          jt=jc
          do ic=y1start(1),y1end(1)
            it=ic-y1start(1)+1
            im=ic-1;ip=ic+1;iu=ic+2
            ArrTemp1(it,jt,kt)=InterpCoe1*ux(im,jc,kc) +InterpCoe2*ux(ic,jc,kc) + &
                               InterpCoe3*ux(ip,jc,kc) +InterpCoe4*ux(iu,jc,kc) +uCRF
          enddo
        enddo
      ENDDO
        
      IF(clcSpectra1D) THEN

        !=============== Spectra in x-dir ===============
        if(FlowType==FT_HC) then
          allocate(arrYplane(y1size(1),y1size(3),1))
        else
          allocate(arrYplane(y1size(1),y1size(3),2))        
        endif
        allocate(arrx1(x1size(1),x1size(2),x1size(3)))
        allocate(arrx2(x1size(1),x1size(2),x1size(3)))
#ifndef OverWriteFFT
        allocate(arrFFT(x1size(1),x1size(2),x1size(3)))
#endif

#ifdef OverWriteFFT
        call transpose_y1_to_x1(ux(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrx1); arrx1=arrx1+uCRF
        call dfftw_execute_r2r(fft_plan_x,arrx1,arrx1)
#else
        call transpose_y1_to_x1(ux(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrFFT); arrFFT=arrFFT+uCRF
        call dfftw_execute_r2r(fft_plan_x,arrFFT,arrx1)
#endif
        call clcEnergySpectraX(arrx1,1)
      
        ! Calculate Linear coherence spectrum (LCS)
        ! Real part
        call transpose_x1_to_y1(arrx1,ArrTemp2)
        IF(FlowType==FT_HC) THEN
          do kc=1,y1size(3)
            do ic=1,y1size(1)
              arrYplane(ic,kc,1)=ArrTemp2(ic,jForLCS,kc)
            enddo
          enddo
        ELSE
          jt=nyc+1-jForLCS
          do kc=1,y1size(3)
            do ic=1,y1size(1)
              arrYplane(ic,kc,1)=ArrTemp2(ic,jForLCS,kc)
              arrYplane(ic,kc,2)=ArrTemp2(ic,jt,kc)
            enddo
          enddo        
        ENDIF
        
        IF(FlowType==FT_HC) THEN
          do kc=1,y1size(3)
            do jc=1,y1size(2)
              do ic=1,y1size(1)
                ArrTemp2(ic,jc,kc)=ArrTemp2(ic,jc,kc)*arrYplane(ic,kc,1)
              enddo
            enddo
          enddo            
        ELSE
          do kc=1,y1size(3)
            do jc=1,nyc/2
              do ic=1,y1size(1)
                ArrTemp2(ic,jc,kc)=ArrTemp2(ic,jc,kc)*arrYplane(ic,kc,1)
              enddo
            enddo
            do jc=nyc/2+1,nyc
              do ic=1,y1size(1)
                ArrTemp2(ic,jc,kc)=ArrTemp2(ic,jc,kc)*arrYplane(ic,kc,2)
              enddo
            enddo
          enddo            
        ENDIF
        call transpose_y1_to_x1(ArrTemp2,arrx2)
        call clcLTS_x_real(arrx2,7)
      
        ! Imaginary part
        arrx2(1,:,:)=arrx1(1,:,:)
        do kc=1,x1size(3)
          do jc=1,x1size(2)
            do ic=2,x1size(1)
              arrx2(ic,jc,kc)=arrx1(nxc+2-ic,jc,kc)
            enddo
          enddo
        enddo
        call transpose_x1_to_y1(arrx2,ArrTemp2)
        IF(FlowType==FT_HC) THEN
          do kc=1,y1size(3)
            do jc=1,y1size(2)
              do ic=1,y1size(1)
                ArrTemp2(ic,jc,kc)=ArrTemp2(ic,jc,kc)*arrYplane(ic,kc,1)
              enddo
            enddo
          enddo
        ELSE
          do kc=1,y1size(3)
            do jc=1,nyc/2
              do ic=1,y1size(1)
                ArrTemp2(ic,jc,kc)=ArrTemp2(ic,jc,kc)*arrYplane(ic,kc,1)
              enddo
            enddo
            do jc=nyc/2+1,nyc
              do ic=1,y1size(1)
                ArrTemp2(ic,jc,kc)=ArrTemp2(ic,jc,kc)*arrYplane(ic,kc,2)
              enddo
            enddo
          enddo        
        ENDIF
        call transpose_y1_to_x1(ArrTemp2,arrx2)
        call clcLTS_x_imag(arrx2,8)

#ifdef OverWriteFFT
        call transpose_y1_to_x1(uy(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrx2)
        call dfftw_execute_r2r(fft_plan_x,arrx2,arrx2)
        call clcEnergySpectraX(arrx2,2)
        
        call transpose_y1_to_x1(uz(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrx2)
        call dfftw_execute_r2r(fft_plan_x,arrx2,arrx2)
        call clcEnergySpectraX(arrx2,3)
        
        call transpose_y1_to_x1(pressure(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrx2)
        call dfftw_execute_r2r(fft_plan_x,arrx2,arrx2)
        call clcEnergySpectraX(arrx2,4)
#else
        call transpose_y1_to_x1(uy(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrFFT)
        call dfftw_execute_r2r(fft_plan_x,arrFFT,arrx2)
        call clcEnergySpectraX(arrx2,2)

        call transpose_y1_to_x1(uz(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrFFT)
        call dfftw_execute_r2r(fft_plan_x,arrFFT,arrx2)
        call clcEnergySpectraX(arrx2,3)
        
        call transpose_y1_to_x1(pressure(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrFFT)
        call dfftw_execute_r2r(fft_plan_x,arrFFT,arrx2)
        call clcEnergySpectraX(arrx2,4)
#endif

        ! Determine uy in cell center
        DO kc=y1start(3),y1end(3)
          kt=kc-y1start(3)+1
          do jc=y1start(2),y1end(2)
            jt=jc
            jp=jc+1
            do ic=y1start(1),y1end(1)
              it=ic-y1start(1)+1
              ArrTemp2(it,jt,kt)=half*(uy(ic,jc,kc)+uy(ic,jp,kc))
            enddo
          enddo
        ENDDO

#ifdef OverWriteFFT
        call transpose_y1_to_x1(ArrTemp2,arrx2)
        call dfftw_execute_r2r(fft_plan_x,arrx2,arrx2)
        call clcCospectraX(arrx1,arrx2,5)
        call transpose_y1_to_x1(ArrTemp1,arrx1)
        call dfftw_execute_r2r(fft_plan_x,arrx1,arrx1)
        call clcCospectraX(arrx1,arrx2,6)
#ifdef ScalarFlow
        call transpose_y1_to_x1(scalar(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrx2)
        call dfftw_execute_r2r(fft_plan_x,arrx2,arrx2)
        call clcEnergySpectraX(arrx2,9)
#endif
        deallocate(arrx1,arrx2)
#else
        call transpose_y1_to_x1(ArrTemp2,arrFFT)
        call dfftw_execute_r2r(fft_plan_x,arrFFT,arrx2)
        call clcCospectraX(arrx1,arrx2,5)
        call transpose_y1_to_x1(ArrTemp1,arrFFT)
        call dfftw_execute_r2r(fft_plan_x,arrFFT,arrx1)
        call clcCospectraX(arrx1,arrx2,6)
#ifdef ScalarFlow
        call transpose_y1_to_x1(scalar(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrFFT)
        call dfftw_execute_r2r(fft_plan_x,arrFFT,arrx2)
        call clcEnergySpectraX(arrx2,9)
#endif
        deallocate(arrx1,arrx2,arrFFT)
#endif
       
        !=============== Spectra in z-dir ===============
        allocate(arrz1(z1size(1),z1size(2),z1size(3)))
        allocate(arrz2(z1size(1),z1size(2),z1size(3)))
#ifndef OverWriteFFT
        allocate(arrFFT(z1size(1),z1size(2),z1size(3)))
#endif

#ifdef OverWriteFFT
        call transpose_y1_to_z1(ux(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrz1); arrz1=arrz1+uCRF
        call dfftw_execute_r2r(fft_plan_z,arrz1,arrz1)
#else
        call transpose_y1_to_z1(ux(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrFFT); arrFFT=arrFFT+uCRF
        call dfftw_execute_r2r(fft_plan_z,arrFFT,arrz1)
#endif
        call clcEnergySpectraZ(arrz1,1)

        ! Calculate Linear coherence spectrum (LCS)
        ! Real part
        call transpose_z1_to_y1(arrz1,ArrTemp2)
        IF(FlowType==FT_HC) THEN
          do kc=1,y1size(3)
            do ic=1,y1size(1)
              arrYplane(ic,kc,1)=ArrTemp2(ic,jForLCS,kc)
            enddo
          enddo
        ELSE
          jt=nyc+1-jForLCS
          do kc=1,y1size(3)
            do ic=1,y1size(1)
              arrYplane(ic,kc,1)=ArrTemp2(ic,jForLCS,kc)
              arrYplane(ic,kc,2)=ArrTemp2(ic,jt,kc)
            enddo
          enddo        
        ENDIF

        IF(FlowType==FT_HC) THEN        
          do kc=1,y1size(3)
            do jc=1,y1size(2)
              do ic=1,y1size(1)
                ArrTemp2(ic,jc,kc)=ArrTemp2(ic,jc,kc)*arrYplane(ic,kc,1)
              enddo
            enddo
          enddo
        ELSE
          do kc=1,y1size(3)
            do jc=1,nyc/2
              do ic=1,y1size(1)
                ArrTemp2(ic,jc,kc)=ArrTemp2(ic,jc,kc)*arrYplane(ic,kc,1)
              enddo
            enddo
            do jc=nyc/2+1,nyc
              do ic=1,y1size(1)
                ArrTemp2(ic,jc,kc)=ArrTemp2(ic,jc,kc)*arrYplane(ic,kc,2)
              enddo
            enddo            
          enddo       
        ENDIF
        call transpose_y1_to_z1(ArrTemp2,arrz2)
        call clcLTS_z_real(arrz2,7)
      
        ! Imaginary part
        arrz2(:,:,1)=arrz1(:,:,1)
        do kc=2,z1size(3)
          do jc=1,z1size(2)
            do ic=1,z1size(1)
              arrz2(ic,jc,kc)=arrz1(ic,jc,nzc+2-kc)
            enddo
          enddo
        enddo
        call transpose_z1_to_y1(arrz2,ArrTemp2)
        IF(FlowType==FT_HC) THEN
          do kc=1,y1size(3)
            do jc=1,y1size(2)
              do ic=1,y1size(1)
                ArrTemp2(ic,jc,kc)=ArrTemp2(ic,jc,kc)*arrYplane(ic,kc,1)
              enddo
            enddo
          enddo
        ELSE
          do kc=1,y1size(3)
            do jc=1,nyc/2
              do ic=1,y1size(1)
                ArrTemp2(ic,jc,kc)=ArrTemp2(ic,jc,kc)*arrYplane(ic,kc,1)
              enddo
            enddo
            do jc=nyc/2+1,nyc
              do ic=1,y1size(1)
                ArrTemp2(ic,jc,kc)=ArrTemp2(ic,jc,kc)*arrYplane(ic,kc,2)
              enddo
            enddo
          enddo      
        ENDIF
        call transpose_y1_to_z1(ArrTemp2,arrz2)
        call clcLTS_z_imag(arrz2,8)

#ifdef OverWriteFFT
        call transpose_y1_to_z1(uy(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrz2)
        call dfftw_execute_r2r(fft_plan_z,arrz2,arrz2)
        call clcEnergySpectraZ(arrz2,2)

        call transpose_y1_to_z1(uz(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrz2)
        call dfftw_execute_r2r(fft_plan_z,arrz2,arrz2)
        call clcEnergySpectraZ(arrz2,3)

        call transpose_y1_to_z1(pressure(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrz2)
        call dfftw_execute_r2r(fft_plan_z,arrz2,arrz2)
        call clcEnergySpectraZ(arrz2,4)
#else
        call transpose_y1_to_z1(uy(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrFFT)
        call dfftw_execute_r2r(fft_plan_z,arrFFT,arrz2)
        call clcEnergySpectraZ(arrz2,2)

        call transpose_y1_to_z1(uz(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrFFT)
        call dfftw_execute_r2r(fft_plan_z,arrFFT,arrz2)
        call clcEnergySpectraZ(arrz2,3)

        call transpose_y1_to_z1(pressure(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrFFT)
        call dfftw_execute_r2r(fft_plan_z,arrFFT,arrz2)
        call clcEnergySpectraZ(arrz2,4)
#endif

        ! Determine uy in cell center
        DO kc=y1start(3),y1end(3)
          kt=kc-y1start(3)+1
          do jc=y1start(2),y1end(2)
            jt=jc
            jp=jc+1
            do ic=y1start(1),y1end(1)
              it=ic-y1start(1)+1
              ArrTemp2(it,jt,kt)=half*(uy(ic,jc,kc)+uy(ic,jp,kc))
            enddo
          enddo
        ENDDO     
#ifdef OverWriteFFT
        call transpose_y1_to_z1(ArrTemp2,arrz2)
        call dfftw_execute_r2r(fft_plan_z,arrz2,arrz2)
        call clcCospectraZ(arrz1,arrz2,5)
        call transpose_y1_to_z1(ArrTemp1,arrz1)
        call dfftw_execute_r2r(fft_plan_z,arrz1,arrz1)
        call clcCospectraZ(arrz1,arrz2,6)     
#ifdef ScalarFlow
        call transpose_y1_to_z1(scalar(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrz2)
        call dfftw_execute_r2r(fft_plan_z,arrz2,arrz2)
        call clcEnergySpectraZ(arrz2,9)
#endif
        deallocate(arrz1,arrz2)
#else
        call transpose_y1_to_z1(ArrTemp2,arrFFT)
        call dfftw_execute_r2r(fft_plan_z,arrFFT,arrz2)
        call clcCospectraZ(arrz1,arrz2,5)
        call transpose_y1_to_z1(ArrTemp1,arrFFT)
        call dfftw_execute_r2r(fft_plan_z,arrFFT,arrz1)
        call clcCospectraZ(arrz1,arrz2,6)     
#ifdef ScalarFlow
        call transpose_y1_to_z1(scalar(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrFFT)
        call dfftw_execute_r2r(fft_plan_z,arrFFT,arrz2)
        call clcEnergySpectraZ(arrz2,9)
#endif
        deallocate(arrz1,arrz2,arrFFT)
#endif

        deallocate(arrYplane)
      ENDIF      

      ! Calculate 2-D Energy Spectra
      IF(clcSpectra2D) THEN
      
        ! Determine uy in cell center
        DO kc=y1start(3),y1end(3)
          kt=kc-y1start(3)+1
          do jc=y1start(2),y1end(2)
            jt=jc
            jp=jc+1
            do ic=y1start(1),y1end(1)
              it=ic-y1start(1)+1
              ArrTemp2(it,jt,kt)=half*(uy(ic,jc,kc)+uy(ic,jp,kc))
            enddo
          enddo
        ENDDO
        call clcEnergySpectra2D(ux(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),1,.false.)
        call clcEnergySpectra2D(uy(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),2,.true.)
        call clcEnergySpectra2D(uz(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),3,.false.)
        call clcEnergySpectra2D(pressure(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),4,.false.)
        call clcCosSpectra2D(ArrTemp1,ArrTemp2,5,.false.)
#ifdef ScalarFlow
        call clcEnergySpectra2D(scalar(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),6,.false.)
#endif
      ENDIF
    
      nSpectime= nSpectime+1
    ENDIF
    if(mod(itime,SaveStat)/=0) return
    
    ! Write statistics
    call MPI_REDUCE(SumStat,SumStatR,NCHASTAT*nyp,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nrank==0) then
      infstime = one/real(nfstime,RK)
      nUnit=GetFileUnit()
      write(filename,"(A,I10.10)") trim(Res_Dir)//'stats',itime
      open(nUnit,file=filename,status='replace',form='formatted',IOSTAT=myistat)
      IF(myistat/=0) THEN
        call MainLog%CheckForError(ErrT_Abort,"clcStat_CH","Cannot open file: "//trim(filename))
      ELSE
        write(nUnit,'(a,I7,a,I7,a,I7)')'  The time step range for this fluid statistics is ', &
                                    itime-(nfstime-1)*ivstats, ':', ivstats, ':', itime
        write(nUnit,'(A)')'  '
        if(IsUxConst) then
          write(nUnit,'(A)')'  Constant velocity in x-dir by adding a pressure gradient.'
          write(nUnit,'(A, ES24.15)')'    time averaged pressure gradient is: ',abs(PrGradsum)*infstime
        else
          write(nUnit,'(A)')'  Variable velocity in x-dir while adding a constant body force.'
        endif
        write(nUnit,'(A)')'  '
        
        Block 
        character(len=50)::FormatStr
        write(FormatStr,'(A,I3,A)')'(',NCHASTAT,'ES24.15)'
        do jc=1,nyp
          write(nUnit,FormatStr)SumStatR(1:NCHASTAT,jc)*infstime
        enddo
        End Block
      ENDIF
      close(nUnit,IOSTAT=myistat)
    endif

    ! Write 1-D Energy Spectra
    IF(clcSpectra1D) THEN
      infstime= one/real(nSpectime,RK)

      ! Spectra in x-dir
      allocate(EnergySpecXR(nxhp,x1size(2),NEnergySpec1D));EnergySpecXR=zero
      call MPI_REDUCE(EnergySpecX,EnergySpecXR,NEnergySpec1D*x1size(2)*nxhp,real_type,MPI_SUM,0,DECOMP_2D_COMM_ROW,ierror)
      call MPI_COMM_RANK(DECOMP_2D_COMM_ROW,nrankX,ierror)
      if(nrankX==0) then
        EnergySpecXR=EnergySpecXR*infstime
        disp_inc=int(mytype_bytes,8)*int(nyc,8)*int(nxhp,8)
        disp=int(mytype_bytes,8)*int(x1start(2)-1,8)*int(nxhp,8)
        
        write(filename,"(A,A,I10.10)") trim(Res_Dir),'SpecX',itime
        call MPI_FILE_OPEN(DECOMP_2D_COMM_COL,filename,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,nUnit,ierror)
        call MPI_FILE_SET_SIZE(nUnit,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
        do kc=1,NEnergySpec1D
          call MPI_FILE_WRITE_AT_ALL(nUnit,disp,EnergySpecXR(:,:,kc),x1size(2)*nxhp,real_type,MPI_STATUS_IGNORE,ierror)
          disp=disp+disp_inc
        enddo
        call MPI_FILE_CLOSE(nUnit,ierror)
      endif
      deallocate(EnergySpecXR)
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)

      ! Spectra in z-dir
      allocate(EnergySpecZR(nzhp,z1size(2),NEnergySpec1D));EnergySpecZR=zero
      call MPI_REDUCE(EnergySpecZ,EnergySpecZR,NEnergySpec1D*z1size(2)*nzhp,real_type,MPI_SUM,0,DECOMP_2D_COMM_COL,ierror)
      call MPI_COMM_RANK(DECOMP_2D_COMM_COL,nrankZ,ierror)
      if(nrankZ==0) then
        EnergySpecZR=EnergySpecZR*infstime
        disp_inc=int(mytype_bytes,8)*int(nyc,8)*int(nzhp,8)
        disp=int(mytype_bytes,8)*int(z1start(2)-1,8)*int(nzhp,8)
        
        write(filename,"(A,A,I10.10)") trim(Res_Dir),'SpecZ',itime
        call MPI_FILE_OPEN(DECOMP_2D_COMM_ROW,filename,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,nUnit,ierror)
        call MPI_FILE_SET_SIZE(nUnit,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
        do kc=1,NEnergySpec1D
          call MPI_FILE_WRITE_AT_ALL(nUnit,disp,EnergySpecZR(:,:,kc),z1size(2)*nzhp,real_type,MPI_STATUS_IGNORE,ierror)
          disp=disp+disp_inc
        enddo
        call MPI_FILE_CLOSE(nUnit,ierror)
      endif
      deallocate(EnergySpecZR)
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    ENDIF

    ! Write 2-D Energy Spectra
    IF(clcSpectra2D) THEN
    
      infstime= one/real(nSpectime,RK)
      EnergySpec2D=EnergySpec2D*infstime
      Block
      integer(MPI_OFFSET_KIND)::disp
      integer::data_type,data_byte,newtype
      integer,dimension(3)::sizes,subsizes,starts
#ifdef SAVE_SINGLE_Spec2D
      real(4),allocatable,dimension(:,:,:,:)::EnergySpecOut
      data_type= MPI_REAL
      allocate(EnergySpecOut(y2size(1),nySpec2D,y2size(3),NEnergySpec2D))
      do kt=1,NEnergySpec2D
        do kc=1,y2size(3)
          do jc=1,nySpec2D
            do ic=1,y2size(1)
              EnergySpecOut(ic,jc,kc,kt)=real(EnergySpec2D(ic,jc,kc,kt),4)
            enddo
          enddo
        enddo
      enddo
#else
      data_type= MPI_DOUBLE_PRECISION
#endif
      call MPI_TYPE_SIZE(data_type,data_byte,ierror)
      sizes(1)= nxc
      sizes(2)= nySpec2D
      sizes(3)= nzc  
      subsizes(1)=y2size(1)
      subsizes(2)=nySpec2D
      subsizes(3)=y2size(3)
      starts(1)=y2start(1)-1
      starts(2)=0
      starts(3)=y2start(3)-1
      call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,data_type,newtype,ierror)
      call MPI_TYPE_COMMIT(newtype,ierror)
      
      disp = 0_MPI_OFFSET_KIND
      write(filename,"(A,A,I10.10)") trim(Res_Dir),'Spec2D',itime
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, nUnit, ierror)
      call MPI_FILE_SET_SIZE(nUnit,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
      do jc=1,NEnergySpec2D
        call MPI_FILE_SET_VIEW(nUnit,disp,data_type,newtype,'native',MPI_INFO_NULL,ierror)
#ifdef SAVE_SINGLE_Spec2D
        call MPI_FILE_WRITE_ALL(nUnit,EnergySpecOut(:,:,:,jc),subsizes(1)*subsizes(2)*subsizes(3),data_type,MPI_STATUS_IGNORE,ierror)
#else
        call MPI_FILE_WRITE_ALL(nUnit,EnergySpec2D(:,:,:,jc), subsizes(1)*subsizes(2)*subsizes(3),data_type,MPI_STATUS_IGNORE,ierror)
#endif
        disp=disp+ int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(data_byte,8) ! Update displacement
      enddo
      call MPI_FILE_CLOSE(nUnit,ierror)
      call MPI_TYPE_FREE(newtype,ierror)
      End Block
    ENDIF

    nfstime=0; SumStat=zero; PrGradsum=zero; nSpectime=0; 
    if(clcSpectra1D) then
      EnergySpecX=zero; EnergySpecZ=zero
    endif
    if(clcSpectra2D) EnergySpec2D=zero
  end subroutine clcStat

  !******************************************************************
  ! clcEnergySpectraX
  !******************************************************************
  subroutine clcEnergySpectraX(arrx,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(x1size(1),x1size(2),x1size(3)),intent(in)::arrx

    ! locals
    integer::ic,jc,kc,icCnter
    real(RK)::normEgy,EgyX(nxhp)

    normEgy= one/(real(nxc,RK)*real(nxc,RK)*real(nzc,RK))
    do jc=1,x1size(2)
      EgyX=zero
      do kc=1,x1size(3)
        EgyX(1)   = EgyX(1)   + arrx(1,jc,kc)   *arrx(1,jc,kc)
        EgyX(nxhp)= EgyX(nxhp)+ arrx(nxhp,jc,kc)*arrx(nxhp,jc,kc)
        do ic=2,nxh
          icCnter=nxc+2-ic
          EgyX(ic)= EgyX(ic) + (arrx(ic,jc,kc)*arrx(ic,jc,kc)+arrx(icCnter,jc,kc)*arrx(icCnter,jc,kc))*two
        enddo
      enddo
      do ic=1,nxhp
        EnergySpecX(ic,jc,m)= EnergySpecX(ic,jc,m) +EgyX(ic)*normEgy
      enddo
    enddo
  end subroutine clcEnergySpectraX

  !******************************************************************
  ! clcLTS_x_real
  !******************************************************************
  subroutine clcLTS_x_real(arrx,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(x1size(1),x1size(2),x1size(3)),intent(in)::arrx

    ! locals
    integer::ic,jc,kc,icCnter
    real(RK)::normEgy,EgyX(nxhp)

    normEgy= one/(real(nxc,RK)*real(nxc,RK)*real(nzc,RK))
    do jc=1,x1size(2)
      EgyX=zero
      do kc=1,x1size(3)
        EgyX(1)   = EgyX(1)   + arrx(1,jc,kc)
        EgyX(nxhp)= EgyX(nxhp)+ arrx(nxhp,jc,kc)
        do ic=2,nxh
          icCnter=nxc+2-ic
          EgyX(ic)= EgyX(ic) + (arrx(ic,jc,kc)+arrx(icCnter,jc,kc)) ! Note, there is NO "*two" here !
        enddo
      enddo
      do ic=1,nxhp
        EnergySpecX(ic,jc,m)= EnergySpecX(ic,jc,m) +EgyX(ic)*normEgy
      enddo
    enddo
  end subroutine clcLTS_x_real

  !******************************************************************
  ! clcLTS_x_imag
  !******************************************************************
  subroutine clcLTS_x_imag(arrx,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(x1size(1),x1size(2),x1size(3)),intent(in)::arrx

    ! locals
    integer::ic,jc,kc,icCnter
    real(RK)::normEgy,EgyX(nxhp)

    normEgy= one/(real(nxc,RK)*real(nxc,RK)*real(nzc,RK))
    do jc=1,x1size(2)
      EgyX=zero
      do kc=1,x1size(3)
        EgyX(1)   = zero
        EgyX(nxhp)= zero
        do ic=2,nxh
          icCnter=nxc+2-ic
          EgyX(ic)= EgyX(ic) + (arrx(ic,jc,kc)-arrx(icCnter,jc,kc)) ! Note here !
        enddo
      enddo
      do ic=1,nxhp
        EnergySpecX(ic,jc,m)= EnergySpecX(ic,jc,m) +EgyX(ic)*normEgy
      enddo
    enddo
  end subroutine clcLTS_x_imag
   
  !******************************************************************
  ! clcCospectraX
  !******************************************************************
  subroutine clcCospectraX(arrx1,arrx2,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(x1size(1),x1size(2),x1size(3)),intent(in)::arrx1,arrx2

    ! locals
    integer::ic,jc,kc,icCnter
    real(RK)::normEgy,EgyX(nxhp)

    normEgy= one/(real(nxc,RK)*real(nxc,RK)*real(nzc,RK))
    do jc=1,x1size(2)
      EgyX=zero
      do kc=1,x1size(3)
        EgyX(1)   = EgyX(1)   + arrx1(1,jc,kc)   *arrx2(1,jc,kc)
        EgyX(nxhp)= EgyX(nxhp)+ arrx1(nxhp,jc,kc)*arrx2(nxhp,jc,kc)
        do ic=2,nxh
          icCnter=nxc+2-ic
          EgyX(ic)= EgyX(ic) + (arrx1(ic,jc,kc)*arrx2(ic,jc,kc)+arrx1(icCnter,jc,kc)*arrx2(icCnter,jc,kc))*two
        enddo
      enddo
      do ic=1,nxhp
        EnergySpecX(ic,jc,m)= EnergySpecX(ic,jc,m) +EgyX(ic)*normEgy
      enddo
    enddo
  end subroutine clcCospectraX

  !******************************************************************
  ! clcEnergySpectraZ
  !******************************************************************
  subroutine clcEnergySpectraZ(arrz,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(z1size(1),z1size(2),z1size(3)),intent(in)::arrz

    ! locals
    integer::ic,jc,kc,kcCnter
    real(RK)::normEgy,EgyZ(nzhp)

    normEgy= one/(real(nxc,RK)*real(nzc,RK)*real(nzc,RK))
    do jc=1,z1size(2)
      EgyZ=zero
      do ic=1,z1size(1)
        EgyZ(1)   = EgyZ(1)   + arrz(ic,jc,1)   *arrz(ic,jc,1)
        EgyZ(nzhp)= EgyZ(nzhp)+ arrz(ic,jc,nzhp)*arrz(ic,jc,nzhp)
        do kc=2,nzh
          kcCnter=nzc+2-kc
          EgyZ(kc)= EgyZ(kc) + (arrz(ic,jc,kc)*arrz(ic,jc,kc)+arrz(ic,jc,kcCnter)*arrz(ic,jc,kcCnter))*two
        enddo
      enddo
      do kc=1,nzhp
        EnergySpecZ(kc,jc,m)= EnergySpecZ(kc,jc,m) +EgyZ(kc)*normEgy
      enddo
    enddo
  end subroutine clcEnergySpectraZ

  !******************************************************************
  ! clcLTS_z_real
  !******************************************************************
  subroutine clcLTS_z_real(arrz,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(z1size(1),z1size(2),z1size(3)),intent(in)::arrz

    ! locals
    integer::ic,jc,kc,kcCnter
    real(RK)::normEgy,EgyZ(nzhp)

    normEgy= one/(real(nxc,RK)*real(nzc,RK)*real(nzc,RK))
    do jc=1,z1size(2)
      EgyZ=zero
      do ic=1,z1size(1)
        EgyZ(1)   = EgyZ(1)   + arrz(ic,jc,1)   
        EgyZ(nzhp)= EgyZ(nzhp)+ arrz(ic,jc,nzhp)
        do kc=2,nzh
          kcCnter=nzc+2-kc
          EgyZ(kc)= EgyZ(kc) + (arrz(ic,jc,kc)+arrz(ic,jc,kcCnter))
        enddo
      enddo
      do kc=1,nzhp
        EnergySpecZ(kc,jc,m)= EnergySpecZ(kc,jc,m) +EgyZ(kc)*normEgy
      enddo
    enddo
  end subroutine clcLTS_z_real

  !******************************************************************
  ! clcLTS_z_imag
  !******************************************************************
  subroutine clcLTS_z_imag(arrz,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(z1size(1),z1size(2),z1size(3)),intent(in)::arrz

    ! locals
    integer::ic,jc,kc,kcCnter
    real(RK)::normEgy,EgyZ(nzhp)

    normEgy= one/(real(nxc,RK)*real(nzc,RK)*real(nzc,RK))
    do jc=1,z1size(2)
      EgyZ=zero
      do ic=1,z1size(1)
        EgyZ(1)   = EgyZ(1)   + zero   
        EgyZ(nzhp)= EgyZ(nzhp)+ zero
        do kc=2,nzh
          kcCnter=nzc+2-kc
          EgyZ(kc)= EgyZ(kc) + (arrz(ic,jc,kc)-arrz(ic,jc,kcCnter))
        enddo
      enddo
      do kc=1,nzhp
        EnergySpecZ(kc,jc,m)= EnergySpecZ(kc,jc,m) +EgyZ(kc)*normEgy
      enddo
    enddo
  end subroutine clcLTS_z_imag
    
  !******************************************************************
  ! clcCospectraZ
  !******************************************************************
  subroutine clcCospectraZ(arrz1,arrz2,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(z1size(1),z1size(2),z1size(3)),intent(in)::arrz1,arrz2

    ! locals
    integer::ic,jc,kc,kcCnter
    real(RK)::normEgy,EgyZ(nzhp)

    normEgy= one/(real(nxc,RK)*real(nzc,RK)*real(nzc,RK))
    do jc=1,z1size(2)
      EgyZ=zero
      do ic=1,z1size(1)
        EgyZ(1)   = EgyZ(1)   + arrz1(ic,jc,1)   *arrz2(ic,jc,1)
        EgyZ(nzhp)= EgyZ(nzhp)+ arrz1(ic,jc,nzhp)*arrz2(ic,jc,nzhp)
        do kc=2,nzh
          kcCnter=nzc+2-kc
          EgyZ(kc)= EgyZ(kc) + (arrz1(ic,jc,kc)*arrz2(ic,jc,kc)+arrz1(ic,jc,kcCnter)*arrz2(ic,jc,kcCnter))*two
        enddo
      enddo
      do kc=1,nzhp
        EnergySpecZ(kc,jc,m)= EnergySpecZ(kc,jc,m) +EgyZ(kc)*normEgy
      enddo
    enddo
  end subroutine clcCospectraZ
  
  !******************************************************************
  ! clcEnergySpectra2D
  !******************************************************************
  subroutine clcEnergySpectra2D(ArrIN,m,IsUy)
    implicit none
    integer,intent(in)::m
    logical,intent(in)::IsUy
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(in)::ArrIn
    
    ! locals
    integer::ic,jc,kc,jct,jcs
    real(RK)::normEgy,RealTemp1,RealTemp2
    real(RK),dimension(:,:,:),allocatable::arry1
#ifdef OverWriteFFT
    real(RK),dimension(:,:,:),allocatable::arrx
    real(RK),dimension(:,:,:),allocatable::arrz
#else
    real(RK),dimension(:,:,:),allocatable::arrx1,arrx2
    real(RK),dimension(:,:,:),allocatable::arrz1,arrz2
#endif
    
    normEgy= one/(real(nxc,RK)*real(nzc,RK))
    normEgy= normEgy*normEgy
    IF(FlowType==FT_CH)normEgy=normEgy*half

#ifdef OverWriteFFT    
    allocate(arrx(x1size(1),x1size(2),x1size(3)))
    allocate(arrz(z2size(1),z2size(2),z2size(3)))
    call transpose_y1_to_x1(ArrIN,arrx)
    call dfftw_execute_r2r(fft_plan_x, arrx, arrx)
    call transpose_x1_to_z2(arrx,arrz)
    deallocate(arrx)
    call dfftw_execute_r2r(fft_plan_z2,arrz,arrz)
    allocate(arry1(y2size(1),y2size(2),y2size(3)))
    call transpose_z2_to_y2(arrz,arry1)
    deallocate(arrz)
#else
    allocate(arrx1(x1size(1),x1size(2),x1size(3)), &
             arrx2(x1size(1),x1size(2),x1size(3))  )
    call transpose_y1_to_x1(ArrIN,arrx1)
    call dfftw_execute_r2r(fft_plan_x, arrx1,arrx2)
    deallocate(arrx1)
    allocate(arrz1(z2size(1),z2size(2),z2size(3)), &
             arrz2(z2size(1),z2size(2),z2size(3))  )  
    call transpose_x1_to_z2(arrx2,arrz1)
    call dfftw_execute_r2r(fft_plan_z2,arrz1,arrz2)
    deallocate(arrx2,arrz1)
    allocate(arry1(y2size(1),y2size(2),y2size(3)))
    call transpose_z2_to_y2(arrz2,arry1)
    deallocate(arrz2)
#endif

    IF(FlowType==FT_HC) THEN
      do jc=1,nySpec2D
        jct=iySpec2D(jc)
        do kc=1,y2size(3)
          do ic=1,y2size(1)
            RealTemp1=arry1(ic,jct,kc)
            EnergySpec2D(ic,jc,kc,m)=EnergySpec2D(ic,jc,kc,m)+normEgy*RealTemp1*RealTemp1
          enddo
        enddo
      enddo
    ELSE 
      do jc=1,nySpec2D
        jct=iySpec2D(jc)
        if(IsUy) then
          jcs=nyp+1-jct
        else
          jcs=nyc+1-jct        
        endif
        do kc=1,y2size(3)
          do ic=1,y2size(1)
            RealTemp1=arry1(ic,jct,kc)
            RealTemp2=arry1(ic,jcs,kc)
            EnergySpec2D(ic,jc,kc,m)=EnergySpec2D(ic,jc,kc,m)+normEgy*(RealTemp1*RealTemp1 +RealTemp2*RealTemp2)
          enddo
        enddo
      enddo
    ENDIF
  end subroutine clcEnergySpectra2D

  !******************************************************************
  ! clcCosSpectra2D
  !******************************************************************
  subroutine clcCosSpectra2D(ArrIN1,ArrIN2,m,IsSymmetry)
    implicit none
    integer,intent(in)::m
    logical,intent(in)::IsSymmetry
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(in)::ArrIn1,ArrIN2
    
    ! locals
    integer::ic,jc,kc,jct,jcs
    real(RK)::normEgy,RealTemp1,RealTemp2,rCoe
#ifdef OverWriteFFT
    real(RK),dimension(:,:,:),allocatable::arrx
    real(RK),dimension(:,:,:),allocatable::arrz
#else
    real(RK),dimension(:,:,:),allocatable::arrx1,arrx2
    real(RK),dimension(:,:,:),allocatable::arrz1,arrz2
#endif
    real(RK),dimension(:,:,:),allocatable::arry1,arry2

    
    normEgy= one/(real(nxc,RK)*real(nzc,RK))
    normEgy= normEgy*normEgy
    IF(FlowType==FT_CH)normEgy=normEgy*half
    
#ifdef OverWriteFFT
    ! arry1
    allocate(arrx(x1size(1),x1size(2),x1size(3)))
    allocate(arrz(z2size(1),z2size(2),z2size(3)))
    call transpose_y1_to_x1(ArrIN1,arrx)
    call dfftw_execute_r2r(fft_plan_x, arrx,arrx)
    call transpose_x1_to_z2(arrx,arrz)
    deallocate(arrx)
    call dfftw_execute_r2r(fft_plan_z2,arrz,arrz)
    allocate(arry1(y2size(1),y2size(2),y2size(3)))
    call transpose_z2_to_y2(arrz,arry1)
    deallocate(arrz)

    ! arry2
    allocate(arrx(x1size(1),x1size(2),x1size(3)))
    allocate(arrz(z2size(1),z2size(2),z2size(3)))
    call transpose_y1_to_x1(ArrIN2,arrx)
    call dfftw_execute_r2r(fft_plan_x, arrx,arrx)
    call transpose_x1_to_z2(arrx,arrz)
    deallocate(arrx)
    call dfftw_execute_r2r(fft_plan_z2,arrz,arrz)
    allocate(arry2(y2size(1),y2size(2),y2size(3)))
    call transpose_z2_to_y2(arrz,arry2)
    deallocate(arrz)
#else

    ! arry1
    allocate(arrx1(x1size(1),x1size(2),x1size(3)), &
             arrx2(x1size(1),x1size(2),x1size(3))  )
    call transpose_y1_to_x1(ArrIN1,arrx1)
    call dfftw_execute_r2r(fft_plan_x, arrx1,arrx2)
    deallocate(arrx1)
    allocate(arrz1(z2size(1),z2size(2),z2size(3)), &
             arrz2(z2size(1),z2size(2),z2size(3))  )  
    call transpose_x1_to_z2(arrx2,arrz1)
    call dfftw_execute_r2r(fft_plan_z2,arrz1,arrz2)
    deallocate(arrx2,arrz1)
    allocate(arry1(y2size(1),y2size(2),y2size(3)))
    call transpose_z2_to_y2(arrz2,arry1)
    deallocate(arrz2)

    ! arry2
    allocate(arrx1(x1size(1),x1size(2),x1size(3)), &
             arrx2(x1size(1),x1size(2),x1size(3))  )
    call transpose_y1_to_x1(ArrIN2,arrx1)
    call dfftw_execute_r2r(fft_plan_x, arrx1,arrx2)
    deallocate(arrx1)
    allocate(arrz1(z2size(1),z2size(2),z2size(3)), &
             arrz2(z2size(1),z2size(2),z2size(3))  )  
    call transpose_x1_to_z2(arrx2,arrz1)
    call dfftw_execute_r2r(fft_plan_z2,arrz1,arrz2)
    deallocate(arrx2,arrz1)
    allocate(arry2(y2size(1),y2size(2),y2size(3)))
    call transpose_z2_to_y2(arrz2,arry2)
    deallocate(arrz2)
#endif

    IF(FlowType==FT_HC) THEN
      do jc=1,nySpec2D
        jct=iySpec2D(jc)
        do kc=1,y2size(3)
          do ic=1,y2size(1)
            RealTemp1=arry1(ic,jct,kc)*arry2(ic,jct,kc)
            EnergySpec2D(ic,jc,kc,m)=EnergySpec2D(ic,jc,kc,m)+normEgy*RealTemp1
          enddo
        enddo
      enddo
    ELSE 
      if(IsSymmetry) then
        rCoe= one      
      else
        rCoe=-one
      endif
      do jc=1,nySpec2D
        jct=iySpec2D(jc)
        jcs=nyc+1-jct
        do kc=1,y2size(3)
          do ic=1,y2size(1)
            RealTemp1=arry1(ic,jct,kc)*arry2(ic,jct,kc)
            RealTemp2=arry1(ic,jcs,kc)*arry2(ic,jcs,kc)
            EnergySpec2D(ic,jc,kc,m)=EnergySpec2D(ic,jc,kc,m)+normEgy*(RealTemp1 +rCoe*RealTemp2)
          enddo
        enddo
      enddo
    ENDIF
  end subroutine clcCosSpectra2D  
#undef NCHASTAT
#undef NEnergySpec1D
#undef NEnergySpec2D
#undef SAVE_SINGLE_Spec2D
end module m_FlowCase
