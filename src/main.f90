program main_channel3d
  use MPI
  use m_decomp2d
  use m_ChannelSystem
  use m_Parameters
  use m_Variables
  use m_IOAndVisu
  use m_LogInfo
  use m_Poisson,only:Destory_Poisson_FFT_Plan
  implicit none
  character(len=64)::ChannelPrm
  integer:: ierror,BcOption(6)

  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)
   
  ! read Channel options
  if(command_argument_count()/=1 .and. nrank==0) write(*,*)'command argument wrong!'
  call get_command_argument(1,ChannelPrm)
  call ReadAndInitParameters(ChannelPrm)

  if(FlowType==FT_CH) then
    BcOption=(/0,0,-1,-1,0,0/)
  elseif(FlowType==FT_HC) then
    BcOption=(/0,0,-1,-2,0,0/)
  endif
  call decomp_2d_init(nxc,nyc,nzc,nproc,p_row,p_col,y_pencil,BcOption)
      
  call ChannelInitialize(ChannelPrm)    ! Topest level initialing for Channel body
  do itime=ifirst, ilast
    call ChannelIterate()
  enddo
  if(nrank==0)call MainLog%OutInfo("Good job! Channel3d_4th finished successfully!",1)

  call Destory_Poisson_FFT_Plan()
  call decomp_2d_finalize()
  call MPI_FINALIZE(ierror)
end program main_channel3d
