program main
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_Decomp2d
  use m_Parameters
  use m_MeshAndMetries
  use m_Interp
  implicit none
  integer::ierror
  type(MatBound)::mbOld
  character(len=64)::interpPrm
  real(RK),dimension(:,:,:),allocatable::ArrOld
  
  ! Initialize MPI
  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)

  ! Read parameters and Initialize Log file
  if(command_argument_count()/=1 .and. nrank==0) write(*,*)'command argument wrong!'
  call get_command_argument(1,interpPrm)
  call ReadAndInitParameters(interpPrm)

  if(nrank==0) call system('mkdir '//Res_Dir//' 2> /dev/null')
  call MainLog%InitLog(Res_Dir,RunName,LF_file_lvl,LF_cmdw_lvl)
  if(nrank==0) call DumpReadedParam()

  ! Initialize decomp2d
  call decomp2d_init(nxcOld,nycOld,nzcOld,p_row,p_col,BcOption)

  ! Read mesh
  call InitMeshAndMetries(interpPrm)

  ! Initialize interpolation
  call interp_init()

  ! Perform tri-linear interpolation
  mbOld%xme=1;  mbOld%xpe=1
  mbOld%yme=1;  mbOld%ype=1
  mbOld%zme=1;  mbOld%zpe=1
  call decomp2d_alloc(ArrOld,mbOld,opt_global=.true.)
  call Interp_uvwp(mbOld,ArrOld)   
  if(nrank==0)call MainLog%OutInfo("Good job! interpolateField finished successfully!",1)

  call MPI_FINALIZE(ierror)
end program main
