module dump_ensight_str
  use precision
  use amrex_fi_mpi
  implicit none
  
  ! Time info
  integer :: nout_time
  !real(WP), dimension(:), allocatable :: out_time
  real(WP), dimension(1) :: out_time
  
  ! SP buffers
  real(SP), dimension(:,:,:), allocatable :: buffer1_SP
  real(SP), dimension(:,:,:,:), allocatable :: buffer3_SP
  
  ! Fileviews
  integer :: fileview,datasize,gdatasize
  
  ! Monitoring
  integer  :: nfiles
  real(WP) :: time_open,time_close,time_write

  !Interface
  integer :: nNodes,nZones
  integer,  dimension(:), allocatable :: nodeList
  real, dimension(:), allocatable :: xNode,yNode,zNode

  ! For curvature calculation
  integer, parameter :: npnorm=3
  double precision, dimension(3,npnorm) :: pnorm
  data pnorm(:,1) /1,0,0/
  data pnorm(:,2) /0,1,0/
  data pnorm(:,3) /0,0,1/

  character(len=300) :: plicdirname,filename
  integer :: nAllocated,nListAllocated
  real, dimension(:), allocatable :: tmpNode
  integer,  dimension(:), allocatable :: tmpNodeList
end module dump_ensight_str

! ======================================= !
! Dump 3D binary ensight gold data - PLIC !
! This uses VOF2mesh which is in          !
! compressible_plic.f90                   !
! ======================================= !
subroutine dump_ensight_str_plic(prob_lo,prob_hi,step)
  use fileio
  use dump_ensight_str
  !use plic_module
  use string
  
  implicit none
  
  integer :: iunit,ierr,i,rank
  character(len=str_medium) :: file
  character(len=80) :: cbuffer,str
  integer :: ibuffer
  real(SP) :: rbuffer
 
  double precision :: prob_lo(3), prob_hi(3)
  integer :: irank, iroot, nproc

  integer :: dummyint
  integer, intent(in) :: step

  call MPI_Comm_Rank(MPI_COMM_WORLD, irank, ierr) 
  call MPI_Comm_Size(MPI_COMM_WORLD, nproc, ierr)

  !!!!!!!!!!!!! 
  irank=irank+1
  iroot=1
  !!!!!!!!!!!!
  
  ! Generate the geometry file - PLIC
  file = "ensight-3D/" // trim(plicdirname) // "/" // 'plic' // "."
  write(file(len_trim(file)+1:len_trim(file)+6),'(i6.6)') step
 
 
  ! Write node x position
  do rank=1,nproc
     if (rank.eq.irank) then
           call BINARY_FILE_OPEN(iunit,trim(file),"a",ierr)
           ! Part header
           cbuffer = 'part'
           call BINARY_FILE_WRITE_CHAR(iunit,cbuffer,80,kind(cbuffer),ierr)
           ibuffer = irank
           call BINARY_FILE_WRITE_INT(iunit,ibuffer,1,kind(ibuffer),ierr)
           cbuffer = 'NGA 3D PLIC geometry per processor'
           call BINARY_FILE_WRITE_CHAR(iunit,cbuffer,80,kind(cbuffer),ierr)
           ! Position of my nodes
           cbuffer = 'coordinates'
           call BINARY_FILE_WRITE_CHAR(iunit,cbuffer,80,kind(cbuffer),ierr) 
           ibuffer = nNodes
           call BINARY_FILE_WRITE_INT(iunit,ibuffer,1,kind(ibuffer),ierr)
           call BINARY_FILE_WRITE_FLOAT(iunit,xNode(1:nNodes),nNodes,kind(xNode),ierr)
           call BINARY_FILE_WRITE_FLOAT(iunit,yNode(1:nNodes),nNodes,kind(yNode),ierr)
           call BINARY_FILE_WRITE_FLOAT(iunit,zNode(1:nNodes),nNodes,kind(zNode),ierr)
           ! Write header
           cbuffer = 'tria3'
           call BINARY_FILE_WRITE_CHAR(iunit,cbuffer,80,kind(cbuffer),ierr)
           ibuffer = nZones
           call BINARY_FILE_WRITE_INT(iunit,ibuffer,1,kind(ibuffer),ierr)
           ! Write my triangles
           call BINARY_FILE_WRITE(iunit,nodeList(1:nZones),3*nZones,kind(nodeList),ierr)
           ! Close the file
           call BINARY_FILE_CLOSE(iunit,ierr)
     	end if
           call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   end do
  ! Deallocate
  deallocate(xNode,yNode,zNode,nodeList)
  
  return
end subroutine dump_ensight_str_plic


subroutine dump_ensight_str_plic_headers(prob_lo,prob_hi,step)
  use fileio
  use dump_ensight_str
  use string
  use my_amr_module
 
  implicit none
  
  integer :: iunit,ierr,i,rank
  character(len=str_medium) :: file
  character(len=80) :: cbuffer,str
  integer :: ibuffer
  real(SP), allocatable :: rbuffer(:)
 
  double precision :: prob_lo(3), prob_hi(3)
  integer :: irank, iroot, nproc

  integer, intent (in) :: step
  integer :: dummyint
  character(len=20) :: stepstr

  nout_time=1
  out_time=0.0d0
  	
  call MPI_Comm_Rank(MPI_COMM_WORLD, irank, ierr) 
  call MPI_Comm_Size(MPI_COMM_WORLD, nproc, ierr)
 

  !!!!!!!!!!!!! 
  irank=irank+1
  iroot=1
  !!!!!!!!!!!!

  ! Generate the geometry file - PLIC
     write(plicdirname,'(a,i5.5,a)') "plic",step
  if (irank.eq.iroot) then
     call CREATE_FOLDER("ensight-3D/" // trim(plicdirname))
  endif
  
  ! Write - Single proc & parallel => only root writes (in ASCII)
  if (irank.eq.iroot) then
     ! Open the file
     iunit=iopen()
     filename="ensight-3D/" // trim(plicdirname) // "/" // 'plic.case'
     open(iunit,file=filename,form="formatted",iostat=ierr,status="REPLACE")
     ! Write the case
     str='FORMAT'
     write(iunit,'(a80)') str
     str='type: ensight gold'
     write(iunit,'(a80)') str
     str='GEOMETRY'
     write(iunit,'(a80)') str
     str='model: 1 plic.******'
     write(iunit,'(a80)') str
     ! Time section
     str='TIME'
     write(iunit,'(a80)') str
     str='time set: 1'
     write(iunit,'(a80)') str
     str='number of steps:'
     write(iunit,'(a16,x,i12)') str,nout_time
     write (stepstr, *) step
     str='filename start number: '//trim(stepstr)
     write(iunit,'(a80)') str
     str='filename increment: 1'
     write(iunit,'(a80)') str
     str='time values:'
     write(iunit,'(a12,x,10000000(3(ES12.5,x),/))') str,t_new(0)
     ! Close the file
     close(iclose(iunit))
  end if
 
 
  ! Generate the geometry file - PLIC
  !if (irank.eq.iroot) call CREATE_FOLDER("ensight-3D/" // plicdirname)
  file = "ensight-3D/" // trim(plicdirname) // "/" // 'plic' // "."

  write(file(len_trim(file)+1:len_trim(file)+6),'(i6.6)') step
  
    ! Write node x position
        if (irank.eq.1) then
           ! Open the file
           call BINARY_FILE_OPEN(iunit,trim(file),"w",ierr)
           cbuffer = 'C Binary'
           call BINARY_FILE_WRITE_CHAR(iunit,cbuffer,80,kind(cbuffer),ierr)
           cbuffer = 'Ensight Gold Geometry File'
           call BINARY_FILE_WRITE_CHAR(iunit,cbuffer,80,kind(cbuffer),ierr)
           cbuffer = 'PLIC Mesh from NGA'
           call BINARY_FILE_WRITE_CHAR(iunit,cbuffer,80,kind(cbuffer),ierr)
           cbuffer = 'node id off'
           call BINARY_FILE_WRITE_CHAR(iunit,cbuffer,80,kind(cbuffer),ierr)
           cbuffer = 'element id off'
           call BINARY_FILE_WRITE_CHAR(iunit,cbuffer,80,kind(cbuffer),ierr)
           cbuffer = 'extents'
           call BINARY_FILE_WRITE_CHAR(iunit,cbuffer,80,kind(cbuffer),ierr)
           allocate(rbuffer(1))
	   !!!!!!!!!!!!!!!
           rbuffer(1) = prob_lo(1); call BINARY_FILE_WRITE_FLOAT(iunit,rbuffer,1,kind(rbuffer),ierr)
           rbuffer = prob_hi(1); call BINARY_FILE_WRITE_FLOAT(iunit,rbuffer,1,kind(rbuffer),ierr)
           rbuffer = prob_lo(2); call BINARY_FILE_WRITE_FLOAT(iunit,rbuffer,1,kind(rbuffer),ierr)
           rbuffer = prob_hi(2); call BINARY_FILE_WRITE_FLOAT(iunit,rbuffer,1,kind(rbuffer),ierr)
           rbuffer = prob_lo(3); call BINARY_FILE_WRITE_FLOAT(iunit,rbuffer,1,kind(rbuffer),ierr)
           rbuffer = prob_hi(3); call BINARY_FILE_WRITE_FLOAT(iunit,rbuffer,1,kind(rbuffer),ierr)
		
           call BINARY_FILE_CLOSE(iunit,ierr)
      end if

  return
end subroutine dump_ensight_str_plic_headers
