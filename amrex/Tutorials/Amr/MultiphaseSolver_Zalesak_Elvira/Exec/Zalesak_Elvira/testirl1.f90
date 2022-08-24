program testirl
integer(kind=8) :: i

	do i=1, 1000
	print*, i
  	  call sub1
	end do

contains

subroutine sub1

use irl_fortran_interface
type(ObjServer_PlanarLoc_type) :: planar_localizer_allocation
integer(IRL_LargeOffsetIndex_t) :: heapsize
integer :: k

	do k=1,1000	
	   heapsize = int(10,8)*int(10,8)*int(10,8)
	   call new(planar_localizer_allocation,heapsize)
	enddo

end subroutine sub1

end program testirl
