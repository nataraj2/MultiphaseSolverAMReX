program testirl
integer(kind=8) :: i, k

	do i=1, 1000
	  do j = 1, 1000
	   do k = 1, 1000
  	      call sub1
  	   end do
          end do
	end do
      

contains

subroutine sub1

use irl_fortran_interface
type(ObjServer_PlanarLoc_type) :: planar_localizer_allocation
integer(IRL_LargeOffsetIndex_t) :: heapsize
integer :: k

	   heapsize = int(10,8)*int(10,8)*int(10,8)
	   call new(planar_localizer_allocation,heapsize)

end subroutine sub1

end program testirl
