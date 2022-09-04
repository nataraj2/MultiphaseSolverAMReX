module mod_compute_vortmag

  use amrex_base_module
  use amrex_amr_module
  use amr_data_module
  use fillpatch_module

  implicit none

contains

subroutine compute_vortmag(lev)

    implicit none

    integer, intent(in) :: lev
    integer :: borlo(4), borhi(4), clo(4), chi(4)
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: vortmagptr,phiborderptr
    integer :: nlevs
    integer, parameter :: ngrow = 3
    type(amrex_multifab) :: phiborder
    real(amrex_real) :: time
	
    time=0.0d0

       call amrex_multifab_build(phiborder, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, ngrow)

       call fillpatch(lev, time, phiborder,phi_new,phi_new)

       call amrex_mfiter_build(mfi, phi_new(lev), tiling=.true.)
       do while (mfi%next())
          bx = mfi%tilebox()
          vortmagptr      => vortmagmfab(lev)%dataptr(mfi)
	  phiborderptr => phiborder%dataptr(mfi)
	  borlo    = lbound(phiborderptr)
	  borhi    = ubound(phiborderptr)
	  clo    = lbound(vortmagptr)
	  chi    = ubound(vortmagptr)

          call vortmag(lev, bx%lo, bx%hi, phiborderptr, borlo, borhi, vortmagptr, clo, chi, amrex_geom(lev)%dx, amrex_problo)

       end do
       call amrex_mfiter_destroy(mfi)

       call amrex_multifab_destroy(phiborder)


end subroutine compute_vortmag

subroutine vortmag(level, lo, hi, phi, borlo, borhi, vortmagnitude, clo, chi, dx, prob_lo)

  implicit none
  integer, intent(in) :: level, lo(3), hi(3), borlo(3), borhi(3), clo(3), chi(3)
  double precision, dimension(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)) :: vortmagnitude
  double precision, dimension(borlo(1):borhi(1),borlo(2):borhi(2),borlo(3):borhi(3),ncomp) :: phi

  double precision, intent(in) :: dx(3), prob_lo(3)
  integer          :: i,j,k
  double precision :: dUdy,dUdz,dVdx,dVdz,dWdx,dWdy,x,y,z
  double precision, dimension(borlo(1):borhi(1),borlo(2):borhi(2),borlo(3):borhi(3)) :: U, V, W 

  do k=borlo(3),borhi(3)
     do j=borlo(2),borhi(2)
        do i=borlo(1),borhi(1)
   	   U(i,j,k)=((1.0d0-phi(i,j,k,11))*phi(i,j,k,2)+phi(i,j,k,11)*phi(i,j,k,7))/((1.0d0-phi(i,j,k,11))*phi(i,j,k,1)+phi(i,j,k,11)*phi(i,j,k,6))
   	   V(i,j,k)=((1.0d0-phi(i,j,k,11))*phi(i,j,k,3)+phi(i,j,k,11)*phi(i,j,k,8))/((1.0d0-phi(i,j,k,11))*phi(i,j,k,1)+phi(i,j,k,11)*phi(i,j,k,6)) 
   	   W(i,j,k)=((1.0d0-phi(i,j,k,11))*phi(i,j,k,4)+phi(i,j,k,11)*phi(i,j,k,9))/((1.0d0-phi(i,j,k,11))*phi(i,j,k,1)+phi(i,j,k,11)*phi(i,j,k,6))
  	end do
     end do
  end do	

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
	   x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)
           y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
           z = prob_lo(3) + dx(3) * (dble(k)+0.5d0)
		
           !if ((dsqrt(y**2+z**2).le.1.0d0.or.dsqrt(y**2+(z-x*5.0d0*t_new(0))**2).le.1.0d0) .and. x .le. 150.0d0*t_new(0))then
   	   !   vortmagnitude(i,j,k) = 1.0d0
           !else
	   !   vortmagnitude(i,j,k) = 0.0d0
	   !endif
	
  	   dUdy = (U(i,j+1,k)-U(i,j-1,k))/(2.0d0*dx(2))
    	   dUdz = (U(i,j,k+1)-U(i,j,k-1))/(2.0d0*dx(3))
	   dVdx = (V(i+1,j,k)-V(i-1,j,k))/(2.0d0*dx(1))
	   dVdz = (V(i,j,k+1)-V(i,j,k-1))/(2.0d0*dx(3))
	   dWdx = (W(i+1,j,k)-W(i-1,j,k))/(2.0d0*dx(1))
	   dWdy = (W(i,j+1,k)-W(i,j-1,k))/(2.0d0*dx(2))
	   vortmagnitude(i,j,k) = dsqrt((dWdy-dVdz)**2+(dUdz-dWdx)**2+(dVdx-dUdy)**2) 

        end do
     end do
  end do



end subroutine vortmag
 
end module mod_compute_vortmag
