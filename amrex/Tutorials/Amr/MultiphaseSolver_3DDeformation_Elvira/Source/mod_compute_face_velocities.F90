module mod_compute_face_velocities

  use amrex_base_module
  use amrex_amr_module
  use amr_data_module
  use fillpatch_module

  implicit none

contains

subroutine compute_face_velocities

    implicit none

    integer :: ilev
    integer :: borlo(4), borhi(4), fxlo(4), fxhi(4), fylo(4), fyhi(4), fzlo(4), fzhi(4)
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: Ufaceptr, Vfaceptr, Wfaceptr, phiborderptr
    integer :: nlevs
    integer, parameter :: ngrow = 3
    type(amrex_multifab) :: phiborder
    real(amrex_real) :: time
	
    time=0.0d0

    nlevs = amrex_get_finest_level()

    do ilev = 0, nlevs

       call amrex_multifab_build(phiborder, phi_new(ilev)%ba, phi_new(ilev)%dm, ncomp, ngrow)
       call fillpatch (ilev, time, phiborder, phi_mp, phi_mp)

       call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.true.)
       do while (mfi%next())
          bx = mfi%tilebox()
 	  phiborderptr => phiborder%dataptr(mfi)
	  Ufaceptr  => facevelmfab(1,ilev)%dataptr(mfi)
	  Vfaceptr  => facevelmfab(2,ilev)%dataptr(mfi)
	  Wfaceptr  => facevelmfab(3,ilev)%dataptr(mfi)
          borlo = lbound(phiborderptr)
          borhi = ubound(phiborderptr)
	  fxlo = lbound(Ufaceptr)
	  fxhi = ubound(Ufaceptr)
	  fylo = lbound(Vfaceptr)
	  fyhi = ubound(Vfaceptr)
	  fzlo = lbound(Wfaceptr)
	  fzhi = ubound(Wfaceptr)

          call interp_rho(ilev, bx%lo, bx%hi, Ufaceptr,&
  	       fxlo,fxhi,Vfaceptr,fylo,fyhi,Wfaceptr,fzlo,fzhi,phiborderptr(:,:,:,1),borlo,borhi,amrex_geom(ilev)%dx, amrex_problo)

       end do
       call amrex_mfiter_destroy(mfi)
     
       call amrex_multifab_destroy(phiborder)

    end do

end subroutine compute_face_velocities

subroutine face_velocities(level, lo, hi, RHOx, fxlo, fxhi, RHOy, fylo, fyhi, RHOz, fzlo, fzhi, RHO, borlo, borhi, dx, prob_lo)

  implicit none
  integer, intent(in) :: level, lo(3), hi(3), borlo(3), borhi(3), fxlo(3), fxhi(3), fylo(3), fyhi(3), fzlo(3), fzhi(3)
  double precision, dimension(borlo(1):borhi(1),borlo(2):borhi(2),borlo(3):borhi(3)) :: RHO
  double precision, dimension(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3)) :: RHOx
  double precision, dimension(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3)) :: RHOy
  double precision, dimension(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3)) :: RHOz

  double precision, intent(in) :: dx(3), prob_lo(3)
  integer          :: i,j,k
  double precision :: x,y,z

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)+1
    	   RHOx(i,j,k) = (RHO(i-1,j,k)+RHO(i,j,k))/2.0d0
        end do
     end do
  end do

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)+1
        do i=lo(1),hi(1)
    	   RHOy(i,j,k) = (RHO(i,j-1,k)+RHO(i,j,k))/2.0d0
        end do
     end do
  end do

  do k=lo(3),hi(3)+1
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
    	   RHOz(i,j,k) = (RHO(i,j,k-1)+RHO(i,j,k))/2.0d0
        end do
     end do
  end do

end subroutine face_velocities

end module mod_compute_face_velocities
