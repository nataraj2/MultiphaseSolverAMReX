module tagging_module

  use iso_c_binding

  use amrex_fort_module
  use amr_data_module

  implicit none

  private

  public :: tag_phi_error

contains

  subroutine tag_phi_error (level, time, lo, hi, phi, philo, phihi, tag, taglo, taghi, phierr, &
       settag, cleartag,dx,prob_lo)
    integer, intent(in) :: level, lo(3), hi(3), philo(4), phihi(4), taglo(4), taghi(4)
    real(amrex_real) , intent(in   ) :: phi(philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3),ncomp)
    character(kind=c_char), intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(amrex_real), intent(in) :: time, phierr
    character(kind=c_char), intent(in) :: settag, cleartag

    integer :: i,j,k

    real(amrex_real) :: phierr_this
    real(amrex_real) :: dx(3), prob_lo(3), x, y, z, dircoord

    phierr_this = phierr

    if (level >= 1) then
       !
       !  This is here for testing only! 
       !  Remove this if you use this as a template.
       !
       if (time > 0.75_amrex_real .and. time < 1.0_amrex_real) then
          phierr_this = 2.0_amrex_real * phierr
       end if
    end if

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
		x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)
		y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
		z = prob_lo(3) + dx(3) * (dble(k)+0.5d0)
		dircoord = y
                !if (phi(i,j,k,1) .ge. 0.2d0 .and. phi(i,j,k,1).le.0.5d0) then
		if(dircoord .ge. 550.0d0*t_new(0)-0.1d0 .and. dircoord .le. 550.0d0*t_new(0)+0.1d0)then	
                !tag(i,j,k) = settag
		endif
		if(dircoord .ge. 300.0d0*t_new(0)-0.2d0 .and. dircoord .le. 300.0d0*t_new(0)+0.2d0)then	
                !tag(i,j,k) = settag
		endif
		if(dircoord .ge. -0.1d0 .and. dircoord .le. 0.1d0)then	
                !tag(i,j,k) = settag
		endif

		if(dsqrt(x**2+y**2).le.347.0d0*t_new(0)+2.0d0)then
		!tag(i,j,k) = settag
		endif

		if(dsqrt(y**2+z**2).le.347.0d0*t_new(0)+2.0d0)then
		!tag(i,j,k) = settag
		endif

		if(dsqrt(x**2+z**2).le.347.0d0*t_new(0)+2.0d0)then
		!tag(i,j,k) = settag
		endif
	
		if(dsqrt(x**2+y**2+z**2).le.347.0d0*t_new(0)+2.0d0)then
		!tag(i,j,k) = settag
		endif

                !if ((dsqrt(y**2+z**2).le.1.0d0.or.dsqrt(y**2+(z-x*5.0d0*t_new(0))**2).le.1.0d0) .and. x .le. 150.0d0*t_new(0))then
                if (dsqrt(y**2+z**2).ge.0.3d0.and.dsqrt(y**2+z**2).le.0.8d0 .and. x .le. 300.0d0*t_new(0))then
		tag(i,j,k) = settag
		endif



          enddo
       enddo
    enddo
  end subroutine tag_phi_error

end module tagging_module

