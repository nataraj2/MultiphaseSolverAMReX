module face_velocity_module

  implicit none

contains

  subroutine get_face_velocity(RHO,rhoU,rhoV,clo,chi, &
       vx, vxlo, vxhi, &
       vy, vylo, vyhi, &
       dx, prob_lo)
    use amrex_base_module

    integer, intent(in) :: vxlo(2), vxhi(2), vylo(2), vyhi(2)
    real(amrex_real), intent(out) :: vx(vxlo(1):vxhi(1),vxlo(2):vxhi(2))
    real(amrex_real), intent(out) :: vy(vylo(1):vyhi(1),vylo(2):vyhi(2))
    real(amrex_real), intent(in) :: dx(2), prob_lo(2)

    integer :: i, j
    real(amrex_real) :: x, y
    
    integer, intent(in) :: clo(3), chi(3)
    double precision, intent(in), dimension(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)) :: RHO,rhoU,rhoV

    ! x velocity
    do j = vxlo(2), vxhi(2)
       y = (dble(j)+0.5d0) * dx(2) + prob_lo(2)
       do i = vxlo(1), vxhi(1)
          x = dble(i) * dx(1) + prob_lo(1)
          vx(i,j) = (rhoU(i-1,j,1)+rhoU(i,j,1))/(RHO(i-1,j,1)+RHO(i,j,1))
       end do
    end do
    
    ! y velocity
    do j = vylo(2), vyhi(2)
       y = dble(j) * dx(2) + prob_lo(2)
       do i = vylo(1), vyhi(1)
          x = (dble(i)+0.5d0) * dx(1) + prob_lo(1)
          vy(i,j) = (rhoV(i,j-1,1)+rhoV(i,j,1))/(RHO(i,j-1,1)+RHO(i,j,1))
       end do
    end do

  end subroutine get_face_velocity

end module face_velocity_module
