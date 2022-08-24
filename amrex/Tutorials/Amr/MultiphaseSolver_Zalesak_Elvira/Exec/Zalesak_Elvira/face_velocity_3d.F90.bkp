module face_velocity_module

  implicit none

contains

subroutine get_face_velocity(RHO,rhoU,rhoV,rhoW,clo,chi, &
     vx, vxlo, vxhi, &
     vy, vylo, vyhi, &
#if BL_SPACEDIM == 3
     vz, vzlo, vzhi, &
#endif
     dx, prob_lo) bind(C, name="get_face_velocity")

use amrex_base_module

  implicit none

  integer, intent(in) :: vxlo(3),vxhi(3)
  double precision, intent(out) :: vx(vxlo(1):vxhi(1),vxlo(2):vxhi(2),vxlo(3):vxhi(3))
  integer, intent(in) :: vylo(3),vyhi(3)
  double precision, intent(out) :: vy(vylo(1):vyhi(1),vylo(2):vyhi(2),vylo(3):vyhi(3))
#if BL_SPACEDIM == 3
  integer, intent(in) :: vzlo(3),vzhi(3)
  double precision, intent(out) :: vz(vzlo(1):vzhi(1),vzlo(2):vzhi(2),vzlo(3):vzhi(3))
#endif
  double precision, intent(in) :: dx(3), prob_lo(3)

  integer :: i, j, k
  double precision :: x, y, z
  integer :: vx_l1,vx_l2,vx_l3,vy_l1,vy_l2,vy_l3,vz_l1,vz_l2,vz_l3
  integer :: vx_h1,vx_h2,vx_h3,vy_h1,vy_h2,vy_h3,vz_h1,vz_h2,vz_h3

  integer, intent(in) :: clo(3), chi(3)
  double precision, intent(in), dimension(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)) :: RHO,rhoU,rhoV,rhoW

  vx_l1=vxlo(1);vx_l2=vxlo(2);vx_l3=vxlo(3)
  vx_h1=vxhi(1);vx_h2=vxhi(2);vx_h3=vxhi(3)
  vy_l1=vylo(1);vy_l2=vylo(2);vy_l3=vylo(3)
  vy_h1=vyhi(1);vy_h2=vyhi(2);vy_h3=vyhi(3)
#if BL_SPACEDIM == 3
  vz_l1=vzlo(1);vz_l2=vzlo(2);vz_l3=vzlo(3) 
  vz_h1=vzhi(1);vz_h2=vzhi(2);vz_h3=vzhi(3)
#endif

  ! x velocity

  do k = vx_l3, vx_h3
     z = (dble(k)+0.5d0) * dx(3) + prob_lo(3)
  do j = vx_l2, vx_h2
     y = (dble(j)+0.5d0) * dx(2) + prob_lo(2)
     do i = vx_l1, vx_h1
        x = dble(i) * dx(1) + prob_lo(1)
#if BL_SPACEDIM == 2
        vx(i,j,k) = (rhoU(i-1,j,1)+rhoU(i,j,1))/(RHO(i-1,j,1)+RHO(i,j,1))
#endif
#if BL_SPACEDIM == 3 
        vx(i,j,k) = (rhoU(i-1,j,k)+rhoU(i,j,k))/(RHO(i-1,j,k)+RHO(i,j,k))
#endif
     end do
  end do
  end do

  ! y velocity
  do k = vy_l3, vy_h3
  do j = vy_l2, vy_h2
     y = dble(j) * dx(2) + prob_lo(2)
     do i = vy_l1, vy_h1
        x = (dble(i)+0.5d0) * dx(1) + prob_lo(1)
#if BL_SPACEDIM == 2
        vy(i,j,k) = (rhoV(i,j-1,1)+rhoV(i,j,1))/(RHO(i,j-1,1)+RHO(i,j,1)) 
#endif
#if BL_SPACEDIM == 3
        vy(i,j,k) = (rhoV(i,j-1,k)+rhoV(i,j,k))/(RHO(i,j-1,k)+RHO(i,j,k)) 
#endif

     end do
  end do
  end do

#if BL_SPACEDIM == 3
 do k = vz_l3, vz_h3
  do j = vz_l2, vz_h2
     y = dble(j) * dx(2) + prob_lo(2)
     do i = vz_l1, vz_h1
	  vz(i,j,k) = (rhoW(i,j,k-1)+rhoW(i,j,k))/(RHO(i,j,k-1)+RHO(i,j,k))
     enddo
  enddo
 enddo	
#endif

end subroutine get_face_velocity

end module face_velocity_module
