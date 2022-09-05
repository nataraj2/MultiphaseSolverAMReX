module face_velocity_module

  implicit none

contains

subroutine get_face_velocity(Grho,GrhoU,GrhoV,GrhoW,Lrho,LrhoU,LrhoV,LrhoW,VOF,clo,chi, &
     vx, vxlo, vxhi, &
     vy, vylo, vyhi, &
     vz, vzlo, vzhi, &
     Gvol, Lvol, vollo, volhi,&
     dx, prob_lo) bind(C, name="get_face_velocity")

use amrex_base_module

  implicit none

  integer, intent(in) :: vxlo(3),vxhi(3),vollo(3),volhi(3)
  double precision, intent(out) :: vx(vxlo(1):vxhi(1),vxlo(2):vxhi(2),vxlo(3):vxhi(3))
  integer, intent(in) :: vylo(3),vyhi(3)
  double precision, intent(out) :: vy(vylo(1):vyhi(1),vylo(2):vyhi(2),vylo(3):vyhi(3))
  integer, intent(in) :: vzlo(3),vzhi(3)
  double precision, intent(out) :: vz(vzlo(1):vzhi(1),vzlo(2):vzhi(2),vzlo(3):vzhi(3))
  double precision, intent(in) :: Gvol(vollo(1):volhi(1),vollo(2):volhi(2),vollo(3):volhi(3),6)
  double precision, intent(in) :: Lvol(vollo(1):volhi(1),vollo(2):volhi(2),vollo(3):volhi(3),6)
  double precision, intent(in) :: dx(3), prob_lo(3)

  integer :: i, j, k
  double precision :: x, y, z
  integer :: vx_l1,vx_l2,vx_l3,vy_l1,vy_l2,vy_l3,vz_l1,vz_l2,vz_l3
  integer :: vx_h1,vx_h2,vx_h3,vy_h1,vy_h2,vy_h3,vz_h1,vz_h2,vz_h3

  integer, intent(in) :: clo(3), chi(3)
  double precision, dimension(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)) :: Grho,GrhoU,GrhoV,GrhoW,Lrho,LrhoU,LrhoV,LrhoW!,VOF
  double precision, dimension(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)) :: VOF
  double precision :: rho_l,rho_r,vol_l,vol_r
  double precision :: U_l,U_r,V_l,V_r,W_l,W_r

  vx_l1=vxlo(1);vx_l2=vxlo(2);vx_l3=vxlo(3)
  vx_h1=vxhi(1);vx_h2=vxhi(2);vx_h3=vxhi(3)
  vy_l1=vylo(1);vy_l2=vylo(2);vy_l3=vylo(3)
  vy_h1=vyhi(1);vy_h2=vyhi(2);vy_h3=vyhi(3)
  vz_l1=vzlo(1);vz_l2=vzlo(2);vz_l3=vzlo(3) 
  vz_h1=vzhi(1);vz_h2=vzhi(2);vz_h3=vzhi(3)

  ! x velocity
  do k = vx_l3, vx_h3
     do j = vx_l2, vx_h2
        do i = vx_l1, vx_h1

		!if(prob_lo(1)+dble(i)*dx(1).eq.0.0d0)then
		!	GrhoU(i-1,j,k)=0.0d0
		!endif
   	   U_l=((1.0d0-VOF(i-1,j,k))*GrhoU(i-1,j,k)+VOF(i-1,j,k)*LrhoU(i-1,j,k))/((1.0d0-VOF(i-1,j,k))*Grho(i-1,j,k)+VOF(i-1,j,k)*Lrho(i-1,j,k))
   	   U_r=((1.0d0-VOF(i  ,j,k))*GrhoU(i  ,j,k)+VOF(i  ,j,k)*LrhoU(i  ,j,k))/((1.0d0-VOF(i  ,j,k))*Grho(i  ,j,k)+VOF(i  ,j,k)*Lrho(i  ,j,k))
	   vol_r=Gvol(i,j,k,1)+Lvol(i,j,k,1)    ; rho_r=Gvol(i,j,k,1)*Grho(i,j,k)+Lvol(i,j,k,1)*Lrho(i,j,k)
           vol_l=Gvol(i-1,j,k,2)+Lvol(i-1,j,k,2); rho_l=Gvol(i-1,j,k,2)*Grho(i-1,j,k)+Lvol(i-1,j,k,2)*Lrho(i-1,j,k)
	   vx(i,j,k)=(rho_l*U_l+rho_r*U_r)/(rho_l+rho_r)
        end do
     end do
  end do

  ! y velocity
  do k = vy_l3, vy_h3
     do j = vy_l2, vy_h2
        do i = vy_l1, vy_h1

	   V_l=((1.0d0-VOF(i,j-1,k))*GrhoV(i,j-1,k)+VOF(i,j-1,k)*LrhoV(i,j-1,k))/((1.0d0-VOF(i,j-1,k))*Grho(i,j-1,k)+VOF(i,j-1,k)*Lrho(i,j-1,k))
   	   V_r=((1.0d0-VOF(i,j  ,k))*GrhoV(i,j  ,k)+VOF(i,j  ,k)*LrhoV(i,j  ,k))/((1.0d0-VOF(i,j  ,k))*Grho(i,j  ,k)+VOF(i,j  ,k)*Lrho(i,j  ,k))
	   vol_r=Gvol(i,j,k,3)+Lvol(i,j,k,3)    ; rho_r=Gvol(i,j,k,3)*Grho(i,j,k)+Lvol(i,j,k,3)*Lrho(i,j  ,k)
           vol_l=Gvol(i,j-1,k,4)+Lvol(i,j-1,k,4); rho_l=Gvol(i,j-1,k,4)*Grho(i,j-1,k)+Lvol(i,j-1,k,4)*Lrho(i,j-1,k)
	   vy(i,j,k)=(rho_l*V_l+rho_r*V_r)/(rho_l+rho_r)	
	  !y=prob_lo(2)+dble(j)*dx(2)
        end do
     end do
  end do

 ! z velocity
 do k = vz_l3, vz_h3
    do j = vz_l2, vz_h2
       do i = vz_l1, vz_h1

	  W_l=((1.0d0-VOF(i,j,k-1))*GrhoW(i,j,k-1)+VOF(i,j,k-1)*LrhoW(i,j,k-1))/((1.0d0-VOF(i,j,k-1))*Grho(i,j,k-1)+VOF(i,j,k-1)*Lrho(i,j,k-1))
   	  W_r=((1.0d0-VOF(i,j,k  ))*GrhoW(i,j,k  )+VOF(i,j,k  )*LrhoW(i,j,k  ))/((1.0d0-VOF(i,j,k  ))*Grho(i,j,k  )+VOF(i,j,k  )*Lrho(i,j,k  ))
	  vol_r=Gvol(i,j,k,5)+Lvol(i,j,k,5)    ; rho_r=Gvol(i,j,k,5)*Grho(i,j,k)+Lvol(i,j,k,5)*Lrho(i,j,k)
          vol_l=Gvol(i,j,k-1,6)+Lvol(i,j,k-1,6); rho_l=Gvol(i,j,k-1,6)*Grho(i,j,k-1)+Lvol(i,j,k-1,6)*Lrho(i,j,k-1)
	  vz(i,j,k)=(rho_l*W_l+rho_r*W_r)/(rho_l+rho_r)
	  if(prob_lo(3)+dble(k)*dx(3).le.-0.0025d0)then
                vz(i,j,k)=700.0d0
           endif
       enddo
    enddo
 enddo	

end subroutine get_face_velocity

end module face_velocity_module
