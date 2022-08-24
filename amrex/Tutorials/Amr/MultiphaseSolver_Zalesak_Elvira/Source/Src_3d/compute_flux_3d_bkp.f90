module compute_flux_module

  use amr_data_module
  use my_amr_module
  use compute_semilagrangian_flux
  use compgeom_lookup

  implicit none

  private

  public :: compute_flux_3d

contains

  subroutine compute_flux_3d(lo, hi, dt, dx,     &
                             phi,ph_lo,ph_hi,    &
                             umac,  u_lo,  u_hi, &
                             vmac,  v_lo,  v_hi, &
                             wmac,  w_lo,  w_hi, &
                             flxx, fx_lo, fx_hi, &
                             flxy, fy_lo, fy_hi, &
                             flxz, fz_lo, fz_hi, &
			     gradxphi_mat,gradyphi_mat,gradzphi_mat,&
			     flo,fhi,x,y,z,clo,chi,xm,ym,zm)

    integer, intent(in) :: lo(3), hi(3)
    double precision, intent(in) :: dt, dx(3)
    integer, intent(in) :: ph_lo(3), ph_hi(3)
    integer, intent(in) ::  u_lo(3),  u_hi(3)
    integer, intent(in) ::  v_lo(3),  v_hi(3)
    integer, intent(in) ::  w_lo(3),  w_hi(3)
    integer, intent(in) :: fx_lo(3), fx_hi(3)
    integer, intent(in) :: fy_lo(3), fy_hi(3)
    integer, intent(in) :: fz_lo(3), fz_hi(3)
    double precision, intent(in   ) :: phi (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3),ncomp)
    double precision, intent(in   ) :: umac( u_lo(1): u_hi(1), u_lo(2): u_hi(2), u_lo(3): u_hi(3))
    double precision, intent(in   ) :: vmac( v_lo(1): v_hi(1), v_lo(2): v_hi(2), v_lo(3): v_hi(3))
    double precision, intent(in   ) :: wmac( w_lo(1): w_hi(1), w_lo(2): w_hi(2), w_lo(3): w_hi(3))
    double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),ncomp)
    double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),ncomp)
    double precision, intent(  out) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),ncomp)
    double precision, intent(in   ) :: gradxphi_mat (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3),ncomp)
    double precision, intent(in   ) :: gradyphi_mat (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3),ncomp)
    double precision, intent(in   ) :: gradzphi_mat (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3),ncomp)
    integer, intent(in) :: flo(3),fhi(3),clo(3),chi(3)
    double precision, intent(in) :: x(flo(1):fhi(1))
    double precision, intent(in) :: y(flo(2):fhi(2))
    double precision, intent(in) :: z(flo(3):fhi(3))
    double precision, intent(in) :: xm(clo(1):chi(1))
    double precision, intent(in) :: ym(clo(2):chi(2))
    double precision, intent(in) :: zm(clo(3):chi(3))
    
    double precision, dimension(3,17) :: pt
    double precision, dimension(3,4) :: tet
    integer,  dimension(3,4) :: ind
    integer :: n
    double precision, dimension(13) :: flux

    integer :: i, j, k, icomp

    !!!!!!!!!!!!!!!!!!!!
    ! final edge states
    !!!!!!!!!!!!!!!!!!!!

    ! update phi on x faces by adding in yz and zy transverse terms
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)+1

    	   ! Construct and project full cell, pt(:,17) is for correction tets
           pt(:, 9)=[x(i),y(j)  ,z(k)  ]; pt(:,1)=pt(:,9)
           pt(:,10)=[x(i),y(j)  ,z(k+1)]; pt(:,2)=pt(:,10)
           pt(:,11)=[x(i),y(j+1),z(k+1)]; pt(:,3)=pt(:,11)
           pt(:,12)=[x(i),y(j+1),z(k)  ]; pt(:,4)=pt(:,12)
           pt(:,13)=[x(i),y(j)  ,z(k)  ]; pt(:,5)=project(pt(:,13),-dt,i,j,k,umac,u_lo,u_hi,&
	   								     vmac,v_lo,v_hi,&
                                                                             wmac,w_lo,w_hi,&
	   								     flo,fhi,x,y,z,&
	   								     clo,chi,xm,ym,zm)
           pt(:,14)=[x(i),y(j)  ,z(k+1)]; pt(:,6)=project(pt(:,14),-dt,i,j,k,umac,u_lo,u_hi,&
	   								      vmac,v_lo,v_hi,&
                                                                              wmac,w_lo,w_hi,&
	   								      flo,fhi,x,y,z,&
	   								      clo,chi,xm,ym,zm)
           pt(:,15)=[x(i),y(j+1),z(k+1)]; pt(:,7)=project(pt(:,15),-dt,i,j,k,umac,u_lo,u_hi,&
	   							              vmac,v_lo,v_hi,&		   
                                                                              wmac,w_lo,w_hi,&
	   								      flo,fhi,x,y,z,&
	   								      clo,chi,xm,ym,zm)
           pt(:,16)=[x(i),y(j+1)  ,z(k)]; pt(:,8)=project(pt(:,16),-dt,i,j,k,umac,u_lo,u_hi,&
	   								     vmac,v_lo,v_hi,&
                                                                             wmac,w_lo,w_hi,&
	    								     flo,fhi,x,y,z,&
	    								     clo,chi,xm,ym,zm)

          ! Advect the face
           flux=0.0d0
	   do n=1,6
          ! Build tet and corresponding indices

	      tet(:,1)=pt(:,tet_map2(1,n)); ind(:,1)=get_indices(tet(:,1),[i,j,k],flo,fhi,x,y,z)
              tet(:,2)=pt(:,tet_map2(2,n)); ind(:,2)=get_indices(tet(:,2),[i,j,k],flo,fhi,x,y,z)
              tet(:,3)=pt(:,tet_map2(3,n)); ind(:,3)=get_indices(tet(:,3),[i,j,k],flo,fhi,x,y,z)
              tet(:,4)=pt(:,tet_map2(4,n)); ind(:,4)=get_indices(tet(:,4),[i,j,k],flo,fhi,x,y,z)

              ! Add corresponding flux
              flux=flux+tet_sign(tet)*tet2flux(tet,ind,dx,ph_lo,ph_hi,phi,&
							  gradxphi_mat,gradyphi_mat,gradzphi_mat,&
							  flo,fhi,x,y,z,clo,chi,xm,ym,zm)
           end do
             ! compute final x-fluxes
	   do icomp=1,ncomp	
              flxx(i,j,k,icomp) = flux(icomp)/(dt*dx(2)*dx(3))
              !flxx(i,j,k,icomp) = 0.0d0
	   enddo		

          end do
       end do
    end do

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)+1
          do i = lo(1), hi(1)

    	   ! Construct and project full cell, pt(:,17) is for correction tets
           !pt(:, 9)=[x(i+1),y(j)  ,z(k)  ]; pt(:,1)=pt(:,9)
           !pt(:,10)=[x(i+1),y(j)  ,z(k+1)]; pt(:,2)=pt(:,10)
           !pt(:,11)=[x(i),y(j),z(k+1)]; pt(:,3)=pt(:,11)
           !pt(:,12)=[x(i),y(j),z(k)  ]; pt(:,4)=pt(:,12)

	   !pt(:,13)=[x(i+1),y(j)  ,z(k)  ]; pt(:,5)=project(pt(:,13),-dt,i,j,k,umac,u_lo,u_hi,&
           !                                                                 vmac,v_lo,v_hi,&
           !                                                                 wmac,w_lo,w_hi,&
           !                                                                 flo,fhi,x,y,z,&
           !                                                                 clo,chi,xm,ym,zm)
           !pt(:,14)=[x(i+1),y(j)  ,z(k+1)]; pt(:,6)=project(pt(:,14),-dt,i,j,k,umac,u_lo,u_hi,&
           !                                                                  vmac,v_lo,v_hi,&
           !                                                                  wmac,w_lo,w_hi,&
           !                                                                  flo,fhi,x,y,z,&
           !                                                                  clo,chi,xm,ym,zm)
           !pt(:,15)=[x(i),y(j),z(k+1)]; pt(:,7)=project(pt(:,15),-dt,i,j,k,umac,u_lo,u_hi,&
           !                                                                  vmac,v_lo,v_hi,&
           !                                                                  wmac,w_lo,w_hi,&
           !                                                                  flo,fhi,x,y,z,&
           !                                                                  clo,chi,xm,ym,zm)
           !pt(:,16)=[x(i),y(j)  ,z(k)]; pt(:,8)=project(pt(:,16),-dt,i,j,k,umac,u_lo,u_hi,&
           !                                                                 vmac,v_lo,v_hi,&
           !                                                                 wmac,w_lo,w_hi,&
           !                                                                 flo,fhi,x,y,z,&
           !                                                                 clo,chi,xm,ym,zm)

           !! Advect the face
           !flux=0.0d0
	   !do n=1,6
           !   ! Build tet and corresponding indices

	   !   tet(:,1)=pt(:,tet_map2(1,n)); ind(:,1)=get_indices(tet(:,1),[i,j,k],flo,fhi,x,y,z)
           !   tet(:,2)=pt(:,tet_map2(2,n)); ind(:,2)=get_indices(tet(:,2),[i,j,k],flo,fhi,x,y,z)
           !   tet(:,3)=pt(:,tet_map2(3,n)); ind(:,3)=get_indices(tet(:,3),[i,j,k],flo,fhi,x,y,z)
           !   tet(:,4)=pt(:,tet_map2(4,n)); ind(:,4)=get_indices(tet(:,4),[i,j,k],flo,fhi,x,y,z)

           !     ! Add corresponding flux
           !     flux=flux+tet_sign(tet)*tet2flux(tet,ind,dx,ph_lo,ph_hi,phi,&
	   !						     gradxphi_mat,gradyphi_mat,gradzphi_mat,&
	   !						     flo,fhi,x,y,z,clo,chi,xm,ym,zm)
           !end do
             ! compute final x-fluxes
	   do icomp=1,ncomp	
              !flxy(i,j,k,icomp) = flux(icomp)/(dt*dx(1)*dx(3))
              flxy(i,j,k,icomp) = 0.0d0!flux(icomp)/(dt*dx(1)*dx(3))
	   enddo		

	  end do
       end do
    end do

    do       k = lo(3), hi(3)+1
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)

          ! Construct and project full cell, pt(:,17) is for correction tets
           !pt(:, 9)=[x(i+1),y(j), z(k)]; pt(:,1)=pt(:,9)
           !pt(:,10)=[x(i),  y(j), z(k)]; pt(:,2)=pt(:,10)
           !pt(:,11)=[x(i),y(j+1),z(k)]; pt(:,3)=pt(:,11)
           !pt(:,12)=[x(i+1),y(j+1),z(k)]; pt(:,4)=pt(:,12)

           !pt(:,13)=[x(i+1), y(j), z(k)]; pt(:,5)=project(pt(:,13),-dt,i,j,k,umac,u_lo,u_hi,&
           !                                                                 vmac,v_lo,v_hi,&
           !                                                                 wmac,w_lo,w_hi,&
           !                                                                 flo,fhi,x,y,z,&
           !                                                                 clo,chi,xm,ym,zm)
           !pt(:,14)=[x(i), y(j) ,z(k)]; pt(:,6)=project(pt(:,14),-dt,i,j,k,umac,u_lo,u_hi,&
           !                                                                  vmac,v_lo,v_hi,&
           !                                                                  wmac,w_lo,w_hi,&
           !                                                                  flo,fhi,x,y,z,&
           !                                                                  clo,chi,xm,ym,zm)
           !pt(:,15)=[x(i),y(j+1),z(k)]; pt(:,7)=project(pt(:,15),-dt,i,j,k,umac,u_lo,u_hi,&
           !                                                                  vmac,v_lo,v_hi,&
           !                                                                  wmac,w_lo,w_hi,&
           !                                                                  flo,fhi,x,y,z,&
           !                                                                  clo,chi,xm,ym,zm)
           !pt(:,16)=[x(i+1),y(j+1),z(k)]; pt(:,8)=project(pt(:,16),-dt,i,j,k,umac,u_lo,u_hi,&
           !                                                                 vmac,v_lo,v_hi,&
           !                                                                 wmac,w_lo,w_hi,&
           !                                                                 flo,fhi,x,y,z,&
           !                                                                 clo,chi,xm,ym,zm)

          ! Advect the face
           !flux=0.0d0
           !do n=1,6
              ! Build tet and corresponding indices

           !   tet(:,1)=pt(:,tet_map2(1,n)); ind(:,1)=get_indices(tet(:,1),[i,j,k],flo,fhi,x,y,z)
           !   tet(:,2)=pt(:,tet_map2(2,n)); ind(:,2)=get_indices(tet(:,2),[i,j,k],flo,fhi,x,y,z)
           !   tet(:,3)=pt(:,tet_map2(3,n)); ind(:,3)=get_indices(tet(:,3),[i,j,k],flo,fhi,x,y,z)
           !   tet(:,4)=pt(:,tet_map2(4,n)); ind(:,4)=get_indices(tet(:,4),[i,j,k],flo,fhi,x,y,z)

           !     ! Add corresponding flux
           !     flux=flux+tet_sign(tet)*tet2flux(tet,ind,dx,ph_lo,ph_hi,phi,&
	   !					           gradxphi_mat,gradyphi_mat,gradzphi_mat,&
      	   !						   flo,fhi,x,y,z,clo,chi,xm,ym,zm)
           !end do
             ! compute final x-fluxes
           do icomp=1,ncomp
              flxz(i,j,k,icomp) = 0.0d0!flux(icomp)/(dt*dx(1)*dx(2))
           enddo

	  end do
       end do
    end do

  end subroutine compute_flux_3d

end module compute_flux_module
