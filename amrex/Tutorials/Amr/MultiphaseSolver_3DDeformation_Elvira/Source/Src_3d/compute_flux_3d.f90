module compute_flux_module

  use amr_data_module
  use my_amr_module
  use compute_semilagrangian_flux
  use compgeom_lookup
  use irl_fortran_interface

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
			     volx, voly, volz, &
			     gradxphi_mat,gradyphi_mat,gradzphi_mat,&
			     plocalizer,pseparator,plocseplink,ploclink,&
			     x,y,z,xm,ym,zm)

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
    double precision, intent(out)  :: volx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3))
    double precision, intent(out)  :: voly(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3))
    double precision, intent(out)  :: volz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1)
    double precision, intent(in   ) :: gradxphi_mat (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3),ncomp)
    double precision, intent(in   ) :: gradyphi_mat (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3),ncomp)
    double precision, intent(in   ) :: gradzphi_mat (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3),ncomp)
    double precision, intent(in) :: x(lo(1)-1:hi(1)+2)
    double precision, intent(in) :: y(lo(2)-1:hi(2)+2)
    double precision, intent(in) :: z(lo(3)-1:hi(3)+2)
    double precision, intent(in) :: xm(lo(1)-1:hi(1)+1)
    double precision, intent(in) :: ym(lo(2)-1:hi(2)+1)
    double precision, intent(in) :: zm(lo(3)-1:hi(3)+1)   

    type(PlanarLoc_type), intent(in)        :: plocalizer(ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3))
    type(PlanarSep_type), intent(in)        :: pseparator(ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3))
    type(LocSepLink_type), intent(inout)    :: plocseplink(ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3))
    type(LocLink_type), intent(inout) 	    :: ploclink(ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3))

    double precision, dimension(3,4) :: tet
    integer,  dimension(3,4) :: ind
    integer :: n
    double precision, dimension(11) :: flux
    double precision :: volume

    integer :: i, j, k, icomp, ipt
    double precision, dimension(3,8) :: original_cell
    real(IRL_double), dimension(3,9) :: projected_cell, transported_cell
    integer,  dimension(3) :: lower_ind, upper_ind
 
  lower_ind(1)=lo(1)-1;lower_ind(2)=lo(2)-1;lower_ind(3)=lo(3)-1
  upper_ind(1)=hi(1)+1;upper_ind(2)=hi(2)+1;upper_ind(3)=hi(3)+1
  call setUpMeshConnectivity(lower_ind, upper_ind,ph_lo,ph_hi,plocseplink,ploclink)

    flux=0.0d0
    volume=0.0d0

  ! Allocate flux_hexahedron in IRL

   call new(flux_hexahedron(1))
   call new(f_moments)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)

	  ! Cell nodal points
	  original_cell(:,1)=[x(i),y(j),    z(k)  ] ;original_cell(:,2)=[x(i+1),y(j),    z(k)];
	  original_cell(:,3)=[x(i+1),y(j),z(k+1) ]  ;original_cell(:,4)=[x(i),y(j),    z(k+1)];
	  original_cell(:,5)=[x(i),y(j+1),  z(k)  ] ;original_cell(:,6)=[x(i+1),y(j+1),  z(k)];
	  original_cell(:,7)=[x(i+1),y(j+1), z(k+1)];original_cell(:,8)=[x(i),y(j+1),  z(k+1)];
	 
    	  ! Project all the cell nodal points just once
	  do ipt=1,8 
          projected_cell(:,ipt) =project(original_cell(:,ipt),-dt,i,j,k,umac,u_lo,u_hi,vmac,v_lo,v_hi,wmac,w_lo,w_hi,lo,hi,x,y,z,xm,ym,zm)
	  enddo
         
	  ! Compute x face flux 
	  transported_cell(:,1)=original_cell(:,4);transported_cell(:,2)=original_cell(:,1)
	  transported_cell(:,3)=original_cell(:,5);transported_cell(:,4)=original_cell(:,8)
	  transported_cell(:,5)=projected_cell(:,4);transported_cell(:,6)=projected_cell(:,1)
	  transported_cell(:,7)=projected_cell(:,5);transported_cell(:,8)=projected_cell(:,8)
	
	  transported_cell(:,9)=0.25d0*[sum(transported_cell(1,1:4)),sum(transported_cell(2,1:4)),sum(transported_cell(3,1:4))]
          transported_cell(:,9)=project(transported_cell(:,9),-dt,i,j,k,umac,u_lo,u_hi,vmac,v_lo,v_hi,wmac,w_lo,w_hi,lo,hi,x,y,z,xm,ym,zm)  

	  call construct(flux_hexahedron(1),transported_cell)
          call adjustCapToMatchVolume(flux_hexahedron(1),dt*umac(i,j,k)*dx(2)*dx(3))
          call getNormMoments(flux_hexahedron(1),plocseplink(i,j,k),f_moments)
          call compute_fluxes(i,j,k,transported_cell,pseparator,plocseplink,ploclink,lo,hi,x,y,z,xm,ym,zm,phi,gradxphi_mat,gradyphi_mat,gradzphi_mat,ph_lo,ph_hi,flux,volume)

          ! compute final x-fluxes
	   do icomp=1,ncomp	
              flxx(i,j,k,icomp) = flux(icomp)
	   enddo		 
	   volx(i,j,k) = volume  

	  ! Compute y face flux
          transported_cell(:,1)=original_cell(:,2);transported_cell(:,2)=original_cell(:,1)
          transported_cell(:,3)=original_cell(:,4);transported_cell(:,4)=original_cell(:,3)
          transported_cell(:,5)=projected_cell(:,2);transported_cell(:,6)=projected_cell(:,1)
          transported_cell(:,7)=projected_cell(:,4);transported_cell(:,8)=projected_cell(:,3)
	  transported_cell(:,9)=0.25d0*[sum(transported_cell(1,1:4)),sum(transported_cell(2,1:4)),sum(transported_cell(3,1:4))]
          transported_cell(:,9)=project(transported_cell(:,9),-dt,i,j,k,umac,u_lo,u_hi,vmac,v_lo,v_hi,wmac,w_lo,w_hi,lo,hi,x,y,z,xm,ym,zm)

          call construct(flux_hexahedron(1),transported_cell)
          call adjustCapToMatchVolume(flux_hexahedron(1),dt*vmac(i,j,k)*dx(1)*dx(3))
          call getNormMoments(flux_hexahedron(1),plocseplink(i,j,k),f_moments)
          call compute_fluxes(i,j,k,transported_cell,pseparator,plocseplink,ploclink,lo,hi,x,y,z,xm,ym,zm,phi,gradxphi_mat,gradyphi_mat,gradzphi_mat,ph_lo,ph_hi,flux,volume)

          ! compute final y-fluxes
           do icomp=1,ncomp
              flxy(i,j,k,icomp) = flux(icomp)
           enddo
	   voly(i,j,k) = volume

	  ! Compute z face flux
          transported_cell(:,1)=original_cell(:,5);transported_cell(:,2)=original_cell(:,1)
          transported_cell(:,3)=original_cell(:,2);transported_cell(:,4)=original_cell(:,6)
          transported_cell(:,5)=projected_cell(:,5);transported_cell(:,6)=projected_cell(:,1)
          transported_cell(:,7)=projected_cell(:,2);transported_cell(:,8)=projected_cell(:,6)
	  transported_cell(:,9)=0.25d0*[sum(transported_cell(1,1:4)),sum(transported_cell(2,1:4)),sum(transported_cell(3,1:4))]
          transported_cell(:,9)=project(transported_cell(:,9),-dt,i,j,k,umac,u_lo,u_hi,vmac,v_lo,v_hi,wmac,w_lo,w_hi,lo,hi,x,y,z,xm,ym,zm)

          call construct(flux_hexahedron(1),transported_cell)
          call adjustCapToMatchVolume(flux_hexahedron(1),dt*wmac(i,j,k)*dx(1)*dx(2))
          call getNormMoments(flux_hexahedron(1),plocseplink(i,j,k),f_moments)
          call compute_fluxes(i,j,k,transported_cell,pseparator,plocseplink,ploclink,lo,hi,x,y,z,xm,ym,zm,phi,gradxphi_mat,gradyphi_mat,gradzphi_mat,ph_lo,ph_hi,flux,volume)

          ! compute final z-fluxes
           do icomp=1,ncomp
              flxz(i,j,k,icomp) = flux(icomp)
           enddo
  	   volz(i,j,k) = volume

          end do
       end do
    end do


   ! Last x face

    do  k = lo(3), hi(3)
        do  j = lo(2), hi(2)
          do i = hi(1)+1, hi(1)+1
	 
	  ! Cell nodal points
	  original_cell(:,1)=[x(i),y(j),    z(k)  ] ;original_cell(:,2)=[x(i+1),y(j),    z(k)];
	  original_cell(:,3)=[x(i+1),y(j),z(k+1) ]  ;original_cell(:,4)=[x(i),y(j),    z(k+1)];
	  original_cell(:,5)=[x(i),y(j+1),  z(k)  ] ;original_cell(:,6)=[x(i+1),y(j+1),  z(k)];
	  original_cell(:,7)=[x(i+1),y(j+1), z(k+1)];original_cell(:,8)=[x(i),y(j+1),  z(k+1)];
	 
    	  ! Project all the cell nodal points just once
          projected_cell(:,1) =project(original_cell(:,1),-dt,i,j,k,umac,u_lo,u_hi,vmac,v_lo,v_hi,wmac,w_lo,w_hi,lo,hi,x,y,z,xm,ym,zm)
          projected_cell(:,4) =project(original_cell(:,4),-dt,i,j,k,umac,u_lo,u_hi,vmac,v_lo,v_hi,wmac,w_lo,w_hi,lo,hi,x,y,z,xm,ym,zm)
          projected_cell(:,8) =project(original_cell(:,8),-dt,i,j,k,umac,u_lo,u_hi,vmac,v_lo,v_hi,wmac,w_lo,w_hi,lo,hi,x,y,z,xm,ym,zm)
          projected_cell(:,5) =project(original_cell(:,5),-dt,i,j,k,umac,u_lo,u_hi,vmac,v_lo,v_hi,wmac,w_lo,w_hi,lo,hi,x,y,z,xm,ym,zm)
         
	  ! Compute x face flux 
	  transported_cell(:,1)=original_cell(:,4);transported_cell(:,2)=original_cell(:,1)
	  transported_cell(:,3)=original_cell(:,5);transported_cell(:,4)=original_cell(:,8)
	  transported_cell(:,5)=projected_cell(:,4);transported_cell(:,6)=projected_cell(:,1)
	  transported_cell(:,7)=projected_cell(:,5);transported_cell(:,8)=projected_cell(:,8)
	  transported_cell(:,9)=0.25d0*[sum(transported_cell(1,1:4)),sum(transported_cell(2,1:4)),sum(transported_cell(3,1:4))]
          transported_cell(:,9)=project(transported_cell(:,9),-dt,i,j,k,umac,u_lo,u_hi,vmac,v_lo,v_hi,wmac,w_lo,w_hi,lo,hi,x,y,z,xm,ym,zm)

          call construct(flux_hexahedron(1),transported_cell)
          call adjustCapToMatchVolume(flux_hexahedron(1),dt*umac(i,j,k)*dx(2)*dx(3))
          call getNormMoments(flux_hexahedron(1),plocseplink(i,j,k),f_moments)
          call compute_fluxes(i,j,k,transported_cell,pseparator,plocseplink,ploclink,lo,hi,x,y,z,xm,ym,zm,phi,gradxphi_mat,gradyphi_mat,gradzphi_mat,ph_lo,ph_hi,flux,volume)

          ! compute final x-fluxes
	   do icomp=1,ncomp	
              flxx(i,j,k,icomp) = flux(icomp)
	   enddo		
	   volx(i,j,k) = volume

          end do
       end do
    end do

   ! Last y face

   do  k = lo(3), hi(3)
       do  j = hi(2)+1, hi(2)+1
          do i = lo(1), hi(1)
	 
	  ! Cell nodal points
	  original_cell(:,1)=[x(i),y(j),    z(k)  ] ;original_cell(:,2)=[x(i+1),y(j),    z(k)];
	  original_cell(:,3)=[x(i+1),y(j),z(k+1) ]  ;original_cell(:,4)=[x(i),y(j),    z(k+1)];
	  original_cell(:,5)=[x(i),y(j+1),  z(k)  ] ;original_cell(:,6)=[x(i+1),y(j+1),  z(k)];
	  original_cell(:,7)=[x(i+1),y(j+1), z(k+1)];original_cell(:,8)=[x(i),y(j+1),  z(k+1)];
	 
    	  ! Project all the cell nodal points just once
          projected_cell(:,1) =project(original_cell(:,1),-dt,i,j,k,umac,u_lo,u_hi,vmac,v_lo,v_hi,wmac,w_lo,w_hi,lo,hi,x,y,z,xm,ym,zm)
          projected_cell(:,2) =project(original_cell(:,2),-dt,i,j,k,umac,u_lo,u_hi,vmac,v_lo,v_hi,wmac,w_lo,w_hi,lo,hi,x,y,z,xm,ym,zm)
          projected_cell(:,3) =project(original_cell(:,3),-dt,i,j,k,umac,u_lo,u_hi,vmac,v_lo,v_hi,wmac,w_lo,w_hi,lo,hi,x,y,z,xm,ym,zm)
          projected_cell(:,4) =project(original_cell(:,4),-dt,i,j,k,umac,u_lo,u_hi,vmac,v_lo,v_hi,wmac,w_lo,w_hi,lo,hi,x,y,z,xm,ym,zm)

	  ! Compute y face flux
          transported_cell(:,1)=original_cell(:,2);transported_cell(:,2)=original_cell(:,1)
          transported_cell(:,3)=original_cell(:,4);transported_cell(:,4)=original_cell(:,3)
          transported_cell(:,5)=projected_cell(:,2);transported_cell(:,6)=projected_cell(:,1)
          transported_cell(:,7)=projected_cell(:,4);transported_cell(:,8)=projected_cell(:,3)
	  transported_cell(:,9)=0.25d0*[sum(transported_cell(1,1:4)),sum(transported_cell(2,1:4)),sum(transported_cell(3,1:4))]
          transported_cell(:,9)=project(transported_cell(:,9),-dt,i,j,k,umac,u_lo,u_hi,vmac,v_lo,v_hi,wmac,w_lo,w_hi,lo,hi,x,y,z,xm,ym,zm)

          call construct(flux_hexahedron(1),transported_cell)
          call adjustCapToMatchVolume(flux_hexahedron(1),dt*vmac(i,j,k)*dx(1)*dx(3))
          call getNormMoments(flux_hexahedron(1),plocseplink(i,j,k),f_moments)
          call compute_fluxes(i,j,k,transported_cell,pseparator,plocseplink,ploclink,lo,hi,x,y,z,xm,ym,zm,phi,gradxphi_mat,gradyphi_mat,gradzphi_mat,ph_lo,ph_hi,flux,volume)

          ! compute final y-fluxes
           do icomp=1,ncomp
              flxy(i,j,k,icomp) = flux(icomp)
           enddo
	   voly(i,j,k) = volume

          end do
       end do
    end do

   ! Last z face

    do  k = hi(3)+1, hi(3)+1
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
	 
	  ! Cell nodal points
	  original_cell(:,1)=[x(i),y(j),    z(k)  ] ;original_cell(:,2)=[x(i+1),y(j),    z(k)];
	  original_cell(:,3)=[x(i+1),y(j),z(k+1) ]  ;original_cell(:,4)=[x(i),y(j),    z(k+1)];
	  original_cell(:,5)=[x(i),y(j+1),  z(k)  ] ;original_cell(:,6)=[x(i+1),y(j+1),  z(k)];
	  original_cell(:,7)=[x(i+1),y(j+1), z(k+1)];original_cell(:,8)=[x(i),y(j+1),  z(k+1)];
	 
    	  ! Project all the cell nodal points just once
          projected_cell(:,1) =project(original_cell(:,1),-dt,i,j,k,umac,u_lo,u_hi,vmac,v_lo,v_hi,wmac,w_lo,w_hi,lo,hi,x,y,z,xm,ym,zm)
          projected_cell(:,5) =project(original_cell(:,5),-dt,i,j,k,umac,u_lo,u_hi,vmac,v_lo,v_hi,wmac,w_lo,w_hi,lo,hi,x,y,z,xm,ym,zm)
          projected_cell(:,6) =project(original_cell(:,6),-dt,i,j,k,umac,u_lo,u_hi,vmac,v_lo,v_hi,wmac,w_lo,w_hi,lo,hi,x,y,z,xm,ym,zm)
          projected_cell(:,2) =project(original_cell(:,2),-dt,i,j,k,umac,u_lo,u_hi,vmac,v_lo,v_hi,wmac,w_lo,w_hi,lo,hi,x,y,z,xm,ym,zm)
         
	  ! Compute z face flux
          transported_cell(:,1)=original_cell(:,5);transported_cell(:,2)=original_cell(:,1)
          transported_cell(:,3)=original_cell(:,2);transported_cell(:,4)=original_cell(:,6)
          transported_cell(:,5)=projected_cell(:,5);transported_cell(:,6)=projected_cell(:,1)
          transported_cell(:,7)=projected_cell(:,2);transported_cell(:,8)=projected_cell(:,6)
	  transported_cell(:,9)=0.25d0*[sum(transported_cell(1,1:4)),sum(transported_cell(2,1:4)),sum(transported_cell(3,1:4))]
          transported_cell(:,9)=project(transported_cell(:,9),-dt,i,j,k,umac,u_lo,u_hi,vmac,v_lo,v_hi,wmac,w_lo,w_hi,lo,hi,x,y,z,xm,ym,zm)

          call construct(flux_hexahedron(1),transported_cell)
          call adjustCapToMatchVolume(flux_hexahedron(1),dt*wmac(i,j,k)*dx(1)*dx(2))
          call getNormMoments(flux_hexahedron(1),plocseplink(i,j,k),f_moments)
          call compute_fluxes(i,j,k,transported_cell,pseparator,plocseplink,ploclink,lo,hi,x,y,z,xm,ym,zm,phi,gradxphi_mat,gradyphi_mat,gradzphi_mat,ph_lo,ph_hi,flux,volume)

          ! compute final z-fluxes
           do icomp=1,ncomp
              flxz(i,j,k,icomp) = flux(icomp)
           enddo
	   volz(i,j,k) = volume

          end do
       end do
    end do

  end subroutine compute_flux_3d

subroutine setUpMeshConnectivity(lower_ind, upper_ind, ph_lo, ph_hi, plocseplink, ploclink)
  use irl_fortran_interface
  implicit none
  integer, intent(in) :: ph_lo(3), ph_hi(3)
  integer, dimension(1:3), intent(in) :: lower_ind
  integer, dimension(1:3), intent(in) :: upper_ind
  type(LocSepLink_type), dimension(ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3)),intent(inout) :: plocseplink
  type(LocLink_type), dimension(ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3)),intent(inout) :: ploclink

  integer :: i,j,k
  integer :: unique_id

  ! Connect all the mesh
  unique_id = 0
 do k = lower_ind(3), upper_ind(3)
    do j = lower_ind(2), upper_ind(2)
       do i = lower_ind(1), upper_ind(1)

      ! Set unique ID
      call setId(plocseplink(i,j,k),unique_id)
      call setId(ploclink(i,j,k),unique_id)

      if(i-1 .ge. lower_ind(1)) then
	call setEdgeConnectivity(plocseplink(i,j,k),0,plocseplink(i-1,j,k))
        call setEdgeConnectivity(ploclink(i,j,k),0,ploclink(i-1,j,k))
      end if

      if(i+1 .le. upper_ind(1)) then
        call setEdgeConnectivity(plocseplink(i,j,k),1,plocseplink(i+1,j,k))
        call setEdgeConnectivity(ploclink(i,j,k),1,ploclink(i+1,j,k))
      end if

      if(j-1 .ge. lower_ind(2)) then
	call setEdgeConnectivity(plocseplink(i,j,k),2,plocseplink(i,j-1,k))
        call setEdgeConnectivity(ploclink(i,j,k),2,ploclink(i,j-1,k))
      end if

      if(j+1 .le. upper_ind(2)) then
	call setEdgeConnectivity(plocseplink(i,j,k),3,plocseplink(i,j+1,k))
        call setEdgeConnectivity(ploclink(i,j,k),3,ploclink(i,j+1,k))
      end if

      if(k-1 .ge. lower_ind(3)) then
	call setEdgeConnectivity(plocseplink(i,j,k),4,plocseplink(i,j,k-1))
        call setEdgeConnectivity(ploclink(i,j,k),4,ploclink(i,j,k-1))
      end if

      if(k+1 .le. upper_ind(3)) then
        call setEdgeConnectivity(plocseplink(i,j,k),5,plocseplink(i,j,k+1))
        call setEdgeConnectivity(ploclink(i,j,k),5,ploclink(i,j,k+1))
      end if

      unique_id = unique_id + 1

      end do
    end do
  end do

 end subroutine setUpMeshConnectivity

subroutine compute_fluxes(i,j,k,transported_cell,pseparator,plocseplink,ploclink,lo,hi,x,y,z,xm,ym,zm,phi,gradxphi_mat,gradyphi_mat,gradzphi_mat,ph_lo,ph_hi,flux,volume)

 use amr_data_module
 use irl_fortran_interface
 implicit none

 integer, intent(in) :: i,j,k,ph_lo(3),ph_hi(3),lo(3),hi(3)
 double precision, dimension(3,8), intent(in) :: transported_cell 
 type(LocSepLink_type), intent(inout) :: plocseplink(ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3))
 type(LocLink_type), intent(inout) :: ploclink(ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3))
 type(PlanarSep_type), intent(in) :: pseparator(ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3))
 double precision, intent(in) :: x(lo(1)-1:hi(1)+2)
 double precision, intent(in) :: y(lo(2)-1:hi(2)+2)
 double precision, intent(in) :: z(lo(3)-1:hi(3)+2)
 double precision, intent(in) :: xm(lo(1)-1:hi(1)+1)
 double precision, intent(in) :: ym(lo(2)-1:hi(2)+1)
 double precision, intent(in) :: zm(lo(3)-1:hi(3)+1)
 double precision, intent(in   ) :: phi (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3),ncomp)
 double precision, intent(in   ) :: gradxphi_mat (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3),ncomp)
 double precision, intent(in   ) :: gradyphi_mat (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3),ncomp)
 double precision, intent(in   ) :: gradzphi_mat (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3),ncomp)
 integer,  dimension(3) :: lower_ind, upper_ind
 double precision :: phi_reconstruct(ncomp)
 integer :: id
 integer :: ii,jj,kk,nlocalizers,icomp, n
 double precision, intent(out) :: flux(11)
 double precision, intent(out) :: volume
 double precision :: dx, dy, dz
 double precision :: centroid1(3),centroid2(3),Gbary(3),Lbary(3)
 integer  :: list_size,uniq_id
 integer,  dimension(3) :: ind
 type(SepVM_type) :: my_SepVM
   
  dx=x(lo(1))-x(lo(1)-1)
  dy=y(lo(2))-y(lo(2)-1)
  dz=z(lo(3))-z(lo(3)-1)
 call setMinimumVolToTrack(dx*dy*dz*1.0d-15)

  flux=0.0d0
  volume=0.0d0

! Get number of tags (elements) in the f_moments
  list_size = getSize(f_moments)
  ! Loop through tags in the list
  do n = 0,list_size-1
     ! Get unique id of current cell
     uniq_id = getTagForIndex(f_moments,n)
     ! Convert unique id to indices ii,jj,kk
     ind = getIndicesFromLexTag(uniq_id)
     ii = ind(1); jj = ind(2); kk = ind(3)
   
     ! Get barycenter and volume of tets in current cell
     call getSepVMAtIndex(f_moments,n,my_SepVM)
     centroid1 = getCentroid(my_SepVM, 0)
     centroid2 = getCentroid(my_SepVM, 1)

     ! Bound barycenters by the cell dimensions
     !my_Lbary(1) = max(x(ii),min(x(ii+1),my_Lbary(1))); my_Gbary(1) = max(x(ii),min(x(ii+1),my_Gbary(1)))
     !my_Lbary(2) = max(y(jj),min(y(jj+1),my_Lbary(2))); my_Gbary(2) = max(y(jj),min(y(jj+1),my_Gbary(2)))
     !my_Lbary(3) = max(z(kk),min(z(kk+1),my_Lbary(3))); my_Gbary(3) = max(z(kk),min(z(kk+1),my_Gbary(3)))

     call compute_barycenters(ii,jj,kk,x,y,z,lo,hi,pseparator,ph_lo,ph_hi,Gbary,Lbary)	

     do icomp=1,5
	phi_reconstruct(icomp)=phi(ii,jj,kk,icomp)+gradxphi_mat(ii,jj,kk,icomp)*(centroid2(1)-Gbary(1))+&
                                                   gradyphi_mat(ii,jj,kk,icomp)*(centroid2(2)-Gbary(2))+&
                                                   gradzphi_mat(ii,jj,kk,icomp)*(centroid2(3)-Gbary(3))
     enddo
     do icomp=6,10
	phi_reconstruct(icomp)=phi(ii,jj,kk,icomp)+gradxphi_mat(ii,jj,kk,icomp)*(centroid1(1)-Lbary(1))+&
                                                   gradyphi_mat(ii,jj,kk,icomp)*(centroid1(2)-Lbary(2))+&
                                                   gradzphi_mat(ii,jj,kk,icomp)*(centroid1(3)-Lbary(3))
     enddo

     ! Quantities in gas phase
      do icomp=1,5
         flux(icomp) = flux(icomp) + getVolume(my_SepVM, 1)*phi_reconstruct(icomp)
      end do
      do icomp=6,10
         flux(icomp) = flux(icomp) + getVolume(my_SepVM, 0)*phi_reconstruct(icomp)
      end do
      flux(ncomp) = flux(ncomp) + getVolume(my_SepVM, 0)
      volume = volume + getVolume(my_SepVM, 0) + getVolume(my_SepVM, 1)

   enddo

contains

subroutine compute_barycenters(ii,jj,kk,x,y,z,lo,hi,pseparator,ph_lo,ph_hi,Gbary,Lbary)

 implicit none

 integer, intent(in) :: ii,jj,kk,ph_lo(3),ph_hi(3),lo(3),hi(3)
 double precision, intent(in) :: x(lo(1)-1:hi(1)+2)
 double precision, intent(in) :: y(lo(2)-1:hi(2)+2)
 double precision, intent(in) :: z(lo(3)-1:hi(3)+2)
 double precision, intent(out) :: Gbary(3), Lbary(3)
 type(RectCub_type) :: cell
 type(SepVM_type) :: separated_volume_moments
 type(PlanarSep_type), intent(in) :: pseparator(ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3))
 double precision :: plane(4)

 call new(cell)
 call new(separated_volume_moments)

 call construct_2pt(cell,(/x(ii),y(jj),z(kk)/), (/x(ii+1),y(jj+1),z(kk+1)/) )
 call getNormMoments(cell, pseparator(ii,jj,kk), separated_volume_moments)
 Lbary = getCentroidPtr(separated_volume_moments,0)
 Gbary = getCentroidPtr(separated_volume_moments,1)

 !plane = getPlane(pseparator(ii,jj,kk),0)
 
 !if(abs(plane(1)).gt.0.0d0.and.abs(plane(2)).gt.0.0d0.and.abs(plane(3)).gt.0.0d0)then
	!print*, x(ii),y(jj),z(kk)
	!print*, x(ii+1),y(jj+1),z(kk+1)
	!print*, Lbary, Gbary
	!stop
 !end if
 

end subroutine compute_barycenters

! =================================== !
  ! Lexicographic indexing of our cells !
  ! =================================== !
  pure function getIndicesFromLexTag(a_tag) result(ind)
    implicit none
    integer(IRL_UnsignedIndex_t), intent(in) :: a_tag
    integer, dimension(3) :: ind
    integer :: nxo_,nyo_,nzo_,imino_,jmino_,kmino_
    nxo_ = hi(1)-lo(1)+3
    nyo_ = hi(2)-lo(2)+3
    nzo_ = hi(3)-lo(3)+3
    imino_ = lo(1)-1
    jmino_ = lo(2)-1
    kmino_ = lo(3)-1
    ind(3) = a_tag/(nxo_*nyo_)
    ind(2) = (a_tag - ind(3)*nxo_*nyo_)/nxo_
    ind(1) = a_tag - nxo_*ind(2) - nxo_*nyo_*ind(3)
    ind = ind + (/imino_,jmino_,kmino_/)
    return
  end function getIndicesFromLexTag

end subroutine compute_fluxes
 
end module compute_flux_module
