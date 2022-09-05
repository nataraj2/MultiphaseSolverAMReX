module advect_module

  use amrex_base_module

  use amr_data_module
  use my_amr_module
  use irl_fortran_interface

  implicit none
  private

  public :: advect

contains

subroutine advect(level,time, lo, hi, &
     &            uin , ui_lo, ui_hi, &
     &            uout, uo_lo, uo_hi, &
     &            vx  , vx_lo, vx_hi, &
     &            vy  , vy_lo, vy_hi, &
     &            vz  , vz_lo, vz_hi, &
     &            flxx, fx_lo, fx_hi, &
     &            flxy, fy_lo, fy_hi, &
     &            flxz, fz_lo, fz_hi, &
     &		  gradxphi_mat,gradyphi_mat,gradzphi_mat, &
     &            plocalizer,pseparator,plocseplink,ploclink,&
     &            x,y,z,xm,ym,zm,dx,prob_lo,dt) !bind(C, name="advect")
  
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use compute_flux_module, only : compute_flux_3d

  implicit none

  integer, intent(in) :: level, lo(3), hi(3)
  double precision, intent(in) :: dx(3), prob_lo(3), dt, time
  integer, intent(in) :: ui_lo(3), ui_hi(3)
  integer, intent(in) :: uo_lo(3), uo_hi(3)
  integer, intent(in) :: vx_lo(3), vx_hi(3)
  integer, intent(in) :: vy_lo(3), vy_hi(3)
  integer, intent(in) :: vz_lo(3), vz_hi(3)
  integer, intent(in) :: fx_lo(3), fx_hi(3)
  integer, intent(in) :: fy_lo(3), fy_hi(3)
  integer, intent(in) :: fz_lo(3), fz_hi(3)
  !double precision, intent(in   ) :: uin (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3),ncomp)
  double precision :: uin (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3),ncomp)
  double precision, intent(inout) :: uout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3),ncomp)
  double precision, intent(in   ) :: vx  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2),vx_lo(3):vx_hi(3))
  double precision, intent(in   ) :: vy  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2),vy_lo(3):vy_hi(3))
  !double precision, intent(in   ) :: vz  (vz_lo(1):vz_hi(1),vz_lo(2):vz_hi(2),vz_lo(3):vz_hi(3))
  double precision :: vz  (vz_lo(1):vz_hi(1),vz_lo(2):vz_hi(2),vz_lo(3):vz_hi(3))
  double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),ncomp)
  double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),ncomp)
  double precision, intent(  out) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),ncomp)
  double precision :: volx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)) 
  double precision :: voly(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3)) 
  double precision :: volz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1) 
  double precision, intent(in   ) :: gradxphi_mat (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3),ncomp)
  double precision, intent(in   ) :: gradyphi_mat (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3),ncomp)
  double precision, intent(in   ) :: gradzphi_mat (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3),ncomp)
  double precision, intent(in) :: x(lo(1)-1:hi(1)+2)
  double precision, intent(in) :: y(lo(2)-1:hi(2)+2)
  double precision, intent(in) :: z(lo(3)-1:hi(3)+2)
  double precision, intent(in) :: xm(lo(1)-1:hi(1)+1)
  double precision, intent(in) :: ym(lo(2)-1:hi(2)+1)
  double precision, intent(in) :: zm(lo(3)-1:hi(3)+1)
  type(PlanarLoc_type), intent(in)        :: plocalizer(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3))
  type(PlanarSep_type), intent(in)        :: pseparator(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3))
  type(LocSepLink_type), intent(inout) :: plocseplink(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3))
  type(LocLink_type), intent(inout) :: ploclink(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3))

  
  integer :: i, j, k, ii, jj, kk
  double precision :: dtdx(3), umax, vmax, wmax
  
  integer :: icomp
  integer :: iloc,jloc,kloc,locarray(3)
  double precision :: xval,yval,zval,rad

  dtdx = dt/dx

  ! We like to allocate these **pointers** here and then pass them to a function
  ! to remove their pointerness for performance, because normally pointers could
  ! be aliasing.  We need to use pointers instead of allocatable arrays because
  ! we like to use BoxLib's bl_allocate to allocate memeory instead of the intrinsic
  ! allocate.  Bl_allocate is much faster than allocate inside OMP.  
  ! Note that one MUST CALL BL_DEALLOCATE.



  
  umax = maxval(abs(vx))
  vmax = maxval(abs(vy))
  wmax = maxval(abs(vz))
  if ( umax*dt .ge. dx(1) .or. &
       vmax*dt .ge. dx(2) .or. &
       wmax*dt .ge. dx(3) ) then
     print *, "umax = ", umax, ", vmax = ", vmax, ", wmax = ", wmax, ", dt = ", dt, ", dx = ", dx
     locarray=maxloc(abs(vx))
     iloc=vx_lo(1)+locarray(1)-1
     jloc=vx_lo(2)+locarray(2)-1
     kloc=vx_lo(3)+locarray(3)-1
     print*,"Location is", prob_lo(1)+(dble(iloc)+0.5)*dx(1),prob_lo(2)+(dble(jloc)+0.5)*dx(2),prob_lo(3)+(dble(kloc)+0.5)*dx(3)
	print*,"Value is",locarray,vx(iloc,jloc,kloc), uin(iloc,jloc,kloc,ncomp)
	print*,uin(iloc-1,jloc,kloc,ncomp),uin(iloc+1,jloc,kloc,ncomp)
	print*,uin(iloc,jloc-1,kloc,ncomp),uin(iloc,jloc+1,kloc,ncomp)
	print*,uin(iloc,jloc,kloc-1,ncomp),uin(iloc,jloc,kloc+1,ncomp)
     call bl_error("CFL violation. Use smaller adv.cfl.")
  end if

  ! call a function to compute flux
  call compute_flux_3d(lo, hi, dt, dx,   &
                       uin, ui_lo, ui_hi,&
                       vx, vx_lo, vx_hi, &
                       vy, vy_lo, vy_hi, &
                       vz, vz_lo, vz_hi, &
                       flxx, fx_lo, fx_hi, &
                       flxy, fy_lo, fy_hi, &
                       flxz, fz_lo, fz_hi, &
		       volx,voly,volz, &
		       gradxphi_mat,gradyphi_mat,gradzphi_mat,&
	               plocalizer,pseparator,plocseplink,ploclink,&
		       x,y,z,xm,ym,zm)

 
  ! Do a conservative update
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)

	   ! First do VOF
           uout(i,j,k,ncomp) = (uin(i,j,k,ncomp)*dx(1)*dx(2)*dx(3) + &
                ( (flxx(i,j,k,ncomp) - flxx(i+1,j,k,ncomp)) &
                + (flxy(i,j,k,ncomp) - flxy(i,j+1,k,ncomp)) &
                + (flxz(i,j,k,ncomp) - flxz(i,j,k+1,ncomp)) ))/(dx(1)*dx(2)*dx(3)+volx(i,j,k)-volx(i+1,j,k)+&
								   		voly(i,j,k)-voly(i,j+1,k)+&
										volz(i,j,k)-volz(i,j,k+1))
		if(uout(i,j,k,ncomp).gt.VOFhi)uout(i,j,k,ncomp)=1.0d0
		if(uout(i,j,k,ncomp).lt.VOFlo)uout(i,j,k,ncomp)=0.0d0
		if(y(j).gt.0.0032d0)uout(i,j,k,ncomp)=0.0d0
		if(x(i).gt.0.0070d0)uout(i,j,k,ncomp)=0.0d0
		if(z(k).gt.0.0115d0)uout(i,j,k,ncomp)=0.0d0
		if(level.lt.amrex_get_finest_level())uout(i,j,k,ncomp)=0.0d0
		
           if(uout(i,j,k,ncomp).lt.VOFlo)then
	   	do icomp=1,5
	           uout(i,j,k,icomp) = uin(i,j,k,icomp)*(1.0d0-uin(i,j,k,ncomp)) + &
        			        ( (flxx(i,j,k,icomp) - flxx(i+1,j,k,icomp)) &
			                + (flxy(i,j,k,icomp) - flxy(i,j+1,k,icomp)) &
			                + (flxz(i,j,k,icomp) - flxz(i,j,k+1,icomp)) )/(dx(1)*dx(2)*dx(3))
    	        enddo
	
	   	do icomp=6,10	
	           uout(i,j,k,icomp)=0.0d0
	        enddo
                uout(i,j,k,ncomp)=0.0d0
           elseif(uout(i,j,k,ncomp).gt.VOFhi)then
	 	 do icomp=1,5
	            uout(i,j,k,icomp)=0.0d0
	   	 enddo
		 do icomp=6,10
           	 	uout(i,j,k,icomp) = uin(i,j,k,icomp)*uin(i,j,k,ncomp) + &
		 	                ( (flxx(i,j,k,icomp) - flxx(i+1,j,k,icomp))  &
		 	                + (flxy(i,j,k,icomp) - flxy(i,j+1,k,icomp))  &
		 	                + (flxz(i,j,k,icomp) - flxz(i,j,k+1,icomp)) )/(dx(1)*dx(2)*dx(3))
		 enddo
                 uout(i,j,k,ncomp)=1.0d0
           else
		do icomp=1,5
	                uout(i,j,k,icomp) = (uin(i,j,k,icomp)*(1.0d0-uin(i,j,k,ncomp)) + &
        			        ( (flxx(i,j,k,icomp) - flxx(i+1,j,k,icomp)) &
			                + (flxy(i,j,k,icomp) - flxy(i,j+1,k,icomp)) &
			                + (flxz(i,j,k,icomp) - flxz(i,j,k+1,icomp)) )/(dx(1)*dx(2)*dx(3)))/(1.0d0-uout(i,j,k,ncomp))
		enddo
		do icomp=6,10
	                uout(i,j,k,icomp) = (uin(i,j,k,icomp)*uin(i,j,k,ncomp) + &
        			        ( (flxx(i,j,k,icomp) - flxx(i+1,j,k,icomp)) &
			                + (flxy(i,j,k,icomp) - flxy(i,j+1,k,icomp)) &
			                + (flxz(i,j,k,icomp) - flxz(i,j,k+1,icomp)) )/(dx(1)*dx(2)*dx(3)))/uout(i,j,k,ncomp)
		enddo
           endif


        enddo
     enddo
  enddo

  

  do icomp=1,ncomp
  ! Scale by face area in order to correctly reflx
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)+1
           flxx(i,j,k,icomp) = flxx(i,j,k,icomp) !* (dt * dx(2)*dx(3))
        enddo
     enddo
  enddo
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)+1 
        do i = lo(1), hi(1)
           flxy(i,j,k,icomp) = flxy(i,j,k,icomp) !* (dt * dx(1)*dx(3))
        enddo
     enddo
  enddo
  do       k = lo(3), hi(3)+1
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           flxz(i,j,k,icomp) = flxz(i,j,k,icomp) !* (dt * dx(1)*dx(2))
        enddo
     enddo
  enddo

  enddo

  do icomp=11,11
  ! Scale by face area in order to correctly reflx
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)+1
           flxx(i,j,k,icomp) = 0.0d0
        enddo
     enddo
  enddo
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)+1 
        do i = lo(1), hi(1)
           flxy(i,j,k,icomp) = 0.0d0
        enddo
     enddo
  enddo
  do       k = lo(3), hi(3)+1
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           flxz(i,j,k,icomp) = 0.0d0
        enddo
     enddo
  enddo

  enddo




end subroutine advect

end module advect_module
