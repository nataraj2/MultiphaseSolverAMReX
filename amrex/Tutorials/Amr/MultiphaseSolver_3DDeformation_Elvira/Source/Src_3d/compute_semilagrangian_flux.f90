module compute_semilagrangian_flux

use compgeom_lookup
use amr_data_module

contains

  ! ==================== !
  ! Get indices of point !
  ! ==================== !
  function get_indices(pt,ind_in,lo,hi,x,y,z) result(ind)
    implicit none
    double precision, dimension(3), intent(in) :: pt
    integer,  dimension(3), intent(in) :: ind_in
    integer, dimension(3) :: ind
    integer, intent(in) :: lo(3),hi(3)
    double precision, intent(in) :: x(lo(1)-1:hi(1)+2)
    double precision, intent(in) :: y(lo(2)-1:hi(2)+2)
    double precision, intent(in) :: z(lo(3)-1:hi(3)+2)

    ! X direction
    ind(1)=ind_in(1)
    do while (pt(1).gt.x(ind(1)+1)); ind(1)=ind(1)+1; end do
    do while (pt(1).lt.x(ind(1)  )); ind(1)=ind(1)-1; end do
    ! Y direction
    ind(2)=ind_in(2)
    do while (pt(2).gt.y(ind(2)+1)); ind(2)=ind(2)+1; end do
    do while (pt(2).lt.y(ind(2)  )); ind(2)=ind(2)-1; end do
    ! Z direction
    ind(3)=ind_in(3)
    do while (pt(3).gt.z(ind(3)+1)); ind(3)=ind(3)+1; end do
    do while (pt(3).lt.z(ind(3)  )); ind(3)=ind(3)-1; end do
    return
   end function get_indices

  ! ====================== !
  ! Computes volume of tet !
  ! ====================== !
  function tet_vol(vert) result(tetvol)
    implicit none
    double precision :: tetvol
    double precision, dimension(3,4), intent(in) :: vert
    double precision, dimension(3) :: a,b,c
    a=vert(:,1)-vert(:,4); b=vert(:,2)-vert(:,4); c=vert(:,3)-vert(:,4)
    tetvol=-(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0d0
    return
  end function tet_vol

  ! ===================== !
  ! Calculate sign of tet !
  ! ===================== !
  function tet_sign(vert) result(s)
    implicit none
    double precision :: s
    double precision, dimension(3,4), intent(in) :: vert
    double precision, dimension(3) :: a,b,c
    a=vert(:,1)-vert(:,4); b=vert(:,2)-vert(:,4); c=vert(:,3)-vert(:,4)
    s=sign(1.0d0,-(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0d0)
    return
  end function tet_sign

  ! ========================================= !
  ! Move point p1 to p2 according to velocity !
  ! ========================================= !
  function project(p1,mydt,myi,myj,myk,umac,u_lo,u_hi,&
				       vmac,v_lo,v_hi,&
				       wmac,w_lo,w_hi,&
				       lo,hi,x,y,z,xm,ym,zm) result(p2)
    implicit none
    double precision, dimension(3) :: p2
    double precision, dimension(3), intent(in) :: p1
    double precision,               intent(in) :: mydt
    integer,                intent(in) :: myi,myj,myk
    double precision, dimension(3) :: v1,v2,v3,v4
    integer, intent(in) ::  u_lo(3),  u_hi(3), v_lo(3), v_hi(3), w_lo(3), w_hi(3)
    double precision, intent(in   ) :: umac( u_lo(1): u_hi(1), u_lo(2): u_hi(2), u_lo(3): u_hi(3))
    double precision, intent(in   ) :: vmac( v_lo(1): v_hi(1), v_lo(2): v_hi(2), v_lo(3): v_hi(3))
    double precision, intent(in   ) :: wmac( w_lo(1): w_hi(1), w_lo(2): w_hi(2), w_lo(3): w_hi(3))
    integer, intent(in) :: lo(3),hi(3)
    double precision, intent(in) :: x(lo(1)-1:hi(1)+2)
    double precision, intent(in) :: y(lo(2)-1:hi(2)+2)
    double precision, intent(in) :: z(lo(3)-1:hi(3)+2)
    double precision, intent(in) :: xm(lo(1)-1:hi(1)+1)
    double precision, intent(in) :: ym(lo(2)-1:hi(2)+1)
    double precision, intent(in) :: zm(lo(3)-1:hi(3)+1)

    v1=get_velocity(p1              ,myi,myj,myk,umac,u_lo,u_hi,&
						 vmac,v_lo,v_hi,&
						 wmac,w_lo,w_hi,&
                                                 lo,hi,x,y,z,xm,ym,zm)
    v2=get_velocity(p1+0.5d0*mydt*v1,myi,myj,myk,umac,u_lo,u_hi,&
						 vmac,v_lo,v_hi,&
						 wmac,w_lo,w_hi,&
                                                 lo,hi,x,y,z,xm,ym,zm)
    v3=get_velocity(p1+0.5d0*mydt*v2,myi,myj,myk,umac,u_lo,u_hi,&
						 vmac,v_lo,v_hi,&
						 wmac,w_lo,w_hi,&
                                                 lo,hi,x,y,z,xm,ym,zm)
    v4=get_velocity(p1+      mydt*v3,myi,myj,myk,umac,u_lo,u_hi,&
						 vmac,v_lo,v_hi,&
						 wmac,w_lo,w_hi,&
                                                 lo,hi,x,y,z,xm,ym,zm)
    p2=p1+mydt/6.0d0*(v1+2.0d0*v2+2.0d0*v3+v4)

  end function project


  ! ================================================ !
  ! Interpolate velocity to pos near cell (i0,j0,k0) !
  ! ================================================ !
  function get_velocity(pos,i0,j0,k0,umac,u_lo,u_hi,&
				     vmac,v_lo,v_hi,&
				     wmac,w_lo,w_hi,&
                                     lo,hi,x,y,z,xm,ym,zm) result(vel)
    implicit none
    double precision, dimension(3) :: vel
    double precision, dimension(3), intent(in) :: pos
    integer, intent(in) :: i0,j0,k0
    integer :: i,j,k
    double precision :: wx1,wy1,wz1
    double precision :: wx2,wy2,wz2
    double precision, dimension(2,2,2) :: ww
    double precision :: x1,x2,y1,y2,z1,z2
    double precision :: xm1,xm2,ym1,ym2,zm1,zm2
    integer :: imino_,imaxo_,jmino_,jmaxo_,kmino_,kmaxo_
    integer, intent(in) ::  u_lo(3),  u_hi(3), v_lo(3), v_hi(3), w_lo(3), w_hi(3)
    double precision, intent(in   ) :: umac( u_lo(1): u_hi(1), u_lo(2): u_hi(2), u_lo(3): u_hi(3))
    double precision, intent(in   ) :: vmac( v_lo(1): v_hi(1), v_lo(2): v_hi(2), v_lo(3): v_hi(3))
    double precision, intent(in   ) :: wmac( w_lo(1): w_hi(1), w_lo(2): w_hi(2), w_lo(3): w_hi(3))

    integer, intent(in) :: lo(3),hi(3)
    double precision, intent(in) :: x(lo(1)-1:hi(1)+2)
    double precision, intent(in) :: y(lo(2)-1:hi(2)+2)
    double precision, intent(in) :: z(lo(3)-1:hi(3)+2)
    double precision, intent(in) :: xm(lo(1)-1:hi(1)+1)
    double precision, intent(in) :: ym(lo(2)-1:hi(2)+1)
    double precision, intent(in) :: zm(lo(3)-1:hi(3)+1)

    vel=0.0d0

    imino_=lo(1)-1;imaxo_=hi(1)+2
    jmino_=lo(2)-1;jmaxo_=hi(2)+2
    kmino_=lo(3)-1;kmaxo_=hi(3)+2

! Interpolate Uface velocity ------------------------------
    ! Find right i index
    i=max(min(imaxo_-1,i0),imino_)
    do while (pos(1)-x (i  ).lt.0.0d0.and.i  .gt.imino_); i=i-1; end do
    do while (pos(1)-x (i+1).ge.0.0d0.and.i+1.lt.imaxo_); i=i+1; end do
    ! Find right j index
    j=max(min(jmaxo_-1,j0),jmino_)
    do while (pos(2)-ym(j  ).lt.0.0d0.and.j  .gt.jmino_); j=j-1; end do
    do while (pos(2)-ym(j+1).ge.0.0d0.and.j+1.lt.jmaxo_); j=j+1; end do
    ! Find right k index
    k=max(min(kmaxo_-1,k0),kmino_)
    do while (pos(3)-zm(k  ).lt.0.0d0.and.k  .gt.kmino_); k=k-1; end do
    do while (pos(3)-zm(k+1).ge.0.0d0.and.k+1.lt.kmaxo_); k=k+1; end do
    ! Prepare tri-linear interpolation coefficients
    wx1=(pos(1)-x (i))/(x (i+1)-x (i)); wx2=1.0d0-wx1
    wy1=(pos(2)-ym(j))/(ym(j+1)-ym(j)); wy2=1.0d0-wy1
    wz1=(pos(3)-zm(k))/(zm(k+1)-zm(k)); wz2=1.0d0-wz1
    ! Recast in array form
    ww(1,1,1)=wx2*wy2*wz2
    ww(2,1,1)=wx1*wy2*wz2
    ww(1,2,1)=wx2*wy1*wz2
    ww(2,2,1)=wx1*wy1*wz2
    ww(1,1,2)=wx2*wy2*wz1
    ww(2,1,2)=wx1*wy2*wz1
    ww(1,2,2)=wx2*wy1*wz1
    ww(2,2,2)=wx1*wy1*wz1
 
   ! Tri-linear interpolation of Uface
    vel(1)=sum(ww*umac(i:i+1,j:j+1,k:k+1))

! Interpolate Vface velocity ------------------------------
    ! Find right i index
    i=max(min(imaxo_-1,i0),imino_)
    do while (pos(1)-xm(i  ).lt.0.0d0.and.i  .gt.imino_); i=i-1; end do
    do while (pos(1)-xm(i+1).ge.0.0d0.and.i+1.lt.imaxo_); i=i+1; end do
    ! Find right j index
    j=max(min(jmaxo_-1,j0),jmino_)
    do while (pos(2)-y (j  ).lt.0.0d0.and.j  .gt.jmino_); j=j-1; end do
    do while (pos(2)-y (j+1).ge.0.0d0.and.j+1.lt.jmaxo_); j=j+1; end do
    ! Find right k index
    k=max(min(kmaxo_-1,k0),kmino_)
    do while (pos(3)-zm(k  ).lt.0.0d0.and.k  .gt.kmino_); k=k-1; end do
    do while (pos(3)-zm(k+1).ge.0.0d0.and.k+1.lt.kmaxo_); k=k+1; end do
    ! Prepare tri-linear interpolation coefficients
    wx1=(pos(1)-xm(i))/(xm(i+1)-xm(i)); wx2=1.0d0-wx1
    wy1=(pos(2)-y (j))/(y (j+1)-y (j)); wy2=1.0d0-wy1
    wz1=(pos(3)-zm(k))/(zm(k+1)-zm(k)); wz2=1.0d0-wz1
    ! Recast in array form
    ww(1,1,1)=wx2*wy2*wz2
    ww(2,1,1)=wx1*wy2*wz2
    ww(1,2,1)=wx2*wy1*wz2
    ww(2,2,1)=wx1*wy1*wz2
    ww(1,1,2)=wx2*wy2*wz1
    ww(2,1,2)=wx1*wy2*wz1
    ww(1,2,2)=wx2*wy1*wz1
    ww(2,2,2)=wx1*wy1*wz1
    ! Tri-linear interpolation of Vface
    vel(2)=sum(ww*vmac(i:i+1,j:j+1,k:k+1))

    ! Interpolate Wface velocity ------------------------------
    ! Find right i index
    i=max(min(imaxo_-1,i0),imino_)
    do while (pos(1)-xm(i  ).lt.0.0d0.and.i  .gt.imino_); i=i-1; end do
    do while (pos(1)-xm(i+1).ge.0.0d0.and.i+1.lt.imaxo_); i=i+1; end do
    ! Find right j index
    j=max(min(jmaxo_-1,j0),jmino_)
    do while (pos(2)-ym(j  ).lt.0.0d0.and.j  .gt.jmino_); j=j-1; end do
    do while (pos(2)-ym(j+1).ge.0.0d0.and.j+1.lt.jmaxo_); j=j+1; end do
    ! Find right k index
    k=max(min(kmaxo_-1,k0),kmino_)
    do while (pos(3)-z (k  ).lt.0.0d0.and.k  .gt.kmino_); k=k-1; end do
    do while (pos(3)-z (k+1).ge.0.0d0.and.k+1.lt.kmaxo_); k=k+1; end do
    ! Prepare tri-linear interpolation coefficients
    wx1=(pos(1)-xm(i))/(xm(i+1)-xm(i)); wx2=1.0d0-wx1
    wy1=(pos(2)-ym(j))/(ym(j+1)-ym(j)); wy2=1.0d0-wy1
    wz1=(pos(3)-z (k))/(z (k+1)-z (k)); wz2=1.0d0-wz1
    ! Recast in array form
    ww(1,1,1)=wx2*wy2*wz2
    ww(2,1,1)=wx1*wy2*wz2
    ww(1,2,1)=wx2*wy1*wz2
    ww(2,2,1)=wx1*wy1*wz2
    ww(1,1,2)=wx2*wy2*wz1
    ww(2,1,2)=wx1*wy2*wz1
    ww(1,2,2)=wx2*wy1*wz1
    ww(2,2,2)=wx1*wy1*wz1
    ! Tri-linear interpolation of Uface
    vel(3)=sum(ww*wmac(i:i+1,j:j+1,k:k+1))


  return
  end function get_velocity


recursive function tet2flux(tet,ind,dx,ph_lo,ph_hi,phi_new_mat,&
				      gradxphi_mat,gradyphi_mat,gradzphi_mat,&
				      lo,hi,x,y,z,xm,ym,zm) result(fluxes)
    implicit none
    double precision, dimension(3,4), intent(in) :: tet
    double precision, intent(in) :: dx(3)
    integer,  dimension(3,4), intent(in) :: ind
    double precision, dimension(13) :: fluxes
    integer :: dir,cut_ind
    integer :: n,nn,case
    integer :: v1,v2
    double precision, dimension(4) :: d
    double precision, dimension(3,8) :: vert
    integer,  dimension(3,8,2) :: vert_ind
    double precision :: mu
    double precision, dimension(3,4) :: newtet
    integer,  dimension(3,4) :: newind

    integer, intent(in) :: ph_lo(3), ph_hi(3)
    double precision, intent(in   ) :: phi_new_mat (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3),ncomp)
    double precision, intent(in   ) :: gradxphi_mat (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3),ncomp)
    double precision, intent(in   ) :: gradyphi_mat (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3),ncomp)
    double precision, intent(in   ) :: gradzphi_mat (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3),ncomp)

    integer, intent(in) :: lo(3),hi(3)
    double precision, intent(in) :: x(lo(1)-1:hi(1)+2)
    double precision, intent(in) :: y(lo(2)-1:hi(2)+2)
    double precision, intent(in) :: z(lo(3)-1:hi(3)+2)
    double precision, intent(in) :: xm(lo(1)-1:hi(1)+1)
    double precision, intent(in) :: ym(lo(2)-1:hi(2)+1)
    double precision, intent(in) :: zm(lo(3)-1:hi(3)+1)

    ! Cut by x planes
    if (maxval(ind(1,:))-minval(ind(1,:)).gt.0) then
       dir=1	
       cut_ind=maxval(ind(1,:))
       d(:)=tet(1,:)-x(cut_ind)
    ! Cut by y planes
    else if (maxval(ind(2,:))-minval(ind(2,:)).gt.0) then
       dir=2
       cut_ind=maxval(ind(2,:))
       d(:)=tet(2,:)-y(cut_ind)
    ! Cut by z planes
    else if (maxval(ind(3,:))-minval(ind(3,:)).gt.0) then
       dir=3
       cut_ind=maxval(ind(3,:))
       d(:)=tet(3,:)-z(cut_ind)
    ! Cut by interface and compute fluxes
    else
       fluxes=0.0d0
       fluxes=tet2flux_plic(tet,ind(1,1),ind(2,1),ind(3,1),ph_lo,ph_hi,phi_new_mat,&
							   gradxphi_mat,gradyphi_mat,gradzphi_mat,&
							   lo,hi,xm,ym,zm)
       return
    end if

    ! Find case of cut
    case=1+int(0.5d0+sign(0.5d0,d(1)))+&
         2*int(0.5d0+sign(0.5d0,d(2)))+&
         4*int(0.5d0+sign(0.5d0,d(3)))+&
         8*int(0.5d0+sign(0.5d0,d(4)))

	!print*, "Case is", case
	!stop

 ! Get vertices and indices of tet
    do n=1,4
       vert    ( : ,n  )=tet(:,n)
       vert_ind( : ,n,1)=ind(:,n)
       vert_ind( : ,n,2)=ind(:,n)
       vert_ind(dir,n,1)=min(vert_ind(dir,n,1),cut_ind-1) ! Enforce boundedness
       vert_ind(dir,n,2)=max(vert_ind(dir,n,1),cut_ind  )
    end do

    ! Create interpolated vertices on cut plane
    do n=1,cut_nvert(case)
       v1=cut_v1(n,case); v2=cut_v2(n,case)
       mu=min(1.0d0,max(0.0d0,-d(v1)/(sign(abs(d(v2)-d(v1))+tiny(1.0d0),d(v2)-d(v1)))))
       vert(:,4+n)=(1.0d0-mu)*vert(:,v1)+mu*vert(:,v2)
       ! Get index for interpolated vertex
       vert_ind(:,4+n,1)=get_indices(vert(:,4+n),vert_ind(:,v1,1),lo,hi,x,y,z)
       ! Enforce boundedness
       vert_ind(:,4+n,1)=max(vert_ind(:,4+n,1),min(vert_ind(:,v1,1),vert_ind(:,v2,1)))
       vert_ind(:,4+n,1)=min(vert_ind(:,4+n,1),max(vert_ind(:,v1,1),vert_ind(:,v2,1)))
       ! Set +/- indices in cut direction
       vert_ind(:,4+n,2)=vert_ind(:,4+n,1)
       vert_ind(dir,4+n,1)=cut_ind-1
       vert_ind(dir,4+n,2)=cut_ind
    end do

    ! Create new tets
    fluxes=0.0d0
    do n=1,cut_ntets(case)
       do nn=1,4
          newtet(:,nn)=vert    (:,cut_vtet(nn,n,case))
          newind(:,nn)=vert_ind(:,cut_vtet(nn,n,case),cut_side(n,case))
       end do
       ! Check for zero-volume tet
       !if (abs(tet_vol(newtet)).lt.volume_epsilon) cycle
       if (abs(tet_vol(newtet)).lt.1.0e-15*(min(dx(1),dx(2),dx(3)))**3) cycle
       ! Cut by next plane
       fluxes=fluxes+tet2flux(newtet,newind,dx,ph_lo,ph_hi,phi_new_mat,&
					       gradxphi_mat,gradyphi_mat,gradzphi_mat,&
					       lo,hi,x,y,z,xm,ym,zm)
    end do

    return
  end function tet2flux

function tet2flux_plic(tet,i,j,k,ph_lo,ph_hi,phi_new_mat,&
			         gradxphi_mat,gradyphi_mat,gradzphi_mat,&
				 lo,hi,xm,ym,zm) result(fluxes)
    implicit none
    double precision, dimension(3,4), intent(in) :: tet
    integer, intent(in) :: i,j,k
    double precision, dimension(13) :: fluxes
    integer :: n,case,v1,v2
    double precision :: mu,my_vol,my_rho,my_rhoE,my_rhoU,my_rhoV,my_rhoW,my_PA
    double precision, dimension(4) :: d
    double precision, dimension(3,8) :: vert
    double precision, dimension(3) :: a,b,c,bary
    double precision :: phi_reconstruct(ncomp)
    integer, intent(in) :: ph_lo(3), ph_hi(3)
    double precision, intent(in), dimension (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3),ncomp) :: phi_new_mat, gradxphi_mat, gradyphi_mat, gradzphi_mat
    integer :: icomp  

    integer, intent(in) :: lo(3),hi(3)
    double precision, intent(in) :: xm(lo(1)-1:hi(1)+1)
    double precision, intent(in) :: ym(lo(2)-1:hi(2)+1)
    double precision, intent(in) :: zm(lo(3)-1:hi(3)+1)


    ! Check if wall
    !if (vol(i,j,k).eq.0.0d0) then
    !   fluxes=0.0d0
    !   return
    !end if

    ! Cut by PLIC
    !d(:)=oldxnorm(i,j,k)*tet(1,:)+oldynorm(i,j,k)*tet(2,:)+oldznorm(i,j,k)*tet(3,:)-olddist(i,j,k)
     d(:)=0.0d0

    ! Find cut case
    case=1+int(0.5d0+sign(0.5d0,d(1)))+&
         2*int(0.5d0+sign(0.5d0,d(2)))+&
         4*int(0.5d0+sign(0.5d0,d(3)))+&
         8*int(0.5d0+sign(0.5d0,d(4)))

    ! Copy vertices
    do n=1,4
       vert(:,n)=tet(:,n)
    end do

    ! Create interpolated vertices on cut plane
    do n=1,cut_nvert(case)
       v1=cut_v1(n,case); v2=cut_v2(n,case)
       mu=min(1.0d0,max(0.0d0,-d(v1)/(sign(abs(d(v2)-d(v1))+tiny(1.0d0),d(v2)-d(v1)))))
       vert(:,4+n)=(1.0d0-mu)*vert(:,v1)+mu*vert(:,v2)
    end do

    ! Zero fluxes
    fluxes=0.0d0

   ! Analyze gas tets
    do n=1,cut_nntet(case)-1
       ! Compute volume
       a=vert(:,cut_vtet(1,n,case))-vert(:,cut_vtet(4,n,case)); b=vert(:,cut_vtet(2,n,case))-vert(:,cut_vtet(4,n,case)); c=vert(:,cut_vtet(3,n,case))-vert(:,cut_vtet(4,n,case))
       my_vol=abs(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0d0
       ! Compute barycenter
       bary=0.25d0*(vert(:,cut_vtet(1,n,case))+vert(:,cut_vtet(2,n,case))+vert(:,cut_vtet(3,n,case))+vert(:,cut_vtet(4,n,case)))-[xm(i),ym(j),zm(k)]
       ! Compute tet variables
	!print*, bary
	
	do icomp=1,ncomp
	   phi_reconstruct(icomp)=phi_new_mat(i,j,k,icomp)+gradxphi_mat(i,j,k,icomp)*bary(1)+&
	   					         gradyphi_mat(i,j,k,icomp)*bary(2)+&
	   						 gradzphi_mat(i,j,k,icomp)*bary(3)
	enddo

      ! Quantities in gas phase
       do icomp=1,ncomp
          fluxes(icomp) = fluxes(icomp) + my_vol*phi_reconstruct(icomp)
       end do
    end do

    return
  end function tet2flux_plic

end module compute_semilagrangian_flux
