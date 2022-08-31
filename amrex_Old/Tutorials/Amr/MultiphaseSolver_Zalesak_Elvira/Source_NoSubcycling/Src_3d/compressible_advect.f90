module compressible_advect
  use data
  use compressible
  use compgeom_lookup
  implicit none
  
  ! Small reference volume
  real(WP) :: volume_epsilon
  
  ! Advective fluxes
  real(WP), dimension(:,:,:), pointer :: F_VOL,F_VOF,F_PA
  real(WP), dimension(:,:,:), pointer :: F_Grho,F_GrhoE,F_GrhoU,F_GrhoV,F_GrhoW
  real(WP), dimension(:,:,:), pointer :: F_Lrho,F_LrhoE,F_LrhoU,F_LrhoV,F_LrhoW
  
  ! Gradients
  real(WP), dimension(:,:,:,:), allocatable :: gradGrho,gradGrhoE,gradGrhoU,gradGrhoV,gradGrhoW
  real(WP), dimension(:,:,:,:), allocatable :: gradLrho,gradLrhoE,gradLrhoU,gradLrhoV,gradLrhoW
  
  ! Barycenter for each phase
  real(WP), dimension(:,:,:,:), allocatable :: Gbary,Lbary
  
contains
  
  ! ==================== !
  ! Get indices of point !
  ! ==================== !
  function get_indices(pt,ind_in) result(ind)
    implicit none
    real(WP), dimension(3), intent(in) :: pt
    integer,  dimension(3), intent(in) :: ind_in
    integer, dimension(3) :: ind
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
    real(WP) :: tetvol
    real(WP), dimension(3,4), intent(in) :: vert
    real(WP), dimension(3) :: a,b,c
    a=vert(:,1)-vert(:,4); b=vert(:,2)-vert(:,4); c=vert(:,3)-vert(:,4)
    tetvol=-(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
    return
  end function tet_vol
    
  ! ===================== !
  ! Calculate sign of tet !
  ! ===================== !
  function tet_sign(vert) result(s)
    implicit none
    real(WP) :: s
    real(WP), dimension(3,4), intent(in) :: vert
    real(WP), dimension(3) :: a,b,c
    a=vert(:,1)-vert(:,4); b=vert(:,2)-vert(:,4); c=vert(:,3)-vert(:,4)
    s=sign(1.0_WP,-(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP)
    return
  end function tet_sign
    
  ! ========================================= !
  ! Move point p1 to p2 according to velocity !
  ! ========================================= !
  function project(p1,mydt,myi,myj,myk) result(p2)
    implicit none
    real(WP), dimension(3) :: p2
    real(WP), dimension(3), intent(in) :: p1
    real(WP),               intent(in) :: mydt
    integer,                intent(in) :: myi,myj,myk
    real(WP), dimension(3) :: v1,v2,v3,v4
    v1=get_velocity(p1               ,myi,myj,myk)
    v2=get_velocity(p1+0.5_WP*mydt*v1,myi,myj,myk)
    v3=get_velocity(p1+0.5_WP*mydt*v2,myi,myj,myk)
    v4=get_velocity(p1+       mydt*v3,myi,myj,myk)
    p2=p1+mydt/6.0_WP*(v1+2.0_WP*v2+2.0_WP*v3+v4)
  end function project
  
  ! ================================================ !
  ! Interpolate velocity to pos near cell (i0,j0,k0) !
  ! ================================================ !
  function get_velocity(pos,i0,j0,k0) result(vel)
    implicit none
    real(WP), dimension(3) :: vel
    real(WP), dimension(3), intent(in) :: pos
    integer, intent(in) :: i0,j0,k0
    integer :: i,j,k
    real(WP) :: wx1,wy1,wz1
    real(WP) :: wx2,wy2,wz2
    real(WP), dimension(2,2,2) :: ww
    
    ! ! Interpolate Uface velocity ------------------------------
    ! ! Find right i index
    ! i=max(min(imaxo_-1,i0),imino_)
    ! do while (pos(1)-xm(i  ).lt.0.0_WP.and.i  .gt.imino_); i=i-1; end do
    ! do while (pos(1)-xm(i+1).ge.0.0_WP.and.i+1.lt.imaxo_); i=i+1; end do
    ! ! Find right j index
    ! j=max(min(jmaxo_-1,j0),jmino_)
    ! do while (pos(2)-ym(j  ).lt.0.0_WP.and.j  .gt.jmino_); j=j-1; end do
    ! do while (pos(2)-ym(j+1).ge.0.0_WP.and.j+1.lt.jmaxo_); j=j+1; end do
    ! ! Find right k index
    ! k=max(min(kmaxo_-1,k0),kmino_)
    ! do while (pos(3)-zm(k  ).lt.0.0_WP.and.k  .gt.kmino_); k=k-1; end do
    ! do while (pos(3)-zm(k+1).ge.0.0_WP.and.k+1.lt.kmaxo_); k=k+1; end do
    ! ! Prepare tri-linear interpolation coefficients
    ! wx1=(pos(1)-xm(i))/(xm(i+1)-xm(i)); wx2=1.0_WP-wx1
    ! wy1=(pos(2)-ym(j))/(ym(j+1)-ym(j)); wy2=1.0_WP-wy1
    ! wz1=(pos(3)-zm(k))/(zm(k+1)-zm(k)); wz2=1.0_WP-wz1
    ! ! Recast in array form
    ! ww(1,1,1)=wx2*wy2*wz2
    ! ww(2,1,1)=wx1*wy2*wz2
    ! ww(1,2,1)=wx2*wy1*wz2
    ! ww(2,2,1)=wx1*wy1*wz2
    ! ww(1,1,2)=wx2*wy2*wz1
    ! ww(2,1,2)=wx1*wy2*wz1
    ! ww(1,2,2)=wx2*wy1*wz1
    ! ww(2,2,2)=wx1*wy1*wz1
    ! ! Tri-linear interpolation of U/V/W
    ! vel(1)=sum(ww*U(i:i+1,j:j+1,k:k+1))
    ! vel(2)=sum(ww*V(i:i+1,j:j+1,k:k+1))
    ! vel(3)=sum(ww*W(i:i+1,j:j+1,k:k+1))
    
    ! Interpolate Uface velocity ------------------------------
    ! Find right i index
    i=max(min(imaxo_-1,i0),imino_)
    do while (pos(1)-x (i  ).lt.0.0_WP.and.i  .gt.imino_); i=i-1; end do
    do while (pos(1)-x (i+1).ge.0.0_WP.and.i+1.lt.imaxo_); i=i+1; end do
    ! Find right j index
    j=max(min(jmaxo_-1,j0),jmino_)
    do while (pos(2)-ym(j  ).lt.0.0_WP.and.j  .gt.jmino_); j=j-1; end do
    do while (pos(2)-ym(j+1).ge.0.0_WP.and.j+1.lt.jmaxo_); j=j+1; end do
    ! Find right k index
    k=max(min(kmaxo_-1,k0),kmino_)
    do while (pos(3)-zm(k  ).lt.0.0_WP.and.k  .gt.kmino_); k=k-1; end do
    do while (pos(3)-zm(k+1).ge.0.0_WP.and.k+1.lt.kmaxo_); k=k+1; end do
    ! Prepare tri-linear interpolation coefficients
    wx1=(pos(1)-x (i))/(x (i+1)-x (i)); wx2=1.0_WP-wx1
    wy1=(pos(2)-ym(j))/(ym(j+1)-ym(j)); wy2=1.0_WP-wy1
    wz1=(pos(3)-zm(k))/(zm(k+1)-zm(k)); wz2=1.0_WP-wz1
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
    vel(1)=sum(ww*Uface(i:i+1,j:j+1,k:k+1))
    
    ! Interpolate Vface velocity ------------------------------
    ! Find right i index
    i=max(min(imaxo_-1,i0),imino_)
    do while (pos(1)-xm(i  ).lt.0.0_WP.and.i  .gt.imino_); i=i-1; end do
    do while (pos(1)-xm(i+1).ge.0.0_WP.and.i+1.lt.imaxo_); i=i+1; end do
    ! Find right j index
    j=max(min(jmaxo_-1,j0),jmino_)
    do while (pos(2)-y (j  ).lt.0.0_WP.and.j  .gt.jmino_); j=j-1; end do
    do while (pos(2)-y (j+1).ge.0.0_WP.and.j+1.lt.jmaxo_); j=j+1; end do
    ! Find right k index
    k=max(min(kmaxo_-1,k0),kmino_)
    do while (pos(3)-zm(k  ).lt.0.0_WP.and.k  .gt.kmino_); k=k-1; end do
    do while (pos(3)-zm(k+1).ge.0.0_WP.and.k+1.lt.kmaxo_); k=k+1; end do
    ! Prepare tri-linear interpolation coefficients
    wx1=(pos(1)-xm(i))/(xm(i+1)-xm(i)); wx2=1.0_WP-wx1
    wy1=(pos(2)-y (j))/(y (j+1)-y (j)); wy2=1.0_WP-wy1
    wz1=(pos(3)-zm(k))/(zm(k+1)-zm(k)); wz2=1.0_WP-wz1
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
    vel(2)=sum(ww*Vface(i:i+1,j:j+1,k:k+1))
        
    ! Interpolate Wface velocity ------------------------------
    ! Find right i index
    i=max(min(imaxo_-1,i0),imino_)
    do while (pos(1)-xm(i  ).lt.0.0_WP.and.i  .gt.imino_); i=i-1; end do
    do while (pos(1)-xm(i+1).ge.0.0_WP.and.i+1.lt.imaxo_); i=i+1; end do
    ! Find right j index
    j=max(min(jmaxo_-1,j0),jmino_)
    do while (pos(2)-ym(j  ).lt.0.0_WP.and.j  .gt.jmino_); j=j-1; end do
    do while (pos(2)-ym(j+1).ge.0.0_WP.and.j+1.lt.jmaxo_); j=j+1; end do
    ! Find right k index
    k=max(min(kmaxo_-1,k0),kmino_)
    do while (pos(3)-z (k  ).lt.0.0_WP.and.k  .gt.kmino_); k=k-1; end do
    do while (pos(3)-z (k+1).ge.0.0_WP.and.k+1.lt.kmaxo_); k=k+1; end do
    ! Prepare tri-linear interpolation coefficients
    wx1=(pos(1)-xm(i))/(xm(i+1)-xm(i)); wx2=1.0_WP-wx1
    wy1=(pos(2)-ym(j))/(ym(j+1)-ym(j)); wy2=1.0_WP-wy1
    wz1=(pos(3)-z (k))/(z (k+1)-z (k)); wz2=1.0_WP-wz1
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
    vel(3)=sum(ww*Wface(i:i+1,j:j+1,k:k+1))
    
    return
  end function get_velocity
  
  ! =============================================== !
  ! Cut tet by computational mesh to compute fluxes !
  ! =============================================== !
  recursive function tet2flux(tet,ind) result(fluxes)
    implicit none
    real(WP), dimension(3,4), intent(in) :: tet
    integer,  dimension(3,4), intent(in) :: ind
    real(WP), dimension(13) :: fluxes
    integer :: dir,cut_ind
    integer :: n,nn,case
    integer :: v1,v2
    real(WP), dimension(4) :: d
    real(WP), dimension(3,8) :: vert
    integer,  dimension(3,8,2) :: vert_ind
    real(WP) :: mu
    real(WP), dimension(3,4) :: newtet
    integer,  dimension(3,4) :: newind
    
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
       fluxes=tet2flux_plic(tet,ind(1,1),ind(2,1),ind(3,1))    
       return
    end if
    
    ! Find case of cut
    case=1+int(0.5_WP+sign(0.5_WP,d(1)))+&
         2*int(0.5_WP+sign(0.5_WP,d(2)))+&
         4*int(0.5_WP+sign(0.5_WP,d(3)))+&
         8*int(0.5_WP+sign(0.5_WP,d(4)))   
    
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
       mu=min(1.0_WP,max(0.0_WP,-d(v1)/(sign(abs(d(v2)-d(v1))+tiny(1.0_WP),d(v2)-d(v1)))))
       vert(:,4+n)=(1.0_WP-mu)*vert(:,v1)+mu*vert(:,v2)
       ! Get index for interpolated vertex
       vert_ind(:,4+n,1)=get_indices(vert(:,4+n),vert_ind(:,v1,1))
       ! Enforce boundedness
       vert_ind(:,4+n,1)=max(vert_ind(:,4+n,1),min(vert_ind(:,v1,1),vert_ind(:,v2,1)))
       vert_ind(:,4+n,1)=min(vert_ind(:,4+n,1),max(vert_ind(:,v1,1),vert_ind(:,v2,1)))
       ! Set +/- indices in cut direction
       vert_ind(:,4+n,2)=vert_ind(:,4+n,1)
       vert_ind(dir,4+n,1)=cut_ind-1
       vert_ind(dir,4+n,2)=cut_ind
    end do
    
    ! Create new tets
    fluxes=0.0_WP
    do n=1,cut_ntets(case)
       do nn=1,4
          newtet(:,nn)=vert    (:,cut_vtet(nn,n,case))
          newind(:,nn)=vert_ind(:,cut_vtet(nn,n,case),cut_side(n,case))
       end do
       ! Check for zero-volume tet
       if (abs(tet_vol(newtet)).lt.volume_epsilon) cycle
       ! Cut by next plane
       fluxes=fluxes+tet2flux(newtet,newind)
    end do
    
    return
  end function tet2flux
  
  ! ========================= !
  ! Cut tet by PLIC interface !
  ! ========================= !
  function tet2flux_plic(tet,i,j,k) result(fluxes)
    implicit none
    real(WP), dimension(3,4), intent(in) :: tet
    integer, intent(in) :: i,j,k
    real(WP), dimension(13) :: fluxes
    integer :: n,case,v1,v2
    real(WP) :: mu,my_vol,my_rho,my_rhoE,my_rhoU,my_rhoV,my_rhoW,my_PA
    real(WP), dimension(4) :: d
    real(WP), dimension(3,8) :: vert
    real(WP), dimension(3) :: a,b,c,bary
    
    ! Check if wall
    if (vol(i,j,k).eq.0.0_WP) then
       fluxes=0.0_WP
       return
    end if
    
    ! Cut by PLIC
    d(:)=oldxnorm(i,j,k)*tet(1,:)+oldynorm(i,j,k)*tet(2,:)+oldznorm(i,j,k)*tet(3,:)-olddist(i,j,k)
    
    ! Find cut case
    case=1+int(0.5_WP+sign(0.5_WP,d(1)))+&
         2*int(0.5_WP+sign(0.5_WP,d(2)))+&
         4*int(0.5_WP+sign(0.5_WP,d(3)))+&
         8*int(0.5_WP+sign(0.5_WP,d(4)))
    
    ! Copy vertices
    do n=1,4
       vert(:,n)=tet(:,n)
    end do
    
    ! Create interpolated vertices on cut plane
    do n=1,cut_nvert(case)
       v1=cut_v1(n,case); v2=cut_v2(n,case)
       mu=min(1.0_WP,max(0.0_WP,-d(v1)/(sign(abs(d(v2)-d(v1))+tiny(1.0_WP),d(v2)-d(v1)))))
       vert(:,4+n)=(1.0_WP-mu)*vert(:,v1)+mu*vert(:,v2)
    end do
    
    ! Zero fluxes
    fluxes=0.0_WP
        
    ! Analyze gas tets
    do n=1,cut_nntet(case)-1
       ! Compute volume
       a=vert(:,cut_vtet(1,n,case))-vert(:,cut_vtet(4,n,case)); b=vert(:,cut_vtet(2,n,case))-vert(:,cut_vtet(4,n,case)); c=vert(:,cut_vtet(3,n,case))-vert(:,cut_vtet(4,n,case))
       my_vol=abs(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
       ! Compute barycenter
       bary=0.25_WP*(vert(:,cut_vtet(1,n,case))+vert(:,cut_vtet(2,n,case))+vert(:,cut_vtet(3,n,case))+vert(:,cut_vtet(4,n,case)))-Gbary(:,i,j,k)
       ! Compute tet variables
       my_rho =oldGrho (i,j,k)+sum(gradGrho (:,i,j,k)*bary(:))
       my_rhoE=oldGrhoE(i,j,k)+sum(gradGrhoE(:,i,j,k)*bary(:))
       my_rhoU=oldGrhoU(i,j,k)+sum(gradGrhoU(:,i,j,k)*bary(:))
       my_rhoV=oldGrhoV(i,j,k)+sum(gradGrhoV(:,i,j,k)*bary(:))
       my_rhoW=oldGrhoW(i,j,k)+sum(gradGrhoW(:,i,j,k)*bary(:))
       my_PA  =oldPA   (i,j,k)!+sum(gradPA   (:,i,j,k)*bary(:))
       ! Quantities in gas phase
       fluxes( 1)=fluxes( 1)+my_vol             ! Total volume
       fluxes( 3)=fluxes( 3)+my_vol*my_rho      ! Gas mass
       fluxes( 4)=fluxes( 4)+my_vol*my_rhoE     ! Gas energy
       fluxes( 5)=fluxes( 5)+my_vol*my_rhoU     ! Gas momentum x
       fluxes( 6)=fluxes( 6)+my_vol*my_rhoV     ! Gas momentum y
       fluxes( 7)=fluxes( 7)+my_vol*my_rhoW     ! Gas momentum z
       fluxes(13)=fluxes(13)+my_vol*my_PA       ! Gas pressure
    end do
    
    ! Analyze liquid tets
    do n=cut_ntets(case),cut_nntet(case),-1
       ! Compute volume
       a=vert(:,cut_vtet(1,n,case))-vert(:,cut_vtet(4,n,case)); b=vert(:,cut_vtet(2,n,case))-vert(:,cut_vtet(4,n,case)); c=vert(:,cut_vtet(3,n,case))-vert(:,cut_vtet(4,n,case))
       my_vol=abs(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
       ! Compute barycenter
       bary=0.25_WP*(vert(:,cut_vtet(1,n,case))+vert(:,cut_vtet(2,n,case))+vert(:,cut_vtet(3,n,case))+vert(:,cut_vtet(4,n,case)))-Lbary(:,i,j,k)
       ! Compute tet variables
       my_rho =oldLrho (i,j,k)+sum(gradLrho (:,i,j,k)*bary(:))
       my_rhoE=oldLrhoE(i,j,k)+sum(gradLrhoE(:,i,j,k)*bary(:))
       my_rhoU=oldLrhoU(i,j,k)+sum(gradLrhoU(:,i,j,k)*bary(:))
       my_rhoV=oldLrhoV(i,j,k)+sum(gradLrhoV(:,i,j,k)*bary(:))
       my_rhoW=oldLrhoW(i,j,k)+sum(gradLrhoW(:,i,j,k)*bary(:))
       my_PA  =oldPA   (i,j,k)!+sum(gradPA   (:,i,j,k)*bary(:))
       ! Quantities in liquid phase
       fluxes( 1)=fluxes( 1)+my_vol             ! Total volume
       fluxes( 2)=fluxes( 2)+my_vol             ! Liquid volume
       fluxes( 8)=fluxes( 8)+my_vol*my_rho      ! Liquid mass
       fluxes( 9)=fluxes( 9)+my_vol*my_rhoE     ! Liquid energy
       fluxes(10)=fluxes(10)+my_vol*my_rhoU     ! Liquid momentum x
       fluxes(11)=fluxes(11)+my_vol*my_rhoV     ! Liquid momentum y
       fluxes(12)=fluxes(12)+my_vol*my_rhoW     ! Liquid momentum z
       fluxes(13)=fluxes(13)+my_vol*my_PA       ! Liquid pressure
    end do
    
    return
  end function tet2flux_plic
  
end module compressible_advect



! ========================================= !
! Initialize semi-Lagrangian advection code !
! ========================================= !
subroutine compressible_advect_init
  use compressible_advect
  use memory
  implicit none

  integer :: i,j,k
  real(WP) :: rho_l, rho_r
  
  ! Create and store reference small volume
  volume_epsilon=1.0e-15_WP*min_meshsize**3
  
  ! Link flux pointers to temprary storage
  F_VOL  =>tmp1
  F_VOF  =>tmp2
  F_Grho =>tmp3
  F_GrhoE=>tmp4
  F_GrhoU=>tmp5
  F_GrhoV=>tmp6
  F_GrhoW=>tmp7
  F_Lrho =>tmp8
  F_LrhoE=>tmp9
  F_LrhoU=>tmp10
  F_LrhoV=>tmp11
  F_LrhoW=>tmp12
  F_PA   =>tmp13
  
  ! Allocate variables for field reconstruction
  allocate(gradGrho (3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(gradGrhoE(3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(gradGrhoU(3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(gradGrhoV(3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(gradGrhoW(3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(gradLrho (3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(gradLrhoE(3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(gradLrhoU(3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(gradLrhoV(3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(gradLrhoW(3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  
  ! Allocate phase barycenters
  allocate(Gbary(3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(Lbary(3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  
  ! Perform PLIC reconstruction
  call compressible_plic_calc

  return
end subroutine compressible_advect_init


! ===================================== !
! Semi-Lagrangian, flux-based advection !
! ===================================== !
subroutine compressible_advect_step
  use compressible_advect
  use metric_generic
  use time_info
  implicit none
  
  integer :: n,i,j,k,index,i_flux,j_flux,k_flux
  real(WP), dimension(3,17) :: pt
  real(WP), dimension(3,4) :: tet
  integer,  dimension(3,4) :: ind
  real(WP), dimension(13) :: flux
  real(WP) :: volcorr

  real(WP) :: A_sponge,sigma_sponge,A_sponge_g_in,A_sponge_l_in,A_sponge_g_out,A_sponge_l_out,A_sponge_sound
  real(WP) :: pmax,pprime,wavenum,freq
  real(WP) :: rhoprime,ubase,uprime,pressurebase,pressure,vbase,vprime
  real(WP) :: rhosponge,rhoUsponge,rhoVsponge,rhoEsponge,xloc1,xloc2,yloc1,yloc2,lambda
  real(WP) :: Grhosponge,GrhoUsponge,GrhoVsponge,GrhoEsponge,Lrhosponge,LrhoUsponge,LrhoVsponge,LrhoEsponge
  real(WP) :: sponge_length_y

  xloc1=0.01
  xloc2=0.04
  A_sponge_g_in=-0e4
  A_sponge_l_in=-10e4
  A_sponge_g_out=0.0_WP
  A_sponge_l_out=0.0_WP

  !pmax=11800.0_WP/2.0_WP
  pmax=8000.0_WP
  freq=1015.0_WP
  wavenum=2.0_WP*3.1415_WP/(340.0_WP/freq)
  yloc1=-0.050_WP
  yloc2=0.050_WP
  sponge_length_y=0.025_WP
  lambda=340.0_WP/freq

  A_sponge_sound=-300000.0_WP
  !A_sponge_sound=0.0_WP

  ! Loop over the domain and compute fluxes using semi-Lagrangian algorithm
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           
           ! Construct and project full cell, pt(:,17) is for correction tets
           pt(:, 9)=[x(i+1),y(j  ),z(k  )]; pt(:,1)=project(pt(:, 9),-dt,i,j,k)
           pt(:,10)=[x(i+1),y(j  ),z(k+1)]; pt(:,2)=project(pt(:,10),-dt,i,j,k)
           pt(:,11)=[x(i+1),y(j+1),z(k+1)]; pt(:,3)=project(pt(:,11),-dt,i,j,k)
           pt(:,12)=[x(i+1),y(j+1),z(k  )]; pt(:,4)=project(pt(:,12),-dt,i,j,k)
           pt(:,13)=[x(i  ),y(j  ),z(k  )]; pt(:,5)=project(pt(:,13),-dt,i,j,k)
           pt(:,14)=[x(i  ),y(j  ),z(k+1)]; pt(:,6)=project(pt(:,14),-dt,i,j,k)
           pt(:,15)=[x(i  ),y(j+1),z(k+1)]; pt(:,7)=project(pt(:,15),-dt,i,j,k)
           pt(:,16)=[x(i  ),y(j+1),z(k  )]; pt(:,8)=project(pt(:,16),-dt,i,j,k)
           
           ! Advect the whole cell
           flux=0.0_WP
           do n=1,6
              ! Build tet and corresponding indices
              tet(:,1)=pt(:,tet_map2(1,n)); ind(:,1)=get_indices(tet(:,1),[i,j,k])
              tet(:,2)=pt(:,tet_map2(2,n)); ind(:,2)=get_indices(tet(:,2),[i,j,k])
              tet(:,3)=pt(:,tet_map2(3,n)); ind(:,3)=get_indices(tet(:,3),[i,j,k])
              tet(:,4)=pt(:,tet_map2(4,n)); ind(:,4)=get_indices(tet(:,4),[i,j,k])
              ! Add corresponding flux
              flux=flux+tet_sign(tet)*tet2flux(tet,ind)
           end do
           ! Store update
           F_VOL  (i,j,k)=flux( 1)
           F_VOF  (i,j,k)=flux( 2)
           F_Grho (i,j,k)=flux( 3)
           F_GrhoE(i,j,k)=flux( 4)
           F_GrhoU(i,j,k)=flux( 5)
           F_GrhoV(i,j,k)=flux( 6)
           F_GrhoW(i,j,k)=flux( 7)
           F_Lrho (i,j,k)=flux( 8)
           F_LrhoE(i,j,k)=flux( 9)
           F_LrhoU(i,j,k)=flux(10)
           F_LrhoV(i,j,k)=flux(11)
           F_LrhoW(i,j,k)=flux(12)
           F_PA   (i,j,k)=flux(13)
           
           ! Correction for left x(i)
           volcorr=dt*Uface(i,j,k)*dy(j)*dz(k)
           do n=1,6
              tet(:,1) = pt(:,tet_mapxm(1,n))
              tet(:,2) = pt(:,tet_mapxm(2,n))
              tet(:,3) = pt(:,tet_mapxm(3,n))
              tet(:,4) = pt(:,tet_mapxm(4,n))
              volcorr=volcorr-tet_sign(tet)*abs(tet_vol(tet))
           end do
           pt(2,17)=0.25_WP*sum(pt(2,5:8))
           pt(3,17)=0.25_WP*sum(pt(3,5:8))
           pt(2,5:8)=pt(2,5:8)-pt(2,17)
           pt(3,5:8)=pt(3,5:8)-pt(3,17)
           pt(1,17)=(pt(3,6)*(pt(1,8)*pt(2,7) - pt(2,8)*pt(1,7)) - 6.0_WP*volcorr - &
                     pt(2,6)*(pt(1,8)*pt(3,7) - pt(3,8)*pt(1,7)) + pt(1,6)*(pt(2,8)*pt(3,7) - &
                     pt(3,8)*pt(2,7)) + pt(3,5)*(pt(1,8)*pt(2,6) - pt(2,8)*pt(1,6)) - &
                     pt(2,5)*(pt(1,8)*pt(3,6) - pt(3,8)*pt(1,6)) + pt(1,5)*(pt(2,8)*pt(3,6) - &
                     pt(3,8)*pt(2,6)))/(pt(2,8)*pt(3,7) - pt(3,8)*pt(2,7) + pt(2,8)*pt(3,6) - &
                     pt(3,8)*pt(2,6) - pt(3,6)*(pt(2,8) - pt(2,7)) + pt(2,6)*(pt(3,8) - pt(3,7)) - &
                     pt(3,5)*(pt(2,8) - pt(2,6)) + pt(2,5)*(pt(3,8) - pt(3,6)))
           pt(2,5:8)=pt(2,5:8)+pt(2,17)
           pt(3,5:8)=pt(3,5:8)+pt(3,17)
           ! Cut each correction tet recursively and calculate flux
           flux=0.0_WP
           i_flux=i-int(0.5_WP*(sign(1.0_WP,Uface(i,j,k))+1.0_WP))
           do n=7,8
              ! Build tet and corresponding indices
              tet(:,1)=pt(:,tet_mapxm(1,n)); ind(:,1)=get_indices(tet(:,1),[i_flux,j,k])
              tet(:,2)=pt(:,tet_mapxm(2,n)); ind(:,2)=get_indices(tet(:,2),[i_flux,j,k])
              tet(:,3)=pt(:,tet_mapxm(3,n)); ind(:,3)=get_indices(tet(:,3),[i_flux,j,k])
              tet(:,4)=pt(:,tet_mapxm(4,n)); ind(:,4)=get_indices(tet(:,4),[i_flux,j,k])
              ! Add corresponding flux
              flux=flux-tet_sign(tet)*tet2flux(tet,ind)
           end do
           ! Store update for right cell from face
           F_VOL  (i,j,k)=F_VOL(i,j,k)   -flux( 1)
           F_VOF  (i,j,k)=F_VOF(i,j,k)   -flux( 2)
           F_Grho (i,j,k)=F_Grho(i,j,k)  -flux( 3)
           F_GrhoE(i,j,k)=F_GrhoE(i,j,k) -flux( 4)
           F_GrhoU(i,j,k)=F_GrhoU(i,j,k) -flux( 5)
           F_GrhoV(i,j,k)=F_GrhoV(i,j,k) -flux( 6)
           F_GrhoW(i,j,k)=F_GrhoW(i,j,k) -flux( 7)
           F_Lrho (i,j,k)=F_Lrho(i,j,k)  -flux( 8)
           F_LrhoE(i,j,k)=F_LrhoE(i,j,k) -flux( 9)
           F_LrhoU(i,j,k)=F_LrhoU(i,j,k) -flux(10)
           F_LrhoV(i,j,k)=F_LrhoV(i,j,k) -flux(11)
           F_LrhoW(i,j,k)=F_LrhoW(i,j,k) -flux(12)
           F_PA   (i,j,k)=F_PA(i,j,k)    -flux(13)
           ! Store update for left cell from face
           F_VOL  (i-1,j,k)=F_VOL  (i-1,j,k) +flux( 1)
           F_VOF  (i-1,j,k)=F_VOF  (i-1,j,k) +flux( 2)
           F_Grho (i-1,j,k)=F_Grho (i-1,j,k) +flux( 3)
           F_GrhoE(i-1,j,k)=F_GrhoE(i-1,j,k) +flux( 4)
           F_GrhoU(i-1,j,k)=F_GrhoU(i-1,j,k) +flux( 5)
           F_GrhoV(i-1,j,k)=F_GrhoV(i-1,j,k) +flux( 6)
           F_GrhoW(i-1,j,k)=F_GrhoW(i-1,j,k) +flux( 7)
           F_Lrho (i-1,j,k)=F_Lrho (i-1,j,k) +flux( 8)
           F_LrhoE(i-1,j,k)=F_LrhoE(i-1,j,k) +flux( 9)
           F_LrhoU(i-1,j,k)=F_LrhoU(i-1,j,k) +flux(10)
           F_LrhoV(i-1,j,k)=F_LrhoV(i-1,j,k) +flux(11)
           F_LrhoW(i-1,j,k)=F_LrhoW(i-1,j,k) +flux(12)
           F_PA   (i-1,j,k)=F_PA   (i-1,j,k) +flux(13)
           
           ! Correction for bottom y(j)
           volcorr=-dt*Vface(i,j,k)*dz(k)*dx(i)
           do n=1,6
              tet(:,1) = pt(:,tet_mapym(1,n))
              tet(:,2) = pt(:,tet_mapym(2,n))
              tet(:,3) = pt(:,tet_mapym(3,n))
              tet(:,4) = pt(:,tet_mapym(4,n))
              volcorr=volcorr-tet_sign(tet)*abs(tet_vol(tet))
           end do
           pt(1,17)=0.25_WP*(pt(1,1)+pt(1,2)+pt(1,5)+pt(1,6))
           pt(3,17)=0.25_WP*(pt(3,1)+pt(3,2)+pt(3,5)+pt(3,6))
           pt(1,1:2)=pt(1,1:2)-pt(1,17);pt(1,5:6)=pt(1,5:6)-pt(1,17);
           pt(3,1:2)=pt(3,1:2)-pt(3,17);pt(3,5:6)=pt(3,5:6)-pt(3,17);
           pt(2,17)=-(pt(3,6)*(pt(1,1)*pt(2,2) - pt(2,1)*pt(1,2)) - 6.0_WP*volcorr - &
                      pt(2,6)*(pt(1,1)*pt(3,2) - pt(3,1)*pt(1,2)) + pt(1,6)*(pt(2,1)*pt(3,2) - &
                      pt(3,1)*pt(2,2)) + pt(3,5)*(pt(1,1)*pt(2,6) - pt(2,1)*pt(1,6)) - &
                      pt(2,5)*(pt(1,1)*pt(3,6) - pt(3,1)*pt(1,6)) + pt(1,5)*(pt(2,1)*pt(3,6) - &
                      pt(3,1)*pt(2,6)))/(pt(1,1)*pt(3,2) - pt(3,1)*pt(1,2) + pt(1,1)*pt(3,6) - &
                      pt(3,1)*pt(1,6) - pt(3,6)*(pt(1,1) - pt(1,2)) + pt(1,6)*(pt(3,1) - &
                      pt(3,2)) - pt(3,5)*(pt(1,1) - pt(1,6)) + pt(1,5)*(pt(3,1) - pt(3,6)))
           pt(1,1:2)=pt(1,1:2)+pt(1,17);pt(1,5:6)=pt(1,5:6)+pt(1,17);
           pt(3,1:2)=pt(3,1:2)+pt(3,17);pt(3,5:6)=pt(3,5:6)+pt(3,17);
           ! Cut each correction tet recursively and calculate flux
           flux=0.0_WP
           j_flux=j-int(0.5_WP*(sign(1.0_WP,Vface(i,j,k))+1.0_WP))
           do n=7,8
              ! Build tet and corresponding indices
              tet(:,1)=pt(:,tet_mapym(1,n)); ind(:,1)=get_indices(tet(:,1),[i,j_flux,k])
              tet(:,2)=pt(:,tet_mapym(2,n)); ind(:,2)=get_indices(tet(:,2),[i,j_flux,k])
              tet(:,3)=pt(:,tet_mapym(3,n)); ind(:,3)=get_indices(tet(:,3),[i,j_flux,k])
              tet(:,4)=pt(:,tet_mapym(4,n)); ind(:,4)=get_indices(tet(:,4),[i,j_flux,k])
              ! Add corresponding flux
              flux=flux+tet_sign(tet)*tet2flux(tet,ind)
           end do
           ! Store update for right cell from face
           F_VOL  (i,j,k)=F_VOL(i,j,k)   -flux( 1)
           F_VOF  (i,j,k)=F_VOF(i,j,k)   -flux( 2)
           F_Grho (i,j,k)=F_Grho(i,j,k)  -flux( 3)
           F_GrhoE(i,j,k)=F_GrhoE(i,j,k) -flux( 4)
           F_GrhoU(i,j,k)=F_GrhoU(i,j,k) -flux( 5)
           F_GrhoV(i,j,k)=F_GrhoV(i,j,k) -flux( 6)
           F_GrhoW(i,j,k)=F_GrhoW(i,j,k) -flux( 7)
           F_Lrho (i,j,k)=F_Lrho(i,j,k)  -flux( 8)
           F_LrhoE(i,j,k)=F_LrhoE(i,j,k) -flux( 9)
           F_LrhoU(i,j,k)=F_LrhoU(i,j,k) -flux(10)
           F_LrhoV(i,j,k)=F_LrhoV(i,j,k) -flux(11)
           F_LrhoW(i,j,k)=F_LrhoW(i,j,k) -flux(12)
           F_PA   (i,j,k)=F_PA(i,j,k)    -flux(13)
           ! Store update for left cell from face
           F_VOL  (i,j-1,k)=F_VOL  (i,j-1,k) +flux( 1)
           F_VOF  (i,j-1,k)=F_VOF  (i,j-1,k) +flux( 2)
           F_Grho (i,j-1,k)=F_Grho (i,j-1,k) +flux( 3)
           F_GrhoE(i,j-1,k)=F_GrhoE(i,j-1,k) +flux( 4)
           F_GrhoU(i,j-1,k)=F_GrhoU(i,j-1,k) +flux( 5)
           F_GrhoV(i,j-1,k)=F_GrhoV(i,j-1,k) +flux( 6)
           F_GrhoW(i,j-1,k)=F_GrhoW(i,j-1,k) +flux( 7)
           F_Lrho (i,j-1,k)=F_Lrho (i,j-1,k) +flux( 8)
           F_LrhoE(i,j-1,k)=F_LrhoE(i,j-1,k) +flux( 9)
           F_LrhoU(i,j-1,k)=F_LrhoU(i,j-1,k) +flux(10)
           F_LrhoV(i,j-1,k)=F_LrhoV(i,j-1,k) +flux(11)
           F_LrhoW(i,j-1,k)=F_LrhoW(i,j-1,k) +flux(12)
           F_PA   (i,j-1,k)=F_PA   (i,j-1,k) +flux(13)
           
           ! Correction for back z(k)
           volcorr=dt*Wface(i,j,k)*dx(i)*dy(j)
           do n=1,6
              tet(:,1) = pt(:,tet_mapzm(1,n))
              tet(:,2) = pt(:,tet_mapzm(2,n))
              tet(:,3) = pt(:,tet_mapzm(3,n))
              tet(:,4) = pt(:,tet_mapzm(4,n))
              volcorr=volcorr-tet_sign(tet)*abs(tet_vol(tet))
           end do
           pt(1,17)=0.25_WP*(pt(1,1)+pt(1,5)+pt(1,8)+pt(1,4))
           pt(2,17)=0.25_WP*(pt(2,1)+pt(2,5)+pt(2,8)+pt(2,4))
           ! Reference to known barycenter plane location
           pt(1,1)=pt(1,1)-pt(1,17); pt(1,5)=pt(1,5)-pt(1,17);
           pt(1,8)=pt(1,8)-pt(1,17); pt(1,4)=pt(1,4)-pt(1,17);
           pt(2,1)=pt(2,1)-pt(2,17); pt(2,5)=pt(2,5)-pt(2,17);
           pt(2,8)=pt(2,8)-pt(2,17); pt(2,4)=pt(2,4)-pt(2,17);
           pt(3,17)=(pt(3,5)*(pt(1,4)*pt(2,8) - pt(2,4)*pt(1,8)) - 6.0_WP*volcorr - &
                     pt(2,5)*(pt(1,4)*pt(3,8) - pt(3,4)*pt(1,8)) + pt(1,5)*(pt(2,4)*pt(3,8) - &
                     pt(3,4)*pt(2,8)) + pt(3,1)*(pt(1,4)*pt(2,5) - pt(2,4)*pt(1,5)) - &
                     pt(2,1)*(pt(1,4)*pt(3,5) - pt(3,4)*pt(1,5)) + pt(1,1)*(pt(2,4)*pt(3,5) - &
                     pt(3,4)*pt(2,5)))/(pt(1,4)*pt(2,8) - pt(2,4)*pt(1,8) + pt(1,4)*pt(2,5) - &
                     pt(2,4)*pt(1,5) - pt(2,5)*(pt(1,4) - pt(1,8)) + pt(1,5)*(pt(2,4) - pt(2,8)) - &
                     pt(2,1)*(pt(1,4) - pt(1,5)) + pt(1,1)*(pt(2,4) - pt(2,5)))
           ! Return points to oriignal value
           pt(1,1)=pt(1,1)+pt(1,17); pt(1,5)=pt(1,5)+pt(1,17);
           pt(1,8)=pt(1,8)+pt(1,17); pt(1,4)=pt(1,4)+pt(1,17);
           pt(2,1)=pt(2,1)+pt(2,17); pt(2,5)=pt(2,5)+pt(2,17);
           pt(2,8)=pt(2,8)+pt(2,17); pt(2,4)=pt(2,4)+pt(2,17);
           ! Cut each tet recursively and calculate flux
           flux=0.0_WP
           k_flux=k-int(0.5_WP*(sign(1.0_WP,Wface(i,j,k))+1.0_WP))
           do n=7,8
              ! Build tet and corresponding indices
              tet(:,1)=pt(:,tet_mapzm(1,n)); ind(:,1)=get_indices(tet(:,1),[i,j,k_flux])
              tet(:,2)=pt(:,tet_mapzm(2,n)); ind(:,2)=get_indices(tet(:,2),[i,j,k_flux])
              tet(:,3)=pt(:,tet_mapzm(3,n)); ind(:,3)=get_indices(tet(:,3),[i,j,k_flux])
              tet(:,4)=pt(:,tet_mapzm(4,n)); ind(:,4)=get_indices(tet(:,4),[i,j,k_flux])
              ! Add corresponding flux
              flux=flux-tet_sign(tet)*tet2flux(tet,ind)
           end do
           ! Store update for right cell from face
           F_VOL  (i,j,k)=F_VOL(i,j,k)   -flux( 1)
           F_VOF  (i,j,k)=F_VOF(i,j,k)   -flux( 2)
           F_Grho (i,j,k)=F_Grho(i,j,k)  -flux( 3)
           F_GrhoE(i,j,k)=F_GrhoE(i,j,k) -flux( 4)
           F_GrhoU(i,j,k)=F_GrhoU(i,j,k) -flux( 5)
           F_GrhoV(i,j,k)=F_GrhoV(i,j,k) -flux( 6)
           F_GrhoW(i,j,k)=F_GrhoW(i,j,k) -flux( 7)
           F_Lrho (i,j,k)=F_Lrho(i,j,k)  -flux( 8)
           F_LrhoE(i,j,k)=F_LrhoE(i,j,k) -flux( 9)
           F_LrhoU(i,j,k)=F_LrhoU(i,j,k) -flux(10)
           F_LrhoV(i,j,k)=F_LrhoV(i,j,k) -flux(11)
           F_LrhoW(i,j,k)=F_LrhoW(i,j,k) -flux(12)
           F_PA   (i,j,k)=F_PA(i,j,k)    -flux(13)
           ! Store update for left cell from face
           F_VOL  (i,j,k-1)=F_VOL  (i,j,k-1) +flux( 1)
           F_VOF  (i,j,k-1)=F_VOF  (i,j,k-1) +flux( 2)
           F_Grho (i,j,k-1)=F_Grho (i,j,k-1) +flux( 3)
           F_GrhoE(i,j,k-1)=F_GrhoE(i,j,k-1) +flux( 4)
           F_GrhoU(i,j,k-1)=F_GrhoU(i,j,k-1) +flux( 5)
           F_GrhoV(i,j,k-1)=F_GrhoV(i,j,k-1) +flux( 6)
           F_GrhoW(i,j,k-1)=F_GrhoW(i,j,k-1) +flux( 7)
           F_Lrho (i,j,k-1)=F_Lrho (i,j,k-1) +flux( 8)
           F_LrhoE(i,j,k-1)=F_LrhoE(i,j,k-1) +flux( 9)
           F_LrhoU(i,j,k-1)=F_LrhoU(i,j,k-1) +flux(10)
           F_LrhoV(i,j,k-1)=F_LrhoV(i,j,k-1) +flux(11)
           F_LrhoW(i,j,k-1)=F_LrhoW(i,j,k-1) +flux(12)
           F_PA   (i,j,k-1)=F_PA   (i,j,k-1) +flux(13)
           
           if (i.eq.imax_) then
              ! Get correction for right face, since it won't be revisited
              volcorr=dt*Uface(i+1,j,k)*dy(j)*dz(k)
              do n=1,6
                 tet(:,1) = pt(:,tet_mapxp(1,n))
                 tet(:,2) = pt(:,tet_mapxp(2,n))
                 tet(:,3) = pt(:,tet_mapxp(3,n))
                 tet(:,4) = pt(:,tet_mapxp(4,n))
                 volcorr=volcorr-tet_sign(tet)*abs(tet_vol(tet))
              end do
              pt(2,17)=0.25_WP*sum(pt(2,1:4))
              pt(3,17)=0.25_WP*sum(pt(3,1:4))
              pt(2,1:4)=pt(2,1:4)-pt(2,17)
              pt(3,1:4)=pt(3,1:4)-pt(3,17)
              pt(1,17)=(pt(3,2)*(pt(1,4)*pt(2,3) - pt(2,4)*pt(1,3)) - 6.0_WP*volcorr - &
                        pt(2,2)*(pt(1,4)*pt(3,3) - pt(3,4)*pt(1,3)) + pt(1,2)*(pt(2,4)*pt(3,3) - &
                        pt(3,4)*pt(2,3)) + pt(3,1)*(pt(1,4)*pt(2,2) - pt(2,4)*pt(1,2)) - &
                        pt(2,1)*(pt(1,4)*pt(3,2) - pt(3,4)*pt(1,2)) + pt(1,1)*(pt(2,4)*pt(3,2) - &
                        pt(3,4)*pt(2,2)))/(pt(2,4)*pt(3,3) - pt(3,4)*pt(2,3) + pt(2,4)*pt(3,2) - &
                        pt(3,4)*pt(2,2) - pt(3,2)*(pt(2,4) - pt(2,3)) + pt(2,2)*(pt(3,4) - &
                        pt(3,3)) - pt(3,1)*(pt(2,4) - pt(2,2)) + pt(2,1)*(pt(3,4) - pt(3,2)))
              pt(2,1:4)=pt(2,1:4)+pt(2,17)
              pt(3,1:4)=pt(3,1:4)+pt(3,17)
              ! Cut each correction tet recursively and calculate flux
              flux=0.0_WP
              i_flux=i-int(0.5_WP*(sign(1.0_WP,Uface(i+1,j,k))+1.0_WP))
              do n=7,8
                 ! Build tet and corresponding indices
                 tet(:,1)=pt(:,tet_mapxp(1,n)); ind(:,1)=get_indices(tet(:,1),[i_flux+1,j,k])
                 tet(:,2)=pt(:,tet_mapxp(2,n)); ind(:,2)=get_indices(tet(:,2),[i_flux+1,j,k])
                 tet(:,3)=pt(:,tet_mapxp(3,n)); ind(:,3)=get_indices(tet(:,3),[i_flux+1,j,k])
                 tet(:,4)=pt(:,tet_mapxp(4,n)); ind(:,4)=get_indices(tet(:,4),[i_flux+1,j,k])
                 ! Add corresponding flux
                 flux=flux-tet_sign(tet)*tet2flux(tet,ind)
              end do
              ! Store update
              F_VOL  (i,j,k)=F_VOL(i,j,k)   +flux( 1)
              F_VOF  (i,j,k)=F_VOF(i,j,k)   +flux( 2)
              F_Grho (i,j,k)=F_Grho(i,j,k)  +flux( 3)
              F_GrhoE(i,j,k)=F_GrhoE(i,j,k) +flux( 4)
              F_GrhoU(i,j,k)=F_GrhoU(i,j,k) +flux( 5)
              F_GrhoV(i,j,k)=F_GrhoV(i,j,k) +flux( 6)
              F_GrhoW(i,j,k)=F_GrhoW(i,j,k) +flux( 7)
              F_Lrho (i,j,k)=F_Lrho(i,j,k)  +flux( 8)
              F_LrhoE(i,j,k)=F_LrhoE(i,j,k) +flux( 9)
              F_LrhoU(i,j,k)=F_LrhoU(i,j,k) +flux(10)
              F_LrhoV(i,j,k)=F_LrhoV(i,j,k) +flux(11)
              F_LrhoW(i,j,k)=F_LrhoW(i,j,k) +flux(12)
              F_PA   (i,j,k)=F_PA(i,j,k)    +flux(13)
           end if
           
           if (j.eq.jmax_) then
              volcorr=-dt*Vface(i,j+1,k)*dz(k)*dx(i)
              do n=1,6
                 tet(:,1) = pt(:,tet_mapyp(1,n))
                 tet(:,2) = pt(:,tet_mapyp(2,n))
                 tet(:,3) = pt(:,tet_mapyp(3,n))
                 tet(:,4) = pt(:,tet_mapyp(4,n))
                 volcorr=volcorr-tet_sign(tet)*abs(tet_vol(tet))
              end do
              pt(1,17)=0.25_WP*(pt(1,4)+pt(1,3)+pt(1,7)+pt(1,8))
              pt(3,17)=0.25_WP*(pt(3,4)+pt(3,3)+pt(3,7)+pt(3,8))
              ! Reference to known barycenter plane location
              pt(1,3:4)=pt(1,3:4)-pt(1,17);pt(1,7:8)=pt(1,7:8)-pt(1,17)
              pt(3,3:4)=pt(3,3:4)-pt(3,17);pt(3,7:8)=pt(3,7:8)-pt(3,17)
              pt(2,17)=-(pt(3,7)*(pt(1,4)*pt(2,3) - pt(2,4)*pt(1,3)) - 6.0_WP*volcorr - &
                         pt(2,7)*(pt(1,4)*pt(3,3) - pt(3,4)*pt(1,3)) + pt(1,7)*(pt(2,4)*pt(3,3) - &
                         pt(3,4)*pt(2,3)) + pt(3,8)*(pt(1,4)*pt(2,7) - pt(2,4)*pt(1,7)) - &
                         pt(2,8)*(pt(1,4)*pt(3,7) - pt(3,4)*pt(1,7)) + pt(1,8)*(pt(2,4)*pt(3,7) - &
                         pt(3,4)*pt(2,7)))/(pt(1,4)*pt(3,3) - pt(3,4)*pt(1,3) + pt(1,4)*pt(3,7) - &
                         pt(3,4)*pt(1,7) - pt(3,7)*(pt(1,4) - pt(1,3)) + pt(1,7)*(pt(3,4) - pt(3,3)) - &
                         pt(3,8)*(pt(1,4) - pt(1,7)) + pt(1,8)*(pt(3,4) - pt(3,7)))
              ! Return points to oriignal value
              pt(1,3:4)=pt(1,3:4)+pt(1,17);pt(1,7:8)=pt(1,7:8)+pt(1,17)
              pt(3,3:4)=pt(3,3:4)+pt(3,17);pt(3,7:8)=pt(3,7:8)+pt(3,17)
              ! Cut each correction tet recursively and calculate flux
              flux=0.0_WP
              j_flux=j-int(0.5_WP*(sign(1.0_WP,Vface(i,j+1,k))+1.0_WP))
              do n=7,8
                 ! Build tet and corresponding indices
                 tet(:,1)=pt(:,tet_mapyp(1,n)); ind(:,1)=get_indices(tet(:,1),[i,j_flux+1,k])
                 tet(:,2)=pt(:,tet_mapyp(2,n)); ind(:,2)=get_indices(tet(:,2),[i,j_flux+1,k])
                 tet(:,3)=pt(:,tet_mapyp(3,n)); ind(:,3)=get_indices(tet(:,3),[i,j_flux+1,k])
                 tet(:,4)=pt(:,tet_mapyp(4,n)); ind(:,4)=get_indices(tet(:,4),[i,j_flux+1,k])
                 ! Add corresponding flux
                 flux=flux+tet_sign(tet)*tet2flux(tet,ind)
              end do
              ! Store update
              F_VOL  (i,j,k)=F_VOL(i,j,k)   +flux( 1)
              F_VOF  (i,j,k)=F_VOF(i,j,k)   +flux( 2)
              F_Grho (i,j,k)=F_Grho(i,j,k)  +flux( 3)
              F_GrhoE(i,j,k)=F_GrhoE(i,j,k) +flux( 4)
              F_GrhoU(i,j,k)=F_GrhoU(i,j,k) +flux( 5)
              F_GrhoV(i,j,k)=F_GrhoV(i,j,k) +flux( 6)
              F_GrhoW(i,j,k)=F_GrhoW(i,j,k) +flux( 7)
              F_Lrho (i,j,k)=F_Lrho(i,j,k)  +flux( 8)
              F_LrhoE(i,j,k)=F_LrhoE(i,j,k) +flux( 9)
              F_LrhoU(i,j,k)=F_LrhoU(i,j,k) +flux(10)
              F_LrhoV(i,j,k)=F_LrhoV(i,j,k) +flux(11)
              F_LrhoW(i,j,k)=F_LrhoW(i,j,k) +flux(12)
              F_PA   (i,j,k)=F_PA(i,j,k)    +flux(13)
           end if
           
           if (k.eq.kmax_) then
              volcorr=dt*Wface(i,j,k+1)*dx(i)*dy(j)
              do n=1,6
                 tet(:,1) = pt(:,tet_mapzp(1,n))
                 tet(:,2) = pt(:,tet_mapzp(2,n))
                 tet(:,3) = pt(:,tet_mapzp(3,n))
                 tet(:,4) = pt(:,tet_mapzp(4,n))
                 volcorr=volcorr-tet_sign(tet)*abs(tet_vol(tet))
              end do
              pt(1,17)=0.25_WP*(pt(1,2)+pt(1,6)+pt(1,7)+pt(1,3))
              pt(2,17)=0.25_WP*(pt(2,2)+pt(2,6)+pt(2,7)+pt(2,3))
              ! Reference to known barycenter plane location
              pt(1,2:3)=pt(1,2:3)-pt(1,17);pt(1,6:7)=pt(1,6:7)-pt(1,17);
              pt(2,2:3)=pt(2,2:3)-pt(2,17);pt(2,6:7)=pt(2,6:7)-pt(2,17);
              pt(3,17)=(pt(3,6)*(pt(1,3)*pt(2,7) - pt(2,3)*pt(1,7)) - 6.0_WP*volcorr - &
                        pt(2,6)*(pt(1,3)*pt(3,7) - pt(3,3)*pt(1,7)) + pt(1,6)*(pt(2,3)*pt(3,7) - &
                        pt(3,3)*pt(2,7)) + pt(3,2)*(pt(1,3)*pt(2,6) - pt(2,3)*pt(1,6)) - &
                        pt(2,2)*(pt(1,3)*pt(3,6) - pt(3,3)*pt(1,6)) + pt(1,2)*(pt(2,3)*pt(3,6) - &
                        pt(3,3)*pt(2,6)))/(pt(1,3)*pt(2,7) - pt(2,3)*pt(1,7) + pt(1,3)*pt(2,6) - &
                        pt(2,3)*pt(1,6) - pt(2,6)*(pt(1,3) - pt(1,7)) + pt(1,6)*(pt(2,3) - &
                        pt(2,7)) - pt(2,2)*(pt(1,3) - pt(1,6)) + pt(1,2)*(pt(2,3) - pt(2,6)))
              ! Return points to oriignal value
              pt(1,2:3)=pt(1,2:3)+pt(1,17);pt(1,6:7)=pt(1,6:7)+pt(1,17);
              pt(2,2:3)=pt(2,2:3)+pt(2,17);pt(2,6:7)=pt(2,6:7)+pt(2,17);
              ! Cut each tet recursively and calculate flux
              flux=0.0_WP
              k_flux=k-int(0.5_WP*(sign(1.0_WP,Wface(i,j,k+1))+1.0_WP))
              do n=7,8
                 ! Build tet and corresponding indices
                 tet(:,1)=pt(:,tet_mapzp(1,n)); ind(:,1)=get_indices(tet(:,1),[i,j,k_flux+1])
                 tet(:,2)=pt(:,tet_mapzp(2,n)); ind(:,2)=get_indices(tet(:,2),[i,j,k_flux+1])
                 tet(:,3)=pt(:,tet_mapzp(3,n)); ind(:,3)=get_indices(tet(:,3),[i,j,k_flux+1])
                 tet(:,4)=pt(:,tet_mapzp(4,n)); ind(:,4)=get_indices(tet(:,4),[i,j,k_flux+1])
                 ! Add corresponding flux
                 flux=flux-tet_sign(tet)*tet2flux(tet,ind)
              end do
              ! Store update
              F_VOL  (i,j,k)=F_VOL(i,j,k)   +flux( 1)
              F_VOF  (i,j,k)=F_VOF(i,j,k)   +flux( 2)
              F_Grho (i,j,k)=F_Grho(i,j,k)  +flux( 3)
              F_GrhoE(i,j,k)=F_GrhoE(i,j,k) +flux( 4)
              F_GrhoU(i,j,k)=F_GrhoU(i,j,k) +flux( 5)
              F_GrhoV(i,j,k)=F_GrhoV(i,j,k) +flux( 6)
              F_GrhoW(i,j,k)=F_GrhoW(i,j,k) +flux( 7)
              F_Lrho (i,j,k)=F_Lrho(i,j,k)  +flux( 8)
              F_LrhoE(i,j,k)=F_LrhoE(i,j,k) +flux( 9)
              F_LrhoU(i,j,k)=F_LrhoU(i,j,k) +flux(10)
              F_LrhoV(i,j,k)=F_LrhoV(i,j,k) +flux(11)
              F_LrhoW(i,j,k)=F_LrhoW(i,j,k) +flux(12)
              F_PA   (i,j,k)=F_PA(i,j,k)    +flux(13)
           end if
           
        end do
     end do
  end do
 
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           
           ! Cycle wall cells
           if (vol(i,j,k).eq.0.0_WP) cycle
           
           ! Update VOF
           VOF(i,j,k)=F_VOF(i,j,k)/F_VOL(i,j,k)
	   !if(xm(i).le.xloc1)then
		!sigma_sponge=(xloc1-xm(i))/0.01_WP
		!VOF(i,j,k)=VOF(i,j,k)-100000.0*sigma_sponge**2*(VOF(i,j,k)-0.0_WP)
	   !endif

	   if(xm(i).ge.xloc2)then
		sigma_sponge=(xm(i)-xloc2)/0.01_WP
		VOF(i,j,k)=VOF(i,j,k)!-5.0*sigma_sponge**2*(VOF(i,j,k)-0.0_WP)
	   endif
	   if(ym(j).le.yloc1)then
		sigma_sponge=(-ym(j)+yloc1)/sponge_length_y
		VOF(i,j,k)=VOF(i,j,k)!-10000.0*sigma_sponge**2*(VOF(i,j,k)-0.0_WP)
	   endif
  	   if(ym(j).ge.yloc2)then
		sigma_sponge=(ym(j)-yloc2)/sponge_length_y	
		VOF(i,j,k)=VOF(i,j,k)!-10000.0*sigma_sponge**2*(VOF(i,j,k)-0.0_WP)
	   endif
           ! Update PA
           PA(i,j,k)=F_PA(i,j,k)/F_VOL(i,j,k)
           ! Update phase density, energy and momentum

	   call SpongeData(j,k,Grhosponge,GrhoUsponge,GrhoVsponge,GrhoEsponge,Lrhosponge,LrhoUsponge,LrhoVsponge,LrhoEsponge)

           if      (VOF(i,j,k).lt.VOFlo) then	 	

	     if(xm(i).gt.xloc1.and.ym(j).gt.yloc1.and.ym(j).lt.yloc2)then	
	     !if(ym(j).gt.yloc1.and.ym(j).lt.yloc2)then	
		A_sponge=0.0
		sigma_sponge=0.0
		Grho (i,j,k)=F_Grho (i,j,k)/vol(i,j,k)+dt*A_sponge*sigma_sponge**2*(Grho(i,j,k)-Grhosponge)
                GrhoE(i,j,k)=F_GrhoE(i,j,k)/vol(i,j,k)+dt*A_sponge*sigma_sponge**2*(GrhoE(i,j,k)-GrhoEsponge)
                GrhoU(i,j,k)=F_GrhoU(i,j,k)/vol(i,j,k)+dt*A_sponge*sigma_sponge**2*(GrhoU(i,j,k)-GrhoUsponge)
                GrhoV(i,j,k)=F_GrhoV(i,j,k)/vol(i,j,k)+dt*A_sponge*sigma_sponge**2*(GrhoV(i,j,k)-GrhoVsponge)
                GrhoW(i,j,k)=F_GrhoW(i,j,k)/vol(i,j,k)!+dt*A_sponge*sigma_sponge**2*(GrhoW(i,j,k)-0.0_WP)

	     elseif(xm(i).ge.xloc2)then
		!A_sponge=A_sponge_g_out
		!sigma_sponge=(xm(i)-xloc2)/(0.01_WP)
		!Grho (i,j,k)=F_Grho (i,j,k)/vol(i,j,k)!+dt*A_sponge*sigma_sponge**2*(Grho(i,j,k)-Grhosponge)
                !GrhoE(i,j,k)=F_GrhoE(i,j,k)/vol(i,j,k)!+dt*A_sponge*sigma_sponge**2*(GrhoE(i,j,k)-GrhoEsponge)
                !GrhoU(i,j,k)=F_GrhoU(i,j,k)/vol(i,j,k)!+dt*A_sponge*sigma_sponge**2*(GrhoU(i,j,k)-GrhoUsponge)
                !GrhoV(i,j,k)=F_GrhoV(i,j,k)/vol(i,j,k)!+dt*A_sponge*sigma_sponge**2*(GrhoV(i,j,k)-GrhoVsponge)
                !GrhoW(i,j,k)=F_GrhoW(i,j,k)/vol(i,j,k)!+dt*A_sponge*sigma_sponge**2*(GrhoW(i,j,k)-0.0_WP)

	     elseif(xm(i).le.xloc1)then
		A_sponge=A_sponge_l_in
                sigma_sponge=(xloc1-xm(i))/(0.01_WP)
	        Grho (i,j,k)=F_Grho (i,j,k)/vol(i,j,k)+dt*A_sponge*sigma_sponge**2*(Grho(i,j,k)-Grhosponge)
                GrhoE(i,j,k)=F_GrhoE(i,j,k)/vol(i,j,k)+dt*A_sponge*sigma_sponge**2*(GrhoE(i,j,k)-GrhoEsponge)
                GrhoU(i,j,k)=F_GrhoU(i,j,k)/vol(i,j,k)+dt*A_sponge*sigma_sponge**2*(GrhoU(i,j,k)-GrhoUsponge)
                GrhoV(i,j,k)=F_GrhoV(i,j,k)/vol(i,j,k)+dt*A_sponge*sigma_sponge**2*(GrhoV(i,j,k)-GrhoVsponge)
                GrhoW(i,j,k)=F_GrhoW(i,j,k)/vol(i,j,k)+dt*A_sponge*sigma_sponge**2*(GrhoW(i,j,k)-0.0_WP)

	     endif

	      if(ym(j).le.yloc1)then

              sigma_sponge=(-ym(j)+yloc1)/sponge_length_y
              ! Create sponge target
              !pprime=-1.0d0*pmax*sin(3.1415d0/4.0d0+wavenum*(-ym(j)-lambda/2.0d0)+2.0_WP*3.1415_WP*freq*time)
              pprime=-1.0d0*pmax*sin(wavenum*(-ym(j)-lambda/2.0d0)+2.0_WP*3.1415_WP*freq*time)
              rhoprime=pprime/sound_speed(i,j,k)**2
              vprime=pprime/(Grho(i,j,k)*sound_speed(i,j,k))

              rhosponge=1.226+rhoprime
              vbase=0.0!GrhoU(i,j,k)/Grho(i,j,k)
              rhoVsponge=rhosponge*(vbase+vprime)
              pressurebase=101325.0_WP
              pressure=pressurebase+pprime
              rhoEsponge=pressure/0.4_WP+0.5_WP*rhoVsponge**2/rhosponge
              rhoUsponge=0.0_WP

              VOF  (i,j,k)=0.0_WP
              Grho (i,j,k)=F_Grho (i,j,k)/vol(i,j,k)+dt*A_sponge_sound*sigma_sponge**2*(Grho(i,j,k)-rhosponge)
              GrhoE(i,j,k)=F_GrhoE(i,j,k)/vol(i,j,k)+dt*A_sponge_sound*sigma_sponge**2*(GrhoE(i,j,k)-rhoEsponge)
              GrhoU(i,j,k)=F_GrhoU(i,j,k)/vol(i,j,k)+dt*A_sponge_sound*sigma_sponge**2*(GrhoU(i,j,k)-rhoUsponge)
              GrhoV(i,j,k)=F_GrhoV(i,j,k)/vol(i,j,k)+dt*A_sponge_sound*sigma_sponge**2*(GrhoV(i,j,k)-rhoVsponge)
              GrhoW(i,j,k)=F_GrhoW(i,j,k)/vol(i,j,k)+dt*A_sponge_sound*sigma_sponge**2*(GrhoW(i,j,k)-0.0_WP)

	      endif	
              if(ym(j).ge.yloc2)then

              sigma_sponge=(ym(j)-yloc2)/sponge_length_y
              ! Create sponge target
              !pprime=pmax*sin(-3.1415d0/4.0d0+wavenum*(ym(j)-lambda/2.0d0)+2.0_WP*3.1415_WP*freq*time)
              pprime=pmax*sin(wavenum*(ym(j)-lambda/2.0d0)+2.0_WP*3.1415_WP*freq*time)
              rhoprime=pprime/sound_speed(i,j,k)**2
              vprime=pprime/(Grho(i,j,k)*sound_speed(i,j,k))

              rhosponge=1.226+rhoprime
              vbase=0.0!GrhoU(i,j,k)/Grho(i,j,k)
              rhoVsponge=-rhosponge*(vbase+vprime)
              pressurebase=101325.0_WP
              pressure=pressurebase+pprime
              rhoEsponge=pressure/0.4_WP+0.5_WP*rhoVsponge**2/rhosponge
              rhoUsponge=0.0_WP

              VOF  (i,j,k)=0.0_WP
              Grho (i,j,k)=F_Grho (i,j,k)/vol(i,j,k)+dt*A_sponge_sound*sigma_sponge**2*(Grho(i,j,k)-rhosponge)
              GrhoE(i,j,k)=F_GrhoE(i,j,k)/vol(i,j,k)+dt*A_sponge_sound*sigma_sponge**2*(GrhoE(i,j,k)-rhoEsponge)
              GrhoU(i,j,k)=F_GrhoU(i,j,k)/vol(i,j,k)+dt*A_sponge_sound*sigma_sponge**2*(GrhoU(i,j,k)-rhoUsponge)
              GrhoV(i,j,k)=F_GrhoV(i,j,k)/vol(i,j,k)+dt*A_sponge_sound*sigma_sponge**2*(GrhoV(i,j,k)-rhoVsponge)
              GrhoW(i,j,k)=F_GrhoW(i,j,k)/vol(i,j,k)+dt*A_sponge_sound*sigma_sponge**2*(GrhoW(i,j,k)-0.0_WP)

	      endif

                Lrho (i,j,k)=0.0_WP
                LrhoE(i,j,k)=0.0_WP
                LrhoU(i,j,k)=0.0_WP
                LrhoV(i,j,k)=0.0_WP
                LrhoW(i,j,k)=0.0_WP	      
                VOF  (i,j,k)=0.0_WP
               
            

             else if (VOF(i,j,k).gt.VOFhi) then
	       if(xm(i).gt.xloc1.and.xm(i).lt.xloc2)then	
	  	 A_sponge=0.0
		 sigma_sponge=0.0
	      Lrho (i,j,k)=F_Lrho (i,j,k)/vol(i,j,k)+dt*A_sponge*sigma_sponge**2*(Lrho(i,j,k)-Lrhosponge)
              LrhoE(i,j,k)=F_LrhoE(i,j,k)/vol(i,j,k)+dt*A_sponge*sigma_sponge**2*(LrhoE(i,j,k)-LrhoEsponge)
              LrhoU(i,j,k)=F_LrhoU(i,j,k)/vol(i,j,k)+dt*A_sponge*sigma_sponge**2*(LrhoU(i,j,k)-LrhoUsponge)
              LrhoV(i,j,k)=F_LrhoV(i,j,k)/vol(i,j,k)+dt*A_sponge*sigma_sponge**2*(LrhoV(i,j,k)-LrhoVsponge)
              LrhoW(i,j,k)=F_LrhoW(i,j,k)/vol(i,j,k)

	       elseif(xm(i).ge.xloc2)then
		 A_sponge=A_sponge_l_out
		 sigma_sponge=(xm(i)-xloc2)/(0.01_WP)
	      Lrho (i,j,k)=F_Lrho (i,j,k)/vol(i,j,k)!+dt*A_sponge*sigma_sponge**2*(Lrho(i,j,k)-Lrhosponge)
              LrhoE(i,j,k)=F_LrhoE(i,j,k)/vol(i,j,k)!+dt*A_sponge*sigma_sponge**2*(LrhoE(i,j,k)-LrhoEsponge)
              LrhoU(i,j,k)=F_LrhoU(i,j,k)/vol(i,j,k)+dt*A_sponge*sigma_sponge**2*(LrhoU(i,j,k)-LrhoUsponge)
              LrhoV(i,j,k)=F_LrhoV(i,j,k)/vol(i,j,k)+dt*A_sponge*sigma_sponge**2*(LrhoV(i,j,k)-LrhoVsponge)
              LrhoW(i,j,k)=F_LrhoW(i,j,k)/vol(i,j,k)

	       elseif(xm(i).le.xloc1)then
		 A_sponge=A_sponge_l_in
                 sigma_sponge=(xloc1-xm(i))/(0.01_WP)
	      Lrho (i,j,k)=F_Lrho (i,j,k)/vol(i,j,k)+dt*A_sponge*sigma_sponge**2*(Lrho(i,j,k)-Lrhosponge)
              LrhoE(i,j,k)=F_LrhoE(i,j,k)/vol(i,j,k)+dt*A_sponge*sigma_sponge**2*(LrhoE(i,j,k)-LrhoEsponge)
              LrhoU(i,j,k)=F_LrhoU(i,j,k)/vol(i,j,k)+dt*A_sponge*sigma_sponge**2*(LrhoU(i,j,k)-LrhoUsponge)
              LrhoV(i,j,k)=F_LrhoV(i,j,k)/vol(i,j,k)+dt*A_sponge*sigma_sponge**2*(LrhoV(i,j,k)-LrhoVsponge)
              LrhoW(i,j,k)=F_LrhoW(i,j,k)/vol(i,j,k)+dt*A_sponge*sigma_sponge**2*(LrhoW(i,j,k)-0.0_WP)

	       endif		

              VOF  (i,j,k)=1.0_WP
              Grho (i,j,k)=0.0_WP
              GrhoE(i,j,k)=0.0_WP
              GrhoU(i,j,k)=0.0_WP
              GrhoV(i,j,k)=0.0_WP
              GrhoW(i,j,k)=0.0_WP
              
           else
	     if(xm(i).gt.xloc1.and.xm(i).lt.xloc2)then	
		A_sponge=0.0
		sigma_sponge=0.0
	      Grho (i,j,k)=F_Grho (i,j,k)/((1.0_WP-VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(Grho(i,j,k)-Grhosponge)
              GrhoE(i,j,k)=F_GrhoE(i,j,k)/((1.0_WP-VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(GrhoE(i,j,k)-GrhoEsponge)
              GrhoU(i,j,k)=F_GrhoU(i,j,k)/((1.0_WP-VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(GrhoU(i,j,k)-GrhoUsponge)
              GrhoV(i,j,k)=F_GrhoV(i,j,k)/((1.0_WP-VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(GrhoV(i,j,k)-GrhoVsponge)
              GrhoW(i,j,k)=F_GrhoW(i,j,k)/((1.0_WP-VOF(i,j,k))*vol(i,j,k))
              Lrho (i,j,k)=F_Lrho (i,j,k)/((       VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(Lrho(i,j,k)-Lrhosponge)
              LrhoE(i,j,k)=F_LrhoE(i,j,k)/((       VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(LrhoE(i,j,k)-LrhoEsponge)
              LrhoU(i,j,k)=F_LrhoU(i,j,k)/((       VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(LrhoU(i,j,k)-LrhoUsponge)
              LrhoV(i,j,k)=F_LrhoV(i,j,k)/((       VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(LrhoV(i,j,k)-LrhoVsponge)
              LrhoW(i,j,k)=F_LrhoW(i,j,k)/((       VOF(i,j,k))*vol(i,j,k))!+dt*A_sponge*sigma_sponge**2*(LrhoW(i,j,k)-LrhoWsponge)


	     elseif(xm(i).ge.xloc2)then
		A_sponge=A_sponge_g_out
		sigma_sponge=(xm(i)-xloc2)/(0.01_WP)
	      Grho (i,j,k)=F_Grho (i,j,k)/((1.0_WP-VOF(i,j,k))*vol(i,j,k))!+dt*A_sponge*sigma_sponge**2*(Grho(i,j,k)-Grhosponge)
              GrhoE(i,j,k)=F_GrhoE(i,j,k)/((1.0_WP-VOF(i,j,k))*vol(i,j,k))!+dt*A_sponge*sigma_sponge**2*(GrhoE(i,j,k)-GrhoEsponge)
              GrhoU(i,j,k)=F_GrhoU(i,j,k)/((1.0_WP-VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(GrhoU(i,j,k)-GrhoUsponge)
              GrhoV(i,j,k)=F_GrhoV(i,j,k)/((1.0_WP-VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(GrhoV(i,j,k)-GrhoVsponge)
              GrhoW(i,j,k)=F_GrhoW(i,j,k)/((1.0_WP-VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(GrhoW(i,j,k)-0.0_WP)

	     elseif(xm(i).le.xloc1)then
		A_sponge=A_sponge_g_in
                sigma_sponge=(xloc1-xm(i))/(0.01)
	      Grho (i,j,k)=F_Grho (i,j,k)/((1.0_WP-VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(Grho(i,j,k)-Grhosponge)
              GrhoE(i,j,k)=F_GrhoE(i,j,k)/((1.0_WP-VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(GrhoE(i,j,k)-GrhoEsponge)
              GrhoU(i,j,k)=F_GrhoU(i,j,k)/((1.0_WP-VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(GrhoU(i,j,k)-GrhoUsponge)
              GrhoV(i,j,k)=F_GrhoV(i,j,k)/((1.0_WP-VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(GrhoV(i,j,k)-GrhoVsponge)
              GrhoW(i,j,k)=F_GrhoW(i,j,k)/((1.0_WP-VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(GrhoW(i,j,k)-0.0_WP)

	     endif		

	      
	     if(xm(i).ge.xloc2)then
		A_sponge=A_sponge_l_out
		sigma_sponge=(xm(i)-xloc2)/(0.01_WP)
	      Lrho (i,j,k)=F_Lrho (i,j,k)/((       VOF(i,j,k))*vol(i,j,k))!+dt*A_sponge*sigma_sponge**2*(Lrho(i,j,k)-Lrhosponge)
              LrhoE(i,j,k)=F_LrhoE(i,j,k)/((       VOF(i,j,k))*vol(i,j,k))!+dt*A_sponge*sigma_sponge**2*(LrhoE(i,j,k)-LrhoEsponge)
              LrhoU(i,j,k)=F_LrhoU(i,j,k)/((       VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(LrhoU(i,j,k)-LrhoUsponge)
              LrhoV(i,j,k)=F_LrhoV(i,j,k)/((       VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(LrhoV(i,j,k)-LrhoVsponge)
              LrhoW(i,j,k)=F_LrhoW(i,j,k)/((       VOF(i,j,k))*vol(i,j,k))!+dt*A_sponge*sigma_sponge**2*(LrhoW(i,j,k)-LrhoWsponge)

	     elseif(xm(i).le.xloc1)then
		A_sponge=A_sponge_l_in
                sigma_sponge=(xloc1-xm(i))/(0.01)
	      Lrho (i,j,k)=F_Lrho (i,j,k)/((       VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(Lrho(i,j,k)-Lrhosponge)
              LrhoE(i,j,k)=F_LrhoE(i,j,k)/((       VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(LrhoE(i,j,k)-LrhoEsponge)
              LrhoU(i,j,k)=F_LrhoU(i,j,k)/((       VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(LrhoU(i,j,k)-LrhoUsponge)
              LrhoV(i,j,k)=F_LrhoV(i,j,k)/((       VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(LrhoV(i,j,k)-LrhoVsponge)
              LrhoW(i,j,k)=F_LrhoW(i,j,k)/((       VOF(i,j,k))*vol(i,j,k))+dt*A_sponge*sigma_sponge**2*(LrhoW(i,j,k)-0.0_WP)

	     endif		

              
           end if

        end do
     end do
  end do
  
  ! Update boundaries
  call compressible_boundary_update(PA)
  call compressible_boundary_update(VOF)
  call compressible_boundary_update(Grho)
  call compressible_boundary_update(GrhoE)
  call compressible_boundary_update(GrhoU)
  call compressible_boundary_update(GrhoV)
  call compressible_boundary_update(GrhoW)
  call compressible_boundary_update(Lrho)
  call compressible_boundary_update(LrhoE)
  call compressible_boundary_update(LrhoU)
  call compressible_boundary_update(LrhoV)
  call compressible_boundary_update(LrhoW)

   do k=kmin_,kmax_
      do j=jmin_,jmax_
         do i=imin_,imax_ 

   	  if(xm(i).ge.xloc2)then
                sigma_sponge=(xm(i)-xloc2)/0.01_WP
                VOF(i,j,k)=0.0d0!VOF(i,j,k)-1000.0*sigma_sponge**2*(VOF(i,j,k)-0.0_WP)
           endif

           if(ym(j).le.yloc1)then
                sigma_sponge=(-ym(j)+yloc1)/sponge_length_y
                VOF(i,j,k)=0.0d0!VOF(i,j,k)-10000.0*sigma_sponge**2*(VOF(i,j,k)-0.0_WP)
           endif

           if(ym(j).ge.yloc2)then
                sigma_sponge=(ym(j)-yloc2)/sponge_length_y
                VOF(i,j,k)=0.0d0!VOF(i,j,k)-10000.0*sigma_sponge**2*(VOF(i,j,k)-0.0_WP)
           endif
	enddo
      enddo
   enddo	
  
  ! Update mixture variables
  RHO =VOF*Lrho +(1.0_WP-VOF)*Grho
  rhoE=VOF*LrhoE+(1.0_WP-VOF)*GrhoE
  rhoU=VOF*LrhoU+(1.0_WP-VOF)*GrhoU
  rhoV=VOF*LrhoV+(1.0_WP-VOF)*GrhoV
  rhoW=VOF*LrhoW+(1.0_WP-VOF)*GrhoW
  U=rhoU/RHO; V=rhoV/RHO; W=rhoW/RHO

   !print*, RHO(imax,5,5)
   !stop
  
  
  ! Perform PLIC reconstruction
  call compressible_plic_calc
  
  return
end subroutine compressible_advect_step


subroutine SpongeData(j,k,Grhosponge,GrhoUsponge,GrhoVsponge,GrhoEsponge,Lrhosponge,LrhoUsponge,LrhoVsponge,LrhoEsponge)

  use compressible
  implicit none
  integer :: j,k
  real(WP), parameter :: Ggamm=1.4_WP
  real(WP), parameter :: Lgamm=4.4_WP,LPref=600000000.0_WP

  real(WP) :: pressure,ulmax,ugmax,yminpert,ymaxpert,deltaypert
  real(WP) :: Grhosponge,GrhoUsponge,GrhoVsponge,GrhoEsponge,Lrhosponge,LrhoUsponge,LrhoVsponge,LrhoEsponge
   
   pressure = 101325.0_WP

   ! Re=2500, We=9
   ulmax=1.1125_WP
   ugmax=16.57_WP

   ! Re=2800, We=80
   !ulmax=1.246_WP
   !ugmax=40.34_WP


   ymaxpert=0.001_WP

   deltaypert=ymaxpert-yminpert

   !if(abs(ym(j)).le.0.002.and.abs(zm(k)).le.0.002)then
   if(dsqrt(ym(j)**2+zm(k)**2).le.0.001)then
   !VOF(i,j,k)=1.0
   Lrhosponge=1000.0_WP
   LrhoUsponge=1000.0_WP*ulmax!*(1.0_WP-(ym(j)/(yminpert+deltaypert))**2)
   LrhoVsponge=0.0_WP
   LrhoEsponge=(pressure+LGamm*LPref)/(Lgamm-1.0)+0.5_WP*LrhoUsponge**2/Lrhosponge
   GrhoUsponge=0.0_WP
   GrhoVsponge=0.0_WP
   endif
   if(dsqrt(ym(j)**2+zm(k)**2).ge.0.001.and.dsqrt(ym(j)**2+zm(k)**2).le.0.002d0)then
   !VOF(i,j,k)=1.0
   Lrhosponge=1000.0_WP
   LrhoUsponge=1000.0_WP*0.0d0!*(1.0_WP-(ym(j)/(yminpert+deltaypert))**2)
   LrhoVsponge=0.0_WP
   LrhoEsponge=(pressure+LGamm*LPref)/(Lgamm-1.0)+0.5_WP*LrhoUsponge**2/Lrhosponge
   Grhosponge=1.226_WP
   GrhoUsponge=1.226_WP*ugmax
   GrhoVsponge=0.0_WP
   GrhoEsponge=pressure/(Ggamm-1.0)+0.5_WP*GrhoUsponge**2/Grhosponge
   endif
   if(dsqrt(ym(j)**2+zm(k)**2).ge.0.002d0)then
   Grhosponge=1.226_WP
   GrhoUsponge=Grhosponge*0.0_WP
   GrhoVsponge=0.0_WP
   GrhoEsponge=pressure/(Ggamm-1.0_WP)+0.5_WP*GrhoUsponge**2/Grhosponge
   LrhoUsponge=0.0
   LrhoVsponge=0.0
   endif

end subroutine SpongeData


! ==================================================== !
! Form conservative linear reconstruction in each cell !
! ==================================================== !
subroutine compressible_advect_reconstruct
  use compressible_advect
  implicit none
  integer :: i,j,k,ip,jp,kp,im,jm,km
  integer :: n,ntet,dir,case,v1,v2
  real(WP) :: idp,idm,mu,my_vol
  real(WP), dimension(3,8) :: vert,pt
  real(WP), dimension(4) :: d
  real(WP), dimension(3) :: a,b,c
  
  ! Zero gradients
  gradGrho =0.0_WP
  gradGrhoE=0.0_WP
  gradGrhoU=0.0_WP
  gradGrhoV=0.0_WP
  gradGrhoW=0.0_WP
  gradLrho =0.0_WP
  gradLrhoE=0.0_WP
  gradLrhoU=0.0_WP
  gradLrhoV=0.0_WP
  gradLrhoW=0.0_WP
  
  ! Form barycenter for each phase
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           ! Choose treatment based on VOF value
           if (VOF(i,j,k).lt.VOFlo.or.VOF(i,j,k).gt.VOFhi) then
              Gbary(:,i,j,k)=[xm(i),ym(j),zm(k)]
              Lbary(:,i,j,k)=[xm(i),ym(j),zm(k)]
           else
              ! Zero barycenters
              Gbary(:,i,j,k)=0.0_WP
              Lbary(:,i,j,k)=0.0_WP
              ! Form cell vertices
              pt(:,1)=[x(i  ),y(j  ),z(k  )]
              pt(:,2)=[x(i+1),y(j  ),z(k  )]
              pt(:,3)=[x(i  ),y(j+1),z(k  )]
              pt(:,4)=[x(i+1),y(j+1),z(k  )]
              pt(:,5)=[x(i  ),y(j  ),z(k+1)]
              pt(:,6)=[x(i+1),y(j  ),z(k+1)]
              pt(:,7)=[x(i  ),y(j+1),z(k+1)]
              pt(:,8)=[x(i+1),y(j+1),z(k+1)]
              ! Cut each sub-tet
              do ntet=1,5
                 ! Create tet and distance
                 do n=1,4
                    vert(:,n)=pt(:,verts2tets(n,ntet))
                    d(n)=xnorm(i,j,k)*vert(1,n)+ynorm(i,j,k)*vert(2,n)+znorm(i,j,k)*vert(3,n)-dist(i,j,k)
                 end do
                 ! Find cut case
                 case=1+int(0.5_WP+sign(0.5_WP,d(1)))+&
                      2*int(0.5_WP+sign(0.5_WP,d(2)))+&
                      4*int(0.5_WP+sign(0.5_WP,d(3)))+&
                      8*int(0.5_WP+sign(0.5_WP,d(4)))
                 ! Create interpolated vertices on cut plane
                 do n=1,cut_nvert(case)
                    v1=cut_v1(n,case); v2=cut_v2(n,case)
                    mu=min(1.0_WP,max(0.0_WP,-d(v1)/(sign(abs(d(v2)-d(v1))+epsilon(1.0_WP),d(v2)-d(v1)))))
                    vert(:,4+n)=(1.0_WP-mu)*vert(:,v1)+mu*vert(:,v2)
                 end do
                 ! Analyze liquid tets
                 do n=cut_ntets(case),cut_nntet(case),-1
                    ! Compute volume
                    a=vert(:,cut_vtet(1,n,case))-vert(:,cut_vtet(4,n,case)); b=vert(:,cut_vtet(2,n,case))-vert(:,cut_vtet(4,n,case)); c=vert(:,cut_vtet(3,n,case))-vert(:,cut_vtet(4,n,case))
                    my_vol=abs(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
                    ! Compute barycenter
                    Lbary(:,i,j,k)=Lbary(:,i,j,k)+0.25_WP*my_vol*(vert(:,cut_vtet(1,n,case))+vert(:,cut_vtet(2,n,case))+vert(:,cut_vtet(3,n,case))+vert(:,cut_vtet(4,n,case)))
                 end do
                 ! Analyze gas tets
                 do n=1,cut_nntet(case)-1
                    ! Compute volume
                    a=vert(:,cut_vtet(1,n,case))-vert(:,cut_vtet(4,n,case)); b=vert(:,cut_vtet(2,n,case))-vert(:,cut_vtet(4,n,case)); c=vert(:,cut_vtet(3,n,case))-vert(:,cut_vtet(4,n,case))
                    my_vol=abs(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
                    ! Compute barycenter
                    Gbary(:,i,j,k)=Gbary(:,i,j,k)+0.25_WP*my_vol*(vert(:,cut_vtet(1,n,case))+vert(:,cut_vtet(2,n,case))+vert(:,cut_vtet(3,n,case))+vert(:,cut_vtet(4,n,case)))
                 end do
              end do
              ! Normalize barycenters
              Gbary(:,i,j,k)=Gbary(:,i,j,k)/((1.0_WP-VOF(i,j,k))*vol(i,j,k))
              Lbary(:,i,j,k)=Lbary(:,i,j,k)/((       VOF(i,j,k))*vol(i,j,k))
           end if
        end do
     end do
  end do
  
  ! Gradients for linear reconstruction
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           do dir=1,3
              select case (dir)
              case (1) ! X gradient
                 ip=i+1; jp=j; kp=k; idp=1.0_WP/(xm(i+1)-xm(i))
                 im=i-1; jm=j; km=k ;idm=1.0_WP/(xm(i)-xm(i-1))
              case (2) ! Y gradient
                 ip=i; jp=j+1; kp=k; idp=1.0_WP/(ym(j+1)-ym(j))
                 im=i; jm=j-1; km=k ;idm=1.0_WP/(ym(j)-ym(j-1))
              case (3) ! Z gradient
                 ip=i; jp=j; kp=k+1; idp=1.0_WP/(zm(k+1)-zm(k))
                 im=i; jm=j; km=k-1 ;idm=1.0_WP/(zm(k)-zm(k-1))
              end select
              gradGrho (dir,i,j,k)=mmgrad((Grho (ip,jp,kp)-Grho (i,j,k))*idp,(Grho (i,j,k)-Grho (im,jm,km))*idm)
              gradGrhoE(dir,i,j,k)=mmgrad((GrhoE(ip,jp,kp)-GrhoE(i,j,k))*idp,(GrhoE(i,j,k)-GrhoE(im,jm,km))*idm)
              gradGrhoU(dir,i,j,k)=mmgrad((GrhoU(ip,jp,kp)-GrhoU(i,j,k))*idp,(GrhoU(i,j,k)-GrhoU(im,jm,km))*idm)
              gradGrhoV(dir,i,j,k)=mmgrad((GrhoV(ip,jp,kp)-GrhoV(i,j,k))*idp,(GrhoV(i,j,k)-GrhoV(im,jm,km))*idm)
              gradGrhoW(dir,i,j,k)=mmgrad((GrhoW(ip,jp,kp)-GrhoW(i,j,k))*idp,(GrhoW(i,j,k)-GrhoW(im,jm,km))*idm)
              gradLrho (dir,i,j,k)=mmgrad((Lrho (ip,jp,kp)-Lrho (i,j,k))*idp,(Lrho (i,j,k)-Lrho (im,jm,km))*idm)
              gradLrhoE(dir,i,j,k)=mmgrad((LrhoE(ip,jp,kp)-LrhoE(i,j,k))*idp,(LrhoE(i,j,k)-LrhoE(im,jm,km))*idm)
              gradLrhoU(dir,i,j,k)=mmgrad((LrhoU(ip,jp,kp)-LrhoU(i,j,k))*idp,(LrhoU(i,j,k)-LrhoU(im,jm,km))*idm)
              gradLrhoV(dir,i,j,k)=mmgrad((LrhoV(ip,jp,kp)-LrhoV(i,j,k))*idp,(LrhoV(i,j,k)-LrhoV(im,jm,km))*idm)
              gradLrhoW(dir,i,j,k)=mmgrad((LrhoW(ip,jp,kp)-LrhoW(i,j,k))*idp,(LrhoW(i,j,k)-LrhoW(im,jm,km))*idm)
           end do
        end do
     end do
  end do
  
  ! Communication
  do dir=1,3
     call compressible_boundary_update(gradGrho (dir,:,:,:))
     call compressible_boundary_update(gradGrhoE(dir,:,:,:))
     call compressible_boundary_update(gradGrhoU(dir,:,:,:))
     call compressible_boundary_update(gradGrhoV(dir,:,:,:))
     call compressible_boundary_update(gradGrhoW(dir,:,:,:))
     call compressible_boundary_update(gradLrho (dir,:,:,:))
     call compressible_boundary_update(gradLrhoE(dir,:,:,:))
     call compressible_boundary_update(gradLrhoU(dir,:,:,:))
     call compressible_boundary_update(gradLrhoV(dir,:,:,:))
     call compressible_boundary_update(gradLrhoW(dir,:,:,:))
  end do
  
contains
  
  ! Minmod gradient
  function mmgrad(g1,g2) result(g)
    implicit none
    real(WP), intent(in) :: g1,g2
    real(WP) :: g
    if (g1*g2.le.0.0_WP) then
       g=0.0_WP
    else
       if (abs(g1).lt.abs(g2)) then
          g=g1
       else
          g=g2
       end if
    end if
  end function mmgrad
  
end subroutine compressible_advect_reconstruct
