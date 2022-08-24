module mod_compute_plic_half
  use amrex_base_module
  use my_amr_module
  use amr_data_module
  use precision 
  implicit none
  
contains

subroutine build_vol_multifab

    use amr_data_module
    integer :: lev, nlevs

    nlevs = amrex_get_finest_level()

    do lev=0, nlevs
       call amrex_multifab_build(Gvolmfab(lev), phi_new(lev)%ba, phi_new(lev)%dm, 6, 1)
       call amrex_multifab_build(Lvolmfab(lev), phi_new(lev)%ba, phi_new(lev)%dm, 6, 1)
    enddo

end subroutine build_vol_multifab

subroutine compute_plic_half
    use fillpatch_module, only : fillpatch!, fillpatchthree
    implicit none
    integer :: ilev
    integer :: clo(4),chi(4),borlo(4),borhi(4),vollo(4),volhi(4)
    integer :: ngrow=3
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    type(amrex_multifab) :: phiborder, null(0:amrex_max_level)
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: phi, phiborderptr,Gvolptr,Lvolptr
    integer :: step

    do ilev = 0, amrex_get_finest_level()

      ! Create and fill solution vector multifab
      call amrex_multifab_build(phiborder, phi_new(ilev)%ba, phi_new(ilev)%dm, ncomp, ngrow)
      call fillpatch(ilev,t_new(0),phiborder,phi_new,phi_old) 

      call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.true.)

       do while (mfi%next())
          bx = mfi%tilebox()
	  phi => phi_new(ilev)%dataptr(mfi)
          phiborderptr => phiborder%dataptr(mfi)
	  Gvolptr => Gvolmfab(ilev)%dataptr(mfi)
	  Lvolptr => Lvolmfab(ilev)%dataptr(mfi)
	  clo = lbound(phi)
	  chi = ubound(phi)
          borlo = lbound(phiborderptr)
          borhi = ubound(phiborderptr)
	  vollo = lbound(Gvolptr)
	  volhi = ubound(Lvolptr)
	  
          call plic_half(bx%lo,bx%hi,clo, chi, borlo, borhi, phiborderptr(:,:,:,ncomp), Gvolptr, Lvolptr, vollo, volhi, amrex_problo, amrex_geom(ilev)%dx)

       end do

       call amrex_mfiter_destroy(mfi)
       call amrex_multifab_destroy(phiborder)
    end do

  end subroutine compute_plic_half
! ====================================== !
! Creates unstructured mesh on interface !
! ====================================== !

  subroutine plic_half(lo, hi, clo, chi, borlo, borhi, VOF, Gvol, Lvol, vollo, volhi, prob_lo, dx)
  use amr_data_module
  use grid_module
  use irl_fortran_interface 
  implicit none

  !!!!!!!!!!!!!!!!!!!!!
  integer, intent(in) :: lo(3), hi(3), clo(4), chi(4), borlo(4), borhi(4), vollo(4), volhi(4)
  double precision, intent(in), dimension (borlo(1):borhi(1),borlo(2):borhi(2),borlo(3):borhi(3)) :: VOF!, xnorm, ynorm, znorm
  double precision, intent(out), dimension (vollo(1):volhi(1),vollo(2):volhi(2),vollo(3):volhi(3),6) :: Gvol, Lvol
  double precision, dimension (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) ::  xnorm, ynorm, znorm, dist
  double precision :: x(lo(1)-1:hi(1)+2), y(lo(2)-1:hi(2)+2), z(lo(3)-1:hi(3)+2) 
  double precision :: prob_lo(3), dx(3)
  !!!!!!!!!!!!!!!!!!!!!

  integer :: i,j,k,ii,jj,kk
  integer :: ilev
  type(PlanarSep_type) :: pseparator(1,1,1)
  double precision, dimension(borlo(1):borhi(1),borlo(2):borhi(2),borlo(3):borhi(3)) :: tmp
  real(8), dimension(4) :: plane
  double precision :: normval
  type(RectCub_type) :: cell
  integer(IRL_LargeOffsetIndex_t) :: heapsize
  type(ObjServer_PlanarSep_type) :: planar_separator_allocation
  type(SepVM_type) :: GandLVols
  double precision :: xm, ym, zm, cell_vol

   call new(cell)
   call new(GandLVols)

 ! Initialize size for IRL
   heapsize = int(1,8)
   call new(planar_separator_allocation,heapsize) 
   call new(pseparator(1,1,1),planar_separator_allocation)

    do i=lo(1)-1,hi(1)+2
       x(i) = dble(i)*dx(1) + prob_lo(1)
    end do
    do j=lo(2)-1,hi(2)+2
       y(j) = dble(j)*dx(2) + prob_lo(2)
    end do
    do k =lo(3)-1,hi(3)+2
       z(k) = dble(k)*dx(3) + prob_lo(3)
    end do

   xnorm=0.0d0;ynorm=0.0d0;znorm=0.0d0;dist=0.0d0

  ! First smooth VOF field
  call smoother(lo,hi,tmp,VOF,borlo,borhi)
 
  ! Normal update using Young's method
  do k=lo(3)-1,hi(3)+1
     do j=lo(2)-1,hi(2)+1
        do i=lo(1)-1,hi(1)+1
           ! VOF gradient
           xnorm(i,j,k)=(tmp(i+1,j,k)-tmp(i-1,j,k))/(2.0d0*dx(1))
           ynorm(i,j,k)=(tmp(i,j+1,k)-tmp(i,j-1,k))/(2.0d0*dx(2))
           znorm(i,j,k)=(tmp(i,j,k+1)-tmp(i,j,k-1))/(2.0d0*dx(3))
           ! Normalize it
           normval=sqrt(xnorm(i,j,k)**2+ynorm(i,j,k)**2+znorm(i,j,k)**2)
           if (normval.gt.0.0d0) then
              xnorm(i,j,k)=-xnorm(i,j,k)/normval
              ynorm(i,j,k)=-ynorm(i,j,k)/normval
              znorm(i,j,k)=-znorm(i,j,k)/normval
           else
              xnorm(i,j,k)=0.0d0
              ynorm(i,j,k)=0.0d0
              znorm(i,j,k)=0.0d0
           end if
	   call setNumberOfPlanes(pseparator(1,1,1),1)
           call setPlane(pseparator(1,1,1),0,(/xnorm(i,j,k),ynorm(i,j,k),znorm(i,j,k)/),0.0d0)
           call construct_2pt(cell,(/x(i),y(j),z(k)/),(/x(i+1),y(j+1),z(k+1)/) )
           call matchVolumeFraction(cell,VOF(i,j,k),pseparator(1,1,1))
	   
	   xm=(x(i)+x(i+1))/2.0d0
	   ym=(y(j)+y(j+1))/2.0d0
	   zm=(z(k)+z(k+1))/2.0d0

	   cell_vol=dx(1)*dx(2)*dx(3)
	   if(VOF(i,j,k).gt.VOFhi)then
	      Lvol(i,j,k,:)=0.5_WP*cell_vol
              Gvol(i,j,k,:)=0.0_WP		
	   else if(VOF(i,j,k).lt.VOFlo)then
	      Lvol(i,j,k,:)=0.0_WP
              Gvol(i,j,k,:)=0.5_WP*cell_vol
	   else if(VOF(i,j,k).gt.VOFlo.and.VOF(i,j,k).lt.VOFhi)then

	   call construct_2pt(cell,(/x(i),y(j),z(k)/),(/xm,y(j+1),z(k+1)/) )
           call getNormMoments(cell,pseparator(1,1,1),GandLVols)
           Lvol(i,j,k,1)=getVolumePtr(GandLVols,0)
           Gvol(i,j,k,1)=getVolumePtr(GandLVols,1)

           call construct_2pt(cell,(/xm,y(j),z(k)/),(/x(i+1),y(j+1),z(k+1)/) )
           call getNormMoments(cell,pseparator(1,1,1),GandLVols)
           Lvol(i,j,k,2)=getVolumePtr(GandLVols,0)
           Gvol(i,j,k,2)=getVolumePtr(GandLVols,1)

           call construct_2pt(cell,(/x(i),y(j),z(k)/),(/x(i+1),ym,z(k+1)/) )
           call getNormMoments(cell,pseparator(1,1,1),GandLVols)
           Lvol(i,j,k,3)=getVolumePtr(GandLVols,0)
           Gvol(i,j,k,3)=getVolumePtr(GandLVols,1)

           call construct_2pt(cell,(/x(i),ym,z(k)/),(/x(i+1),y(j+1),z(k+1)/) )
           call getNormMoments(cell,pseparator(1,1,1),GandLVols)
           Lvol(i,j,k,4)=getVolumePtr(GandLVols,0)
           Gvol(i,j,k,4)=getVolumePtr(GandLVols,1)

           call construct_2pt(cell,(/x(i),y(j),z(k)/),(/x(i+1),y(j+1),zm/) )
           call getNormMoments(cell,pseparator(1,1,1),GandLVols)
           Lvol(i,j,k,5)=getVolumePtr(GandLVols,0)
           Gvol(i,j,k,5)=getVolumePtr(GandLVols,1)

           call construct_2pt(cell,(/x(i),y(j),zm/),(/x(i+1),y(j+1),z(k+1)/) )
           call getNormMoments(cell,pseparator(1,1,1),GandLVols)
           Lvol(i,j,k,6)=getVolumePtr(GandLVols,0)
           Gvol(i,j,k,6)=getVolumePtr(GandLVols,1)
           endif

       end do
     end do
  end do

contains

  ! =========================================== !
  ! Smoothing of VOF to improve Young's normals !
  ! =========================================== !
  subroutine smoother(lo,hi,tmp,VOF,clo,chi)
    use precision
    implicit none
    integer :: i,j,k,ii,jj,kk,count
    integer, intent(in) :: lo(3), hi(3), clo(3), chi(3)
    double precision, dimension(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)), intent(in) :: VOF
    double precision, dimension(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)), intent(out) :: tmp

    ! Perform local average
    tmp=0.0d0
    do k=lo(3)-2,hi(3)+2
       do j=lo(2)-2,hi(2)+2
          do i=lo(1)-2,hi(1)+2
          ! Smooth local value
             count=0
            do kk=k-1,k+1
               do jj=j-1,j+1
                  do ii=i-1,i+1
                     ! Average
                      tmp(i,j,k)=tmp(i,j,k)+VOF(ii,jj,kk)
                      count=count+1
                   end do
                end do
            end do
            tmp(i,j,k)=tmp(i,j,k)/real(count,WP)
           end do
        end do
     end do

  end subroutine smoother

end subroutine plic_half

end module mod_compute_plic_half
