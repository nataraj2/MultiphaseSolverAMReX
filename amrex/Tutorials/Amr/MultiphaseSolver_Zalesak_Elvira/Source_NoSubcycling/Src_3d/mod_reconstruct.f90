module reconstruct_module

use amrex_base_module
use amrex_amr_module
use amr_data_module

contains

subroutine compute_reconstruct(ilev,phiborder)
    implicit none
    integer :: ilev
    integer :: clo(4), chi(4), borlo(4), borhi(4)
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: phiborder_ptr, gradxphi_ptr,gradyphi_ptr,gradzphi_ptr
    type(amrex_multifab), intent (in) :: phiborder

    !do ilev = 0, amrex_get_finest_level()
       !$omp parallel private(rlo,rhi,elo,ehi,bx,mfi,prhs,pexact)
       call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.true.)

       do while (mfi%next())
          bx = mfi%tilebox()
          phiborder_ptr  =>  phiborder%dataptr(mfi)
          gradxphi_ptr   =>  gradxphi(ilev) % dataptr(mfi)
          gradyphi_ptr   =>  gradyphi(ilev) % dataptr(mfi)
          gradzphi_ptr   =>  gradzphi(ilev) % dataptr(mfi)
	  clo = lbound(gradxphi_ptr)
          chi = ubound(gradxphi_ptr)
          borlo = lbound(phiborder_ptr)
          borhi = ubound(phiborder_ptr)
          call reconstruct(clo(1:3), chi(1:3), borlo, borhi, phiborder_ptr,  gradxphi_ptr, gradyphi_ptr, gradzphi_ptr, amrex_problo, amrex_geom(ilev)%dx)
       end do

       call amrex_mfiter_destroy(mfi)
       !$omp end parallel

    !end do

  end subroutine compute_reconstruct

  subroutine reconstruct(clo, chi, borlo, borhi, phiborder, gradxphi_mat, gradyphi_mat, gradzphi_mat, prob_lo, dx)
  implicit none

  integer, intent(in) :: clo(3), chi(3)
  integer, intent(in) :: borlo(3), borhi(3)
  double precision, intent(in) :: phiborder(borlo(1):borhi(1),borlo(2):borhi(2),borlo(3):borhi(3),ncomp)
  double precision, intent(out) :: gradxphi_mat(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3),ncomp)
  double precision, intent(out) :: gradyphi_mat(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3),ncomp)
  double precision, intent(out) :: gradzphi_mat(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3),ncomp)
  double precision, intent(in) :: prob_lo(3), dx(3)

  integer :: i,j,k,ip,jp,kp,im,jm,km
  integer :: n,ntet,dir,case,v1,v2
  double precision :: idp,idm,mu,my_vol
  double precision, dimension(3,8) :: vert,pt
  double precision, dimension(4) :: d
  double precision, dimension(3) :: a,b,c
  integer :: icomp
  double precision :: xm1,xm2,xm3,ym1,ym2,ym3,zm1,zm2,zm3

  double precision, dimension (borlo(1):borhi(1)) :: xm
  double precision, dimension (borlo(2):borhi(2)) :: ym
  double precision, dimension (borlo(3):borhi(3)) :: zm

  do i=borlo(1),borhi(1)
     xm(i) = prob_lo(1)+(dble(i)+0.5d0)*dx(1)
  end do   
  do j=borlo(2),borhi(2)
     ym(j) = prob_lo(2)+(dble(j)+0.5d0)*dx(2)
  end do   
  do k=borlo(3),borhi(3)
     zm(k) = prob_lo(3)+(dble(k)+0.5d0)*dx(3)
  end do   


! Gradients for linear reconstruction
  do k=clo(3),chi(3)
     do j=clo(2),chi(2)
        do i=clo(1),chi(1)
           do dir=1,3
              select case (dir)
              case (1) ! X gradient
		 ip=i+1; jp=j; kp=k; idp=1.0d0/(xm(i+1)-xm(i))
                 im=i-1; jm=j; km=k ;idm=1.0d0/(xm(i)-xm(i-1))
		 do icomp=1, ncomp
                    gradxphi_mat(i,j,k,icomp)=mmgrad((phiborder(ip,jp,kp,icomp)-phiborder(i,j,k,icomp))*idp,(phiborder(i,j,k,icomp)-phiborder(im,jm,km,icomp))*idm)
		 end do
              case (2) ! Y gradient
	 	 ip=i; jp=j+1; kp=k; idp=1.0d0/(ym(j+1)-ym(j))
                 im=i; jm=j-1; km=k ;idm=1.0d0/(ym(j)-ym(j-1))
		 do icomp=1, ncomp
                    gradyphi_mat(i,j,k,icomp)=mmgrad((phiborder(ip,jp,kp,icomp)-phiborder(i,j,k,icomp))*idp,(phiborder(i,j,k,icomp)-phiborder(im,jm,km,icomp))*idm)
		 end do
              case (3) ! Z gradient
		 ip=i; jp=j; kp=k+1; idp=1.0d0/(zm(k+1)-zm(k))
                 im=i; jm=j; km=k-1 ;idm=1.0d0/(zm(k)-zm(k-1))
		 do icomp=1, ncomp
                    gradzphi_mat(i,j,k,icomp)=mmgrad((phiborder(ip,jp,kp,icomp)-phiborder(i,j,k,icomp))*idp,(phiborder(i,j,k,icomp)-phiborder(im,jm,km,icomp))*idm)
		 end do

              end select

           end do
        end do
     end do
  end do

contains

  ! Minmod gradient
  function mmgrad(g1,g2) result(g)
    implicit none
    double precision, intent(in) :: g1,g2
    double precision :: g
    if (g1*g2.le.0.0d0) then
       g=0.0d0
    else
       if (abs(g1).lt.abs(g2)) then
          g=g1
       else
          g=g2
       end if
    end if
  end function mmgrad

end subroutine reconstruct

end module reconstruct_module
