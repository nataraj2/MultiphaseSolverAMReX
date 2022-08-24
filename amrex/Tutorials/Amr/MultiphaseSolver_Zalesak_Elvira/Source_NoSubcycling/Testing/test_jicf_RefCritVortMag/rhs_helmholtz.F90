module rhs_helmholtz
  use amrex_base_module
  use amrex_amr_module
  use amr_data_module
  use fillpatch_module

  implicit none
  private
  public ::  init_prob_start, init_prob_data, compute_helmholtz_rhs

  real(amrex_real), private, parameter :: a = 1.0d0
  real(amrex_real), private, parameter :: b = 1.0d0

contains

  subroutine init_prob_start(ilev)
 
    integer, intent(in) :: ilev
    type(amrex_box) :: bx, gbx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: phi
    integer :: nlevs


       call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.true.)
       
       do while (mfi%next())
          bx  =  mfi%tilebox()
          gbx =  mfi%growntilebox(1)
          phi => phi_new(ilev)%dataptr(mfi)

	  call init_prob_data(ilev, t_new(ilev), bx%lo, bx%hi, phi(:,:,:,1),phi(:,:,:,2),phi(:,:,:,3),phi(:,:,:,4),phi(:,:,:,5),phi(:,:,:,6),&
							       phi(:,:,:,7),phi(:,:,:,8),phi(:,:,:,9),phi(:,:,:,10),phi(:,:,:,11),&
							       lbound(phi), ubound(phi), &
               amrex_geom(ilev)%dx, amrex_problo)

       end do

       call amrex_mfiter_destroy(mfi)

       !call Pmfab(ilev)%setVal(0.0_amrex_real)

  end subroutine init_prob_start

  subroutine init_prob_data(level, time, lo, hi, &
     Grho, GrhoU, GrhoV, GrhoW, GrhoE, &
     Lrho, LrhoU, LrhoV, LrhoW, LrhoE, &
     VOF, var_lo, var_hi,&
     dx, prob_lo) bind(C, name="initdata")

  implicit none
  integer, intent(in) :: level, lo(3), hi(3), var_lo(3), var_hi(3)
  double precision, intent(in) :: time
  double precision, intent(inout), dimension(var_lo(1):var_hi(1),var_lo(2):var_hi(2),var_lo(3):var_hi(3)) :: Grho,GrhoU,GrhoV,GrhoW,GrhoE,Lrho,LrhoU,LrhoV,LrhoW,LrhoE,VOF
  double precision, intent(in) :: dx(3), prob_lo(3)
  integer          :: i,j,k
  double precision :: x,y,z,r2
  double precision :: xcen(5),ycen(5), zcen(5)
  real(amrex_real) :: pi
  real(amrex_real) :: Tinfty, Tj, Vj, Mj, rhoj, Vz, afac, rad, jetrad

     Tinfty=300.0d0
     Tj=300.0d0*1.0d0
     Vj=400.0d0
     Mj=Vj/dsqrt(1.4d0*287.0d0*Tj)
     afac=1.0d0
     jetrad=0.5d0
     rhoj=101325.0d0/(287.0d0*Tj)


  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
        y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
        rad = dsqrt(y**2+(z-5.0d0)**2)
        Vz=0.5d0*(1.0d0+tanh(afac*(jetrad/rad-rad/jetrad)))

        do i=lo(1),hi(1)
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

                if(x.le.1.0d0)then
                Grho(i,j,k) = rhoj*(0.5d0*0.4d0*Vz*(1.0d0-Vz)*Mj**2+Vz+Tinfty/Tj*(1.0d0-Vz))**(-1)
                GrhoU(i,j,k) = Grho(i,j,k)*Vz*Vj
                else
		Grho(i,j,k)=1.176d0
                GrhoU(i,j,k) = 0.0d0
                endif
                GrhoV(i,j,k) = 0.0d0
                GrhoW(i,j,k) = 1.176*30.0d0
                GrhoE(i,j,k) = 101325.0/0.4d0 + 0.5d0*(GrhoU(i,j,k)**2+GrhoV(i,j,k)**2+GrhoW(i,j,k)**2)/Grho(i,j,k)

		Lrho(i,j,k)=0.0d0
		LrhoU(i,j,k)=0.0d0
		LrhoV(i,j,k)=0.0d0
		LrhoW(i,j,k)=0.0d0
		LrhoE(i,j,k)=0.0d0
		VOF(i,j,k)=0.0d0

        end do
     end do
  end do

  end subroutine init_prob_data

  subroutine compute_helmholtz_rhs
   
    implicit none

    integer :: ilev
    integer :: clo(4), chi(4), borlo(4), borhi(4), fxlo(4), fxhi(4), fylo(4), fyhi(4), fzlo(4), fzhi(4), betalo(4), betahi(4)
    type(amrex_box) :: bx, gbx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: phi,phi_nborptr,Ufaceptr,Vfaceptr,Wfaceptr,rhsptr,alphaptr,betaptr,betafabptr
    integer :: nlevs
    integer, parameter :: ngrow = 3
    type(amrex_multifab) :: phi_mp_nbormfab
    real(amrex_real) :: time
    integer :: nc,ng,idim
    logical :: nodal(3)

    time=0.0d0

    nlevs = amrex_get_finest_level()

    do ilev = 0, nlevs

      call amrex_multifab_build(Pmfab(ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, nc=1, ng=1)
      call amrex_multifab_build(exact_solution(ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, nc=1, ng=0)
      call amrex_multifab_build(rhs(ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, nc=1, ng=0)
      call amrex_multifab_build(acoef(ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, nc=1, ng=0)
      call amrex_multifab_build(bcoef(ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, nc=1, ng=1)

      do idim = 1, amrex_spacedim
          nodal = .false.
          nodal(idim) = .true.
          call amrex_multifab_build(beta(idim,ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, 1, 0, nodal)
          call amrex_multifab_build(gradPmfab(idim,ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, 1, 0, nodal)
       end do 
   enddo
 
    do ilev = 0, nlevs

       call amrex_multifab_build(phi_mp_nbormfab, phi_n(ilev)%ba, phi_n(ilev)%dm, ncomp, ngrow)
       call fillpatch (ilev, time, phi_mp_nbormfab, phi_mp_n, phi_mp_n)

       call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.true.)
       do while (mfi%next())
          bx = mfi%tilebox()
	  gbx = mfi%growntilebox(1)
          phi          => phi_mp(ilev)%dataptr(mfi)
          phi_nborptr  => phi_mp_nbormfab%dataptr(mfi)
	  alphaptr     => acoef(ilev)%dataptr(mfi)
	  Ufaceptr   => facevelmfab(1,ilev)%dataptr(mfi)
	  Vfaceptr   => facevelmfab(2,ilev)%dataptr(mfi)
	  Wfaceptr   => facevelmfab(3,ilev)%dataptr(mfi)
	  betaptr      => bcoef(ilev)%dataptr(mfi)
          rhsptr       => rhs(ilev) % dataptr(mfi)
          clo    = lbound(phi)
          chi    = ubound(phi)
	  borlo  = lbound(phi_nborptr)
          borhi  = ubound(phi_nborptr)
	  fxlo    = lbound(Ufaceptr)
          fxhi    = ubound(Ufaceptr)
	  fylo    = lbound(Vfaceptr)
          fyhi    = ubound(Vfaceptr)
	  fzlo    = lbound(Wfaceptr)
          fzhi    = ubound(Wfaceptr)
	  betalo = lbound(betaptr)
          betahi = ubound(betaptr)

          call helmholtz_rhs(ilev, bx%lo, bx%hi, gbx%lo, gbx%hi,&
			     phi(:,:,:,1),phi(:,:,:,2),phi(:,:,:,3),phi(:,:,:,4),phi(:,:,:,5),alphaptr,clo,chi,&
			     phi_nborptr(:,:,:,1),phi_nborptr(:,:,:,2),phi_nborptr(:,:,:,3),phi_nborptr(:,:,:,4),phi_nborptr(:,:,:,5),borlo,borhi,&
			     Ufaceptr,fxlo,fxhi,Vfaceptr,fylo,fyhi,Wfaceptr,fzlo,fzhi,betaptr,betalo,betahi,rhsptr,amrex_geom(ilev)%dx, amrex_problo)

       end do
       call amrex_mfiter_destroy(mfi)

       call amrex_multifab_destroy(phi_mp_nbormfab)

       call Pmfab(ilev)%setVal(0.0_amrex_real)

    end do
  

  end subroutine compute_helmholtz_rhs


  subroutine helmholtz_rhs(level,lo,hi,glo,ghi,RHO,rhoU,rhoV,rhoW,rhoE,alpha,clo,chi,&
			   RHO_n,rhoU_n,rhoV_n,rhoW_n,rhoE_n,borlo,borhi,&
			   Uface,fxlo,fxhi,Vface,fylo,fyhi,Wface,fzlo,fzhi,beta,betalo,betahi,rhs,dx,prob_lo)
 
  implicit none
  integer, intent(in) :: level, lo(3), hi(3), clo(4), chi(4), borlo(4), borhi(4), fxlo(4), fxhi(4), fylo(4), fyhi(4), fzlo(4), fzhi(4), betalo(4), betahi(4), glo(3), ghi(3)
  double precision, dimension(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)) :: RHO,rhoU,rhoV,rhoW,rhoE,alpha,rhs
  double precision, dimension(borlo(1):borhi(1),borlo(2):borhi(2),borlo(3):borhi(3)) :: RHO_n,rhoU_n,rhoV_n,rhoW_n,rhoE_n
  double precision, dimension(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3)) :: Uface
  double precision, dimension(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3)) :: Vface
  double precision, dimension(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3)) :: Wface
  double precision, dimension(betalo(1):betahi(1),betalo(2):betahi(2),betalo(3):betahi(3)) :: beta
  double precision, intent(in) :: dx(3), prob_lo(3)
  integer          :: i,j,k
  double precision :: x,y,z
  real(amrex_real) :: PA, P_n, spdsnd
  real(amrex_real) :: pi

    ascalar=a
    bscalar=b

   !print*,"aaaaa"
  
  pi = 4.d0 * atan(1.d0)

  do k = lo(3), hi(3)
       z = prob_lo(3) + dx(3) * (dble(k)+0.5d0)
       do j = lo(2), hi(2)
          y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
          do i = lo(1), hi(1)
             x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)
                
     	     PA=0.4d0*(rhoE(i,j,k)-0.5*(rhoU(i,j,k)**2+rhoV(i,j,k)**2+rhoW(i,j,k)**2)/RHO(i,j,k))
	     P_n=0.4d0*(rhoE_n(i,j,k)-0.5*(rhoU_n(i,j,k)**2+rhoV_n(i,j,k)**2+rhoW_n(i,j,k)**2)/RHO_n(i,j,k))
	     spdsnd=(1.4d0*P_n/RHO_n(i,j,k))**0.5d0
             alpha(i,j,k) = 1.0d0/(RHO_n(i,j,k)*spdsnd**2*dt_step**2)
	     rhs(i,j,k) = PA/(RHO_n(i,j,k)*spdsnd**2*dt_step**2)-( (Uface(i+1,j,k)-Uface(i,j,k))/(dx(1)*dt_step)+&
						  		   (Vface(i,j+1,k)-Vface(i,j,k))/(dx(2)*dt_step)+&
						  		   (Wface(i,j,k+1)-Wface(i,j,k))/(dx(3)*dt_step) )
	
             !alpha(i,j,k) = y**3
	     !rhs(i,j,k) = 1.0!ascalar*alpha(i,j,k)*cos(pi*x)+bscalar*pi*(x*pi*cos(pi*x)+sin(pi*x))
	     !rhs(i,j,k) = ascalar*alpha(i,j,k)*cos(pi*y)+bscalar*pi*(y*pi*cos(pi*y)+sin(pi*y))

          end do
       end do
    end do

    do k = glo(3), ghi(3)
       z = prob_lo(3) + dx(3) * (dble(k)+0.5d0)
       do j = glo(2), ghi(2)
          y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
          do i = glo(1), ghi(1)
             x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)
             beta(i,j,k) = RHO_n(i,j,k)
             !beta(i,j,k) = 1.0d0/y
          end do
       end do
    end do
 

!   print*, dt_step

  end subroutine helmholtz_rhs

end module rhs_helmholtz
