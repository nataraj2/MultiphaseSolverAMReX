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
  real(amrex_real) :: pi
double precision :: droprad,delxcell,delycell,delzcell,delx,dely,delz,rad,radfine,vofcell,xfine,yfine,zfine
  integer :: icell,jcell,kcell,nfine,xcennum,ycennum,zcennum
  double precision :: xcen(5),ycen(5), zcen(5)

  droprad=0.00051d0
  delx=dx(1)
  dely=dx(2)
  delz=dx(3)

  nfine=10

  ! Create an array of droplets

  xcen(1)=-0.4d0
  xcen(2)=-0.2d0
  xcen(3)=0.00d0
  xcen(4)=0.2d0
  xcen(5)=0.4d0

  ycen(1)=-0.4d0
  ycen(2)=-0.2d0
  ycen(3)=0.0d0
  ycen(4)=0.2d0
  ycen(5)=0.4d0

  zcen(3)=0.002d0

do zcennum=3,3
  do xcennum=3,3
    do ycennum=3,3
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)

                Grho(i,j,k)=0.0d0
                GrhoU(i,j,k)=0.0d0
                GrhoV(i,j,k)=0.0d0
                GrhoW(i,j,k)=0.0d0
                GrhoE(i,j,k)=0.0d0
                Lrho(i,j,k)=0.0d0
                LrhoU(i,j,k)=0.0d0
                LrhoV(i,j,k)=0.0d0
                LrhoW(i,j,k)=0.0d0
                LrhoE(i,j,k)=0.0d0
                VOF(i,j,k)=0.0d0

                x=prob_lo(1)+(dble(i)+0.5d0)*dx(1)
                y=prob_lo(2)+(dble(j)+0.5d0)*dx(2)
                z=prob_lo(3)+(dble(k)+0.5d0)*dx(3)

	    rad=dsqrt((prob_lo(1)+(dble(i)+0.5d0)*dx(1)-xcen(xcennum))**2+(prob_lo(2)+(dble(j)+0.5d0)*dx(2)-ycen(ycennum))**2+(prob_lo(3)+(dble(k)+0.5d0)*dx(3)-zcen(zcennum))**2)
            if(abs(droprad-rad)<=dsqrt(delx**2+dely**2+delz**2))then
                delxcell=dx(1)/nfine
                delycell=dx(2)/nfine
                delzcell=dx(3)/nfine
                vofcell=0.0
              do kcell=1,nfine
                do jcell=1,nfine
                   do icell=1,nfine
                        xfine=prob_lo(1)+dble(i)*dx(1)+(icell-1)*delxcell+delxcell/2.0
                        yfine=prob_lo(2)+dble(j)*dx(2)+(jcell-1)*delycell+delycell/2.0
                        zfine=prob_lo(3)+dble(k)*dx(3)+(kcell-1)*delzcell+delzcell/2.0
                        radfine=dsqrt((xfine-xcen(xcennum))**2+(yfine-ycen(ycennum))**2+(zfine-zcen(zcennum))**2)
                        if(droprad-radfine>=dsqrt(delxcell**2+delycell**2+delzcell**2))then
                                vofcell=vofcell+delxcell*delycell*delzcell
                        endif
                    enddo
                enddo
              enddo
                VOF(i,j,k)=vofcell/(delx*dely*delz)
             endif

             if(droprad-rad>=dsqrt(delx**2+dely**2+delz**2))then
                VOF(i,j,k)=1.0d0
             endif

		!VOF(i,j,k)=1.0d0
		!if(x.le.0.0025d0)VOF(i,j,k)=1.0d0		
		!if(rad.le.droprad)VOF(i,j,k)=1.0d0		

		Lrho(i,j,k)=1000.0d0
		if(rad.le.droprad)then
		LrhoU(i,j,k)=Lrho(i,j,k)*20.0d0*(1.0d0-rad/droprad)**(1.0d0/7.0d0)
		endif
		LrhoW(i,j,k)=Lrho(i,j,k)*0.0d0
		LrhoE(i,j,k)=(2.8d4+4.4d0*6.0d8)/(3.4d0)+0.5d0*(LrhoU(i,j,k)**2+LrhoV(i,j,k)**2+LrhoW(i,j,k)**2)/Lrho(i,j,k)

                Grho(i,j,k)=0.3d0
                GrhoU(i,j,k)=Grho(i,j,k)*0.0d0
		!if(z.le.-0.001d0)then
                GrhoW(i,j,k)=Grho(i,j,k)*700.0d0
		!endif
		GrhoE(i,j,k)=2.8d4/(0.4d0)+0.5d0*(GrhoU(i,j,k)**2+GrhoV(i,j,k)**2+GrhoW(i,j,k)**2)/Grho(i,j,k)

            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
 end subroutine init_prob_data

  subroutine compute_helmholtz_rhs
   
    implicit none

    integer :: ilev
    integer :: clo(4), chi(4), borlo(4), borhi(4), fxlo(4), fxhi(4), fylo(4), fyhi(4), fzlo(4), fzhi(4), betalo(4), betahi(4)
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: phiconsptr_n,phiconsptr_new,phi,phi_nborptr,Ufaceptr,Vfaceptr,Wfaceptr,rhsptr,alphaptr
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

      do idim = 1, amrex_spacedim
          nodal = .false.
          nodal(idim) = .true.
          call amrex_multifab_build(gradPmfab(idim,ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, 1, 0, nodal)
       end do 
   enddo
 
    do ilev = 0, nlevs

       call amrex_multifab_build(phi_mp_nbormfab, phi_n(ilev)%ba, phi_n(ilev)%dm, ncomp, ngrow)
       call fillpatch (ilev, time, phi_mp_nbormfab, phi_mp_n, phi_mp_n)

       call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.true.)
       do while (mfi%next())
          bx = mfi%tilebox()
	  phiconsptr_n   => phi_n(ilev)%dataptr(mfi)
	  phiconsptr_new   => phi_new(ilev)%dataptr(mfi)
          phi          => phi_mp(ilev)%dataptr(mfi)
          !phi_nborptr  => phi_mp_nbormfab%dataptr(mfi)
          phi_nborptr  => phi_mp_n(ilev)%dataptr(mfi)
	  alphaptr     => acoef(ilev)%dataptr(mfi)
	  Ufaceptr   => facevelmfab(1,ilev)%dataptr(mfi)
	  Vfaceptr   => facevelmfab(2,ilev)%dataptr(mfi)
	  Wfaceptr   => facevelmfab(3,ilev)%dataptr(mfi)
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

          call helmholtz_rhs(ilev, nlevs, bx%lo, bx%hi,&
			     phi(:,:,:,1),phi(:,:,:,2),phi(:,:,:,3),phi(:,:,:,4),phi(:,:,:,5),phiconsptr_n(:,:,:,11),phiconsptr_new(:,:,:,11),alphaptr,clo,chi,&
			     phi_nborptr(:,:,:,1),phi_nborptr(:,:,:,2),phi_nborptr(:,:,:,3),phi_nborptr(:,:,:,4),phi_nborptr(:,:,:,5),borlo,borhi,&
			     Ufaceptr,fxlo,fxhi,Vfaceptr,fylo,fyhi,Wfaceptr,fzlo,fzhi,rhsptr,amrex_geom(ilev)%dx, amrex_problo)

       end do
       call amrex_mfiter_destroy(mfi)

       call amrex_multifab_destroy(phi_mp_nbormfab)

       call Pmfab(ilev)%setVal(0.0_amrex_real)

    end do
  

  end subroutine compute_helmholtz_rhs


  subroutine helmholtz_rhs(level,max_level,lo,hi,RHO,rhoU,rhoV,rhoW,rhoE,VOF_n,VOF_new,alpha,clo,chi,&
			   RHO_n,rhoU_n,rhoV_n,rhoW_n,rhoE_n,borlo,borhi,&
			   Uface,fxlo,fxhi,Vface,fylo,fyhi,Wface,fzlo,fzhi,rhs,dx,prob_lo)
 
  implicit none
  integer, intent(in) :: level, max_level, lo(3), hi(3), clo(4), chi(4), borlo(4), borhi(4), fxlo(4), fxhi(4), fylo(4), fyhi(4), fzlo(4), fzhi(4)
  double precision, dimension(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)) :: RHO,rhoU,rhoV,rhoW,rhoE,VOF_n,VOF_new,alpha,rhs
  double precision, dimension(borlo(1):borhi(1),borlo(2):borhi(2),borlo(3):borhi(3)) :: RHO_n,rhoU_n,rhoV_n,rhoW_n,rhoE_n
  double precision, dimension(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3)) :: Uface
  double precision, dimension(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3)) :: Vface
  double precision, dimension(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3)) :: Wface
  double precision, intent(in) :: dx(3), prob_lo(3)
  integer          :: i,j,k
  double precision :: x,y,z
  real(amrex_real) :: PA, P_n, spdsnd
  real(amrex_real) :: pi
  double precision, parameter :: gamma_g=1.4d0, gamma_l=4.4d0, Pref_g=0.0d0, Pref_l=6.0d8
  double precision :: gamma_mix_n, gamma_mix_new, Pref_n, Pref_new, P1, P2, P3, P4, P5, P6, P7, P8, P9, P10, P11, P12
  double precision :: gamma_mix_n1, gamma_mix_n2, gamma_mix_n3, gamma_mix_n4, gamma_mix_n5, gamma_mix_n6

    ascalar=a
    bscalar=b

  pi = 4.d0 * atan(1.d0)

  do k = lo(3), hi(3)
       z = prob_lo(3) + dx(3) * (dble(k)+0.5d0)
       do j = lo(2), hi(2)
          y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
          do i = lo(1), hi(1)
             x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)
               
	     ! Mixing rules from Shukla's 2010 JCP
             gamma_mix_n=1.0d0+(gamma_l-1.0d0)*(gamma_g-1.0d0)/(VOF_n(i,j,k)*(gamma_g-1.0d0)+(1.0d0-VOF_n(i,j,k))*(gamma_l-1.0d0))
             gamma_mix_new=1.0d0+(gamma_l-1.0d0)*(gamma_g-1.0d0)/(VOF_new(i,j,k)*(gamma_g-1.0d0)+(1.0d0-VOF_new(i,j,k))*(gamma_l-1.0d0))

             Pref_n=(gamma_mix_n-1.0d0)/gamma_mix_n*(VOF_n(i,j,k)*gamma_l*Pref_l/(gamma_l-1.0d0)+(1.0d0-VOF_n(i,j,k))*gamma_g*Pref_g/(gamma_g-1.0d0))
             Pref_new=(gamma_mix_new-1.0d0)/gamma_mix_new*(VOF_new(i,j,k)*gamma_l*Pref_l/(gamma_l-1.0d0)+(1.0d0-VOF_new(i,j,k))*gamma_g*Pref_g/(gamma_g-1.0d0))

	     PA=(gamma_mix_new-1.0d0)*(rhoE(i,j,k)-0.5d0*(rhoU(i,j,k)**2+rhoV(i,j,k)**2+rhoW(i,j,k)**2)/RHO(i,j,k))-gamma_mix_new*Pref_new
	     P_n=(gamma_mix_n-1.0d0)*(rhoE_n(i,j,k)-0.5d0*(rhoU_n(i,j,k)**2+rhoV_n(i,j,k)**2+rhoW_n(i,j,k)**2)/RHO_n(i,j,k))-gamma_mix_n*Pref_n

	     spdsnd=(gamma_mix_n*(P_n+Pref_n)/RHO_n(i,j,k))**0.5d0

		if(level.eq.max_level)then
		if(gamma_mix_n*(P_n+Pref_n)/RHO_n(i,j,k).le.0.0d0)then
			print*, rhoE_n(i,j,k), rhoU_n(i,j,k), rhoV_n(i,j,k), rhoW_n(i,j,k), gamma_mix_n, P_n, Pref_n, x, y, z

             		gamma_mix_n1=1.0d0+(gamma_l-1.0d0)*(gamma_g-1.0d0)/(VOF_n(i-1,j,k)*(gamma_g-1.0d0)+(1.0d0-VOF_n(i-1,j,k))*(gamma_l-1.0d0))
			P1=(gamma_mix_n1-1.0d0)*(rhoE_n(i-1,j,k)-0.5d0*(rhoU_n(i-1,j,k)**2+rhoV_n(i-1,j,k)**2+rhoW_n(i-1,j,k)**2)/RHO_n(i-1,j,k))-gamma_mix_n1*Pref_n
             		gamma_mix_n2=1.0d0+(gamma_l-1.0d0)*(gamma_g-1.0d0)/(VOF_n(i+1,j,k)*(gamma_g-1.0d0)+(1.0d0-VOF_n(i+1,j,k))*(gamma_l-1.0d0))
			P2=(gamma_mix_n2-1.0d0)*(rhoE_n(i+1,j,k)-0.5d0*(rhoU_n(i+1,j,k)**2+rhoV_n(i+1,j,k)**2+rhoW_n(i+1,j,k)**2)/RHO_n(i+1,j,k))-gamma_mix_n2*Pref_n
             		gamma_mix_n3=1.0d0+(gamma_l-1.0d0)*(gamma_g-1.0d0)/(VOF_n(i,j-1,k)*(gamma_g-1.0d0)+(1.0d0-VOF_n(i,j-1,k))*(gamma_l-1.0d0))
			P3=(gamma_mix_n3-1.0d0)*(rhoE_n(i,j-1,k)-0.5d0*(rhoU_n(i,j-1,k)**2+rhoV_n(i,j-1,k)**2+rhoW_n(i,j-1,k)**2)/RHO_n(i,j-1,k))-gamma_mix_n3*Pref_n
             		gamma_mix_n4=1.0d0+(gamma_l-1.0d0)*(gamma_g-1.0d0)/(VOF_n(i,j+1,k)*(gamma_g-1.0d0)+(1.0d0-VOF_n(i,j+1,k))*(gamma_l-1.0d0))
			P4=(gamma_mix_n4-1.0d0)*(rhoE_n(i,j+1,k)-0.5d0*(rhoU_n(i,j+1,k)**2+rhoV_n(i,j+1,k)**2+rhoW_n(i,j+1,k)**2)/RHO_n(i,j+1,k))-gamma_mix_n4*Pref_n
             		gamma_mix_n5=1.0d0+(gamma_l-1.0d0)*(gamma_g-1.0d0)/(VOF_n(i,j,k-1)*(gamma_g-1.0d0)+(1.0d0-VOF_n(i,j,k-1))*(gamma_l-1.0d0))
			P5=(gamma_mix_n5-1.0d0)*(rhoE_n(i,j,k-1)-0.5d0*(rhoU_n(i,j,k-1)**2+rhoV_n(i,j,k-1)**2+rhoW_n(i,j,k-1)**2)/RHO_n(i,j,k-1))-gamma_mix_n5*Pref_n
             		gamma_mix_n6=1.0d0+(gamma_l-1.0d0)*(gamma_g-1.0d0)/(VOF_n(i,j,k+1)*(gamma_g-1.0d0)+(1.0d0-VOF_n(i,j,k+1))*(gamma_l-1.0d0))
			P6=(gamma_mix_n6-1.0d0)*(rhoE_n(i,j,k+1)-0.5d0*(rhoU_n(i,j,k+1)**2+rhoV_n(i,j,k+1)**2+rhoW_n(i,j,k+1)**2)/RHO_n(i,j,k+1))-gamma_mix_n6*Pref_n

			!gamma_mix_n=1.0d0+(gamma_l-1.0d0)*(gamma_g-1.0d0)/(VOF_n(i-2,j,k)*(gamma_g-1.0d0)+(1.0d0-VOF_n(i-2,j,k))*(gamma_l-1.0d0))
			!P7=(gamma_mix_n-1.0d0)*(rhoE_n(i-2,j,k)-0.5d0*(rhoU_n(i-2,j,k)**2+rhoV_n(i-2,j,k)**2+rhoW_n(i-2,j,k)**2)/RHO_n(i-2,j,k))-gamma_mix_n*Pref_n
             		!gamma_mix_n=1.0d0+(gamma_l-1.0d0)*(gamma_g-1.0d0)/(VOF_n(i+2,j,k)*(gamma_g-1.0d0)+(1.0d0-VOF_n(i+2,j,k))*(gamma_l-1.0d0))
			!P8=(gamma_mix_n-1.0d0)*(rhoE_n(i+2,j,k)-0.5d0*(rhoU_n(i+2,j,k)**2+rhoV_n(i+2,j,k)**2+rhoW_n(i+2,j,k)**2)/RHO_n(i+2,j,k))-gamma_mix_n*Pref_n
             		!gamma_mix_n=1.0d0+(gamma_l-1.0d0)*(gamma_g-1.0d0)/(VOF_n(i,j-2,k)*(gamma_g-1.0d0)+(1.0d0-VOF_n(i,j-2,k))*(gamma_l-1.0d0))
			!P9=(gamma_mix_n-1.0d0)*(rhoE_n(i,j-2,k)-0.5d0*(rhoU_n(i,j-2,k)**2+rhoV_n(i,j-2,k)**2+rhoW_n(i,j-2,k)**2)/RHO_n(i,j-2,k))-gamma_mix_n*Pref_n
             		!gamma_mix_n=1.0d0+(gamma_l-1.0d0)*(gamma_g-1.0d0)/(VOF_n(i,j+2,k)*(gamma_g-1.0d0)+(1.0d0-VOF_n(i,j+2,k))*(gamma_l-1.0d0))
			!P10=(gamma_mix_n-1.0d0)*(rhoE_n(i,j+2,k)-0.5d0*(rhoU_n(i,j+2,k)**2+rhoV_n(i,j+2,k)**2+rhoW_n(i,j+2,k)**2)/RHO_n(i,j+2,k))-gamma_mix_n*Pref_n
             		!gamma_mix_n=1.0d0+(gamma_l-1.0d0)*(gamma_g-1.0d0)/(VOF_n(i,j,k-2)*(gamma_g-1.0d0)+(1.0d0-VOF_n(i,j,k-2))*(gamma_l-1.0d0))
			!P11=(gamma_mix_n-1.0d0)*(rhoE_n(i,j,k-2)-0.5d0*(rhoU_n(i,j,k-2)**2+rhoV_n(i,j,k-2)**2+rhoW_n(i,j,k-2)**2)/RHO_n(i,j,k-2))-gamma_mix_n*Pref_n
             		!gamma_mix_n=1.0d0+(gamma_l-1.0d0)*(gamma_g-1.0d0)/(VOF_n(i,j,k+2)*(gamma_g-1.0d0)+(1.0d0-VOF_n(i,j,k+2))*(gamma_l-1.0d0))
			!P12=(gamma_mix_n-1.0d0)*(rhoE_n(i,j,k+2)-0.5d0*(rhoU_n(i,j,k+2)**2+rhoV_n(i,j,k+2)**2+rhoW_n(i,j,k+2)**2)/RHO_n(i,j,k+2))-gamma_mix_n*Pref_n

			if(P1.gt.0.0d0)then
			P_n=P1
			elseif(P2.gt.0.0d0)then
			P_n=P2
			elseif(P3.gt.0.0d0)then
			P_n=P3
			elseif(P4.gt.0.0d0)then
			P_n=P4
			elseif(P5.gt.0.0d0)then
			P_n=P5
			elseif(P6.gt.0.0d0)then
			P_n=P6
			!elseif(P7.gt.0.0d0)then
			!P_n=P7
			!elseif(P8.gt.0.0d0)then
			!P_n=P8
			!elseif(P9.gt.0.0d0)then
			!P_n=P9
			!elseif(P10.gt.0.0d0)then
			!P_n=P10
			!elseif(P11.gt.0.0d0)then
			!P_n=P11
			!elseif(P12.gt.0.0d0)then
			!P_n=P12
			endif
			spdsnd=(gamma_mix_n*(P_n+Pref_n)/RHO_n(i,j,k))**0.5d0
			!stop
		endif
		endif

		if(level.eq.max_level)then
                if(gamma_mix_n*(P_n+Pref_n)/RHO_n(i,j,k).le.0.0d0)then
		print*, "Still a problem", rhoE_n(i,j,k), rhoU_n(i,j,k), rhoV_n(i,j,k), rhoW_n(i,j,k), gamma_mix_n, P_n, Pref_n, x, y, z
		!stop
		endif
		endif
	
             !P(i,j,k)=(gamma_mix-1.0d0)*(rhoE(i,j,k)-0.5d0*RHO(i,j,k)*(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2))-gamma_mix*my_Pref
             !spdsnd=sqrt(gamma_mix*(P(i,j,k)+my_Pref)/RHO(i,j,k))
 
             alpha(i,j,k) = 1.0d0/(RHO_n(i,j,k)*spdsnd**2*dt_step**2)
	     rhs(i,j,k) = PA/(RHO_n(i,j,k)*spdsnd**2*dt_step**2)-( (Uface(i+1,j,k)-Uface(i,j,k))/(dx(1)*dt_step)+&
						  		   (Vface(i,j+1,k)-Vface(i,j,k))/(dx(2)*dt_step)+&
						  		   (Wface(i,j,k+1)-Wface(i,j,k))/(dx(3)*dt_step) )

          end do
       end do
    end do

  end subroutine helmholtz_rhs

end module rhs_helmholtz
