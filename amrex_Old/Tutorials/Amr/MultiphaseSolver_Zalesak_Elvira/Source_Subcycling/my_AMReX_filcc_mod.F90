
#include "AMReX_BC_TYPES.H"

module my_amrex_filcc_module

  use amrex_fort_module, only : amrex_real, amrex_spacedim!, get_loop_bounds
  use amrex_constants_module

  implicit none

  interface my_amrex_filcc
     module procedure my_amrex_filcc_n
  end interface my_amrex_filcc

  private
  public :: my_amrex_filcc, my_filccn, my_amrex_hoextraptocc

contains

  subroutine my_amrex_filcc_n(q,qlo,qhi,domlo,domhi,dx,xlo,bclo,bchi)
    integer, intent(in) :: qlo(4), qhi(4)
    integer, dimension(amrex_spacedim), intent(in) :: domlo, domhi
    real(amrex_real), intent(in) :: dx(amrex_spacedim), xlo(amrex_spacedim)
    integer, intent(in) :: bclo(amrex_spacedim,*), bchi(amrex_spacedim,*)
    real(amrex_real), intent(inout) :: q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),qlo(4):qhi(4))
    integer :: i, bc(amrex_spacedim,2)

   
    do i = qlo(4), qhi(4)
       bc(:,1) = bclo(:,i)
       bc(:,2) = bchi(:,i)
#if (AMREX_SPACEDIM == 3)
       call my_filcc(i,q(:,:,:,i),qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3),domlo,domhi,dx,xlo,bc)
#elif (AMREX_SPACEDIM == 2)
       call my_filcc(i,q(:,:,:,i),qlo(1),qlo(2),       qhi(1),qhi(2),       domlo,domhi,dx,xlo,bc)
#else
       call my_filcc(i,q(:,:,:,i),qlo(1),              qhi(1),              domlo,domhi,dx,xlo,bc)
#endif
    end do
  end subroutine my_amrex_filcc_n

  subroutine my_filccn(compno,lo, hi, q, q_lo, q_hi, ncomp, domlo, domhi, dx, xlo, bc)

    implicit none

    integer,          intent(in   ) :: lo(3), hi(3)
    integer,          intent(in   ) :: q_lo(3), q_hi(3)
    integer,          intent(in   ) :: ncomp
    integer,          intent(in   ) :: domlo(amrex_spacedim), domhi(amrex_spacedim)
    real(amrex_real), intent(in   ) :: xlo(amrex_spacedim), dx(amrex_spacedim)
    real(amrex_real), intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),ncomp)
    integer,          intent(in   ) :: bc(amrex_spacedim,2,ncomp)

    integer :: ilo, ihi, jlo, jhi, klo, khi
    integer :: is, ie, js, je, ks, ke
    integer :: i, j, k, n
    integer :: imin, imax, jmin, jmax, kmin, kmax
    integer :: compno

    real(amrex_real) :: x, y, z  
    real(amrex_real) :: Tinfty, Tj, Vj, Mj, rhoj, Vz, afac, rad, jetrad, Grho, GrhoU, GrhoE

    is = max(q_lo(1), domlo(1))
    ie = min(q_hi(1), domhi(1))
    ilo = domlo(1)
    ihi = domhi(1)

#if AMREX_SPACEDIM >= 2
    js = max(q_lo(2), domlo(2))
    je = min(q_hi(2), domhi(2))
    jlo = domlo(2)
    jhi = domhi(2)
#endif

#if AMREX_SPACEDIM == 3
    ks = max(q_lo(3), domlo(3))
    ke = min(q_hi(3), domhi(3))
    klo = domlo(3)
    khi = domhi(3)
#endif

	!print*, lo(1), ilo
	!stop
	!print*, ncomp
	!stop

    do n = 1, ncomp

       if (lo(1) < ilo) then
          imin = lo(1)
          imax = min(hi(1),ilo-1)

		!print*, imin, imax, ilo
		!stop

          if (bc(1,1,n) .eq. EXT_DIR) then

		!print*, "here...", xlo(2), xlo(3), dx(2), dx(3)

	Tinfty=300.0d0
        Tj=300.0d0*1.0d0
        Vj=400.0d0
        Mj=Vj/dsqrt(1.4d0*287.0d0*Tj)
        afac=1.0d0
	jetrad=0.5d0
	rhoj=101325.0d0/(287.0d0*Tj)

             ! Do nothing.
	    do k = lo(3), hi(3)
	       z=0.0d0+(dble(k)+0.5d0) * dx(3)
                do j = lo(2), hi(2)
		   y=-5.0d0+(dble(j)+0.5d0) * dx(2)

		   rad = dsqrt(y**2+(z-5.0d0)**2)
	           Vz=0.5d0*(1.0d0+tanh(afac*(jetrad/rad-rad/jetrad)))

                   !do i = ilo, ilo
                   do i = imin, imax
		  
		      Grho = rhoj*(0.5d0*0.4d0*Vz*(1.0d0-Vz)*Mj**2+Vz+Tinfty/Tj*(1.0d0-Vz))**(-1) 
		      GrhoU = Grho*Vz*Vj
		      GrhoE = 101325.0/0.4d0 + 0.5d0*GrhoU**2/Grho
 
		      if(compno.eq.1.or.compno.eq.6)then
                      q(i,j,k,n) = 2.0d0*Grho - q(ilo,j,k,n)
			elseif(compno.eq.2.or.compno.eq.7)then
			q(i,j,k,n) = 2.0d0*GrhoU - q(ilo,j,k,n)
			elseif(compno.eq.5.or.compno.eq.10)then
		      q(i,j,k,n) = 2.0d0*(101325.0d0/0.4d0+0.5d0*GrhoU**2/Grho) - q(ilo,j,k,n)
			elseif(compno.eq.11)then
			if(rad.le.jetrad)then
			q(i,j,k,n)=2.0d0-q(ilo,j,k,n)
			else
			q(i,j,k,n)=0.0d0
			endif
			else
			q(i,j,k,n) = q(ilo,j,k,n)
			endif

			!if(compno.eq.2)q(i,j,k,n)=10000.0d0

			!if(compno.eq.6)then
                      !q(i,j,k,n) = 2.0d0*Grho - q(ilo,j,k,n)
			!elseif(compno.eq.7)then
			!q(i,j,k,n) = 2.0d0*GrhoU - q(ilo,j,k,n)
			!elseif(compno.eq.10)then
		      !q(i,j,k,n) = 2.0d0*(101325.0d0/0.4d0+0.5d0*GrhoU**2/Grho) - q(ilo,j,k,n)
			!elseif(compno.eq.8.or.compno.eq.9)then
			!q(i,j,k,n) = q(ilo,j,k,n)
			!elseif(compno.le.5)then
			!q(i,j,k,n)=1000.0d0
			!else
			!q(i,j,k,n)=0.0d0
			!endif


                   end do
                end do
             end do

	    !do k = lo(3), hi(3)
             !   do j = lo(2), hi(2)
              !     do i = imin, imax
               !       q(i,j,k,n) = q(ilo,j,k,n)
               !    end do
               ! end do
             !end do
		

          else if (bc(1,1,n) .eq. FOEXTRAP) then
             
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = q(ilo,j,k,n)
                   end do
                end do
             end do
             
          else if (bc(1,1,n) .eq. HOEXTRAP) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax

                      if (i < ilo - 1) then
                         q(i,j,k,n) = q(ilo,j,k,n)
                      else if (i == ilo - 1) then
                         if (ilo+2 <= ie) then
                            q(i,j,k,n) = eighth * (15*q(ilo,j,k,n) - 10*q(ilo+1,j,k,n) + 3*q(ilo+2,j,k,n))
                         else
                            q(i,j,k,n) = half * (3*q(ilo,j,k,n) - q(ilo+1,j,k,n))
                         end if
                      end if

                   end do
                end do
             end do
             
          else if (bc(1,1,n) .eq. REFLECT_EVEN) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = q(ilo+(ilo-i)-1,j,k,n)
                   end do
                end do
             end do

          else if (bc(1,1,n) .eq. REFLECT_ODD) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = -q(ilo+(ilo-i)-1,j,k,n)
                   end do
                end do
             end do
             
          end if

       end if

       if (hi(1) > ihi) then
          imin = max(lo(1),ihi+1)
          imax = hi(1)

          if (bc(1,2,n) .eq. EXT_DIR) then

             ! Do nothing.
	    do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = ihi, ihi
                      q(i,j,k,n) = 0.5!q(ihi,j,k,n)
                   end do
                end do
             end do

          else if (bc(1,2,n) .eq. FOEXTRAP) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = q(ihi,j,k,n)
                   end do
                end do
             end do
             
          else if (bc(1,2,n) .eq. HOEXTRAP) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax

                      if (i > ihi + 1) then
                         q(i,j,k,n) = q(ihi,j,k,n)
                      else if (i == ihi + 1) then
                         if (ihi-2 >= is) then
                            q(i,j,k,n) = eighth * (15*q(ihi,j,k,n) - 10*q(ihi-1,j,k,n) + 3*q(ihi-2,j,k,n))
                         else
                            q(i,j,k,n) = half * (3*q(ihi,j,k,n) - q(ihi-1,j,k,n))
                         end if
                      end if

                   end do
                end do
             end do

          else if (bc(1,2,n) .eq. REFLECT_EVEN) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = q(ihi-(i-ihi)+1,j,k,n)
                   end do
                end do
             end do
             
          else if (bc(1,2,n) .eq. REFLECT_ODD) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = -q(ihi-(i-ihi)+1,j,k,n)
                   end do
                end do
             end do
             
          end if

       end if


#if AMREX_SPACEDIM >= 2

       if (lo(2) < jlo) then
          jmin = lo(2)
          jmax = min(hi(2),jlo-1)

          if (bc(2,1,n) .eq. EXT_DIR) then

             ! Do nothing.
	    do k = lo(3), hi(3)
                do j = jlo, jlo
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = 0.8!q(i,jlo,k,n)
                   end do
                end do
             end do

          else if (bc(2,1,n) .eq. FOEXTRAP) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,jlo,k,n)
                   end do
                end do
             end do
             
          else if (bc(2,1,n) .eq. HOEXTRAP) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)

                      if (j < jlo - 1) then
                         q(i,j,k,n) = q(i,jlo,k,n)
                      else if (j == jlo - 1) then
                         if (jlo+2 <= je) then
                            q(i,j,k,n) = eighth * (15*q(i,jlo,k,n) - 10*q(i,jlo+1,k,n) + 3*q(i,jlo+2,k,n))
                         else
                            q(i,j,k,n) = half * (3*q(i,jlo,k,n) - q(i,jlo+1,k,n))
                         end if
                      end if

                   end do
                end do
             end do

          else if (bc(2,1,n) .eq. REFLECT_EVEN) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,jlo+(jlo-j)-1,k,n)
                   end do
                end do
             end do

          else if (bc(2,1,n) .eq. REFLECT_ODD) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = -q(i,jlo+(jlo-j)-1,k,n)
                   end do
                end do
             end do
             
          end if

       end if

       if (hi(2) > jhi) then
          jmin = max(lo(2),jhi+1)
          jmax = hi(2)

          if (bc(2,2,n) .eq. EXT_DIR) then

             ! Do nothing.
	   do k = lo(3), hi(3)
                do j = jhi, jhi
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = 0.6!q(i,jhi,k,n)
                   end do
                end do
             end do

          else if (bc(2,2,n) .eq. FOEXTRAP) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,jhi,k,n)
                   end do
                end do
             end do
             
          else if (bc(2,2,n) .eq. HOEXTRAP) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)

                      if (j > jhi + 1) then
                         q(i,j,k,n) = q(i,jhi,k,n)
                      else if (j == jhi + 1) then
                         if (jhi-2 >= js) then
                            q(i,j,k,n) = eighth * (15*q(i,jhi,k,n) - 10*q(i,jhi-1,k,n) + 3*q(i,jhi-2,k,n))
                         else
                            q(i,j,k,n) = half * (3*q(i,jhi,k,n) - q(i,jhi-1,k,n))
                         end if
                      end if
                      
                   end do
                end do
             end do

          else if (bc(2,2,n) .eq. REFLECT_EVEN) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,jhi-(j-jhi)+1,k,n)
                   end do
                end do
             end do
             
          else if (bc(2,2,n) .eq. REFLECT_ODD) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = -q(i,jhi-(j-jhi)+1,k,n)
                   end do
                end do
             end do
             
          end if

       end if
#endif



#if AMREX_SPACEDIM == 3

       if (lo(3) < klo) then
          kmin = lo(3)
          kmax = min(hi(3),klo-1)

          if (bc(3,1,n) .eq. EXT_DIR) then

             ! Do nothing.
             do k = klo, klo
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = 0.5!q(i,j,klo,n)
                   end do
                end do
             end do

          else if (bc(3,1,n) .eq. FOEXTRAP) then
             
             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,j,klo,n)
                   end do
                end do
             end do

          else if (bc(3,1,n) .eq. HOEXTRAP) then

             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)

                      if (k < klo - 1) then
                         q(i,j,k,n) = q(i,j,klo,n)
                      else if (k == klo - 1) then
                         if (klo+2 <= ke) then
                            q(i,j,k,n) = eighth * (15*q(i,j,klo,n) - 10*q(i,j,klo+1,n) + 3*q(i,j,klo+2,n))
                         else
                            q(i,j,k,n) = half * (3*q(i,j,klo,n) - q(i,j,klo+1,n))
                         end if
                      end if

                   end do
                end do
             end do
             
          else if (bc(3,1,n) .eq. REFLECT_EVEN) then

             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,j,klo+(klo-k)-1,n)
                   end do
                end do
             end do
             
          else if (bc(3,1,n) .eq. REFLECT_ODD) then

             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = -q(i,j,klo+(klo-k)-1,n)
                   end do
                end do
             end do
             
          end if

       end if

       if (hi(3) > khi) then
          kmin = max(lo(3),khi+1)
          kmax = hi(3)

          if (bc(3,2,n) .eq. EXT_DIR) then

             ! Do nothing.
              do k = khi, khi
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = 0.7!q(i,j,khi,n)
                   end do
                end do
             end do
  
          else if (bc(3,2,n) .eq. FOEXTRAP) then
             
             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,j,khi,n)
                   end do
                end do
             end do

          else if (bc(3,2,n) .eq. HOEXTRAP) then

             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)

                      if (k > khi + 1) then
                         q(i,j,k,n) = q(i,j,khi,n)
                      else if (k == khi + 1) then
                         if (khi-2 >= ks) then
                            q(i,j,k,n) = eighth * (15*q(i,j,khi,n) - 10*q(i,j,khi-1,n) + 3*q(i,j,khi-2,n))
                         else
                            q(i,j,k,n) = half * (3*q(i,j,khi,n) - q(i,j,khi-1,n))
                         end if
                      end if

                   end do
                end do
             end do
             
          else if (bc(3,2,n) .eq. REFLECT_EVEN) then

             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,j,khi-(k-khi)+1,n)
                   end do
                end do
             end do
             
          else if (bc(3,2,n) .eq. REFLECT_ODD) then

             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = -q(i,j,khi-(k-khi)+1,n)
                   end do
                end do
             end do

          end if

       end if
#endif



#if AMREX_SPACEDIM >= 2
       ! Now take care of the higher contributions

       !
       ! First correct the i-j edges and all corners
       !

       if (bc(1,1,n) .eq. HOEXTRAP .and. bc(2,1,n) .eq. HOEXTRAP) then

          if (lo(1) < ilo .and. lo(2) < jlo) then
             imin = lo(1)
             imax = min(hi(1),ilo-1)
             jmin = lo(2)
             jmax = min(hi(2),jlo-1)

             i = ilo-1
             j = jlo-1

             if (i.ge.imin .and. i.le.imax .and. j.ge.jmin .and. j.le.jmax) then

                do k = lo(3), hi(3)
                      
                   if (jlo+2 <= je) then
                      q(i,j,k,n) = half * eighth * (15*q(ilo-1,jlo,k,n) - 10*q(ilo-1,jlo+1,k,n) + 3*q(ilo-1,jlo+2,k,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(ilo-1,jlo,k,n) - q(ilo-1,jlo+1,k,n))
                   end if
                   
                   if (ilo+2 <= ie) then
                      q(i,j,k,n) = q(ilo-1,jlo-1,k,n) + &
                           half * eighth * (15*q(ilo,jlo-1,k,n) - 10*q(ilo+1,jlo-1,k,n) + 3*q(ilo+2,jlo-1,k,n))
                   else
                      q(i,j,k,n) = q(ilo-1,jlo-1,k,n) + half * half * (3*q(ilo,jlo-1,k,n) - q(ilo+1,jlo-1,k,n))
                   end if
                   
#if AMREX_SPACEDIM == 3
                   
                   if (k == klo-1 .and. bc(3,1,n) .eq. HOEXTRAP) then
                      if (klo+2 <= ke) then
                         q(i,j,k,n) = eighth * ( (15*q(ilo-1,jlo-1,klo,n) - 10*q(ilo-1,jlo-1,klo+1,n) + &
                              3*q(ilo-1,jlo-1,klo+2,n)) )
                      else
                         q(i,j,k,n) = half * (3*q(ilo-1,jlo-1,klo,n) - q(ilo-1,jlo-1,klo+1,n))
                      end if
                   end if
                   
                   if (k == khi+1 .and. bc(3,2,n) .eq. HOEXTRAP) then
                      if (khi-2 >= ks) then
                         q(i,j,k,n) = eighth * ( (15*q(ilo-1,jlo-1,khi,n) - 10*q(ilo-1,jlo-1,khi-1,n) + &
                              3*q(ilo-1,jlo-1,khi-2,n)) )
                      else
                         q(i,j,k,n) = half * (3*q(ilo-1,jlo-1,khi,n) - q(ilo-1,jlo-1,khi-1,n))
                      end if
                   end if
#endif
                   
                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,1,n) .eq. HOEXTRAP .and. bc(2,2,n) .eq. HOEXTRAP) then

          if (lo(1) < ilo .and. hi(2) > jhi) then
             imin = lo(1)
             imax = min(hi(1),ilo-1)
             jmin = max(lo(2),jhi+1)
             jmax = hi(2)

             i = ilo-1
             j = jhi+1

             if (i.ge.imin .and. i.le.imax .and. j.ge.jmin .and. j.le.jmax) then

                do k = lo(3), hi(3)

                   if (jhi-2 >= js) then
                      q(i,j,k,n) = half * eighth * (15*q(ilo-1,jhi,k,n) - 10*q(ilo-1,jhi-1,k,n) + 3*q(ilo-1,jhi-2,k,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(ilo-1,jhi,k,n) - q(ilo-1,jhi-1,k,n))
                   end if
                   
                   if (ilo+2 <= ie) then
                      q(i,j,k,n) = q(ilo-1,jhi+1,k,n) + &
                           half * eighth * (15*q(ilo,jhi+1,k,n) - 10*q(ilo+1,jhi+1,k,n) + 3*q(ilo+2,jhi+1,k,n))
                   else
                      q(i,j,k,n) = q(ilo-1,jhi+1,k,n) + half * half * (3*q(ilo,jhi+1,k,n) - q(ilo+1,jhi+1,k,n))
                   end if

#if (AMREX_SPACEDIM == 3)
                   if (k == klo-1 .and. bc(3,1,n) .eq. HOEXTRAP) then
                      if (klo+2 <= ke) then
                         q(i,j,k,n) = eighth * ( (15*q(ilo-1,jhi+1,klo,n) - 10*q(ilo-1,jhi+1,klo+1,n) + &
                              3*q(ilo-1,jhi+1,klo+2,n)) )
                      else
                         q(i,j,k,n) = half * (3*q(ilo-1,jhi+1,klo,n) - q(ilo-1,jhi+1,klo+1,n))
                      end if
                   end if
                   
                   if (k == khi+1 .and. bc(3,2,n) .eq. HOEXTRAP) then
                      if (khi-2 >= ks) then
                         q(i,j,k,n) = eighth * ( (15*q(ilo-1,jhi+1,khi,n) - 10*q(ilo-1,jhi+1,khi-1,n) + &
                              3*q(ilo-1,jhi+1,khi-2,n)) )
                      else
                         q(i,j,k,n) = half * (3*q(ilo-1,jhi+1,khi,n) - q(ilo-1,jhi+1,khi-1,n))
                      end if
                   end if

#endif

                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,2,n) .eq. HOEXTRAP .and. bc(2,1,n) .eq. HOEXTRAP) then

          if (hi(1) > ihi .and. lo(2) < jlo) then
             imin = max(lo(1),ihi+1)
             imax = hi(1)
             jmin = lo(2)
             jmax = min(hi(2),jlo-1)

             i = ihi+1
             j = jlo-1

             if (i.ge.imin .and. i.le.imax .and. j.ge.jmin .and. j.le.jmax) then

                do k = lo(3), hi(3)

                   if (jlo+2 <= je) then
                      q(i,j,k,n) = half * eighth * (15*q(ihi+1,jlo,k,n) - 10*q(ihi+1,jlo+1,k,n) + 3*q(ihi+1,jlo+2,k,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(ihi+1,jlo,k,n) - q(ihi+1,jlo+1,k,n))
                   end if
                   
                   if (ihi-2 >= is) then
                      q(i,j,k,n) = q(ihi+1,jlo-1,k,n) + &
                           half * eighth * (15*q(ihi,jlo-1,k,n) - 10*q(ihi-1,jlo-1,k,n) + 3*q(ihi-2,jlo-1,k,n))
                   else
                      q(i,j,k,n) = q(ihi+1,jlo-1,k,n) + half * half * (3*q(ihi,jlo-1,k,n) - q(ihi-1,jlo-1,k,n))
                   end if
                   
#if (AMREX_SPACEDIM == 3)
                   if (k == klo-1 .and. bc(3,1,n) .eq. HOEXTRAP) then
                      if (klo+2 <= ke) then
                         q(i,j,k,n) = eighth * (15*q(ihi+1,jlo-1,klo,n) - 10*q(ihi+1,jlo-1,klo+1,n) + 3*q(ihi+1,jlo-1,klo+2,n))
                      else
                         q(i,j,k,n) = half * (3*q(ihi+1,jlo-1,klo,n) - q(ihi+1,jlo-1,klo+1,n))
                      end if
                   end if
                   
                   if (k == khi+1 .and. bc(3,2,n) .eq. HOEXTRAP) then
                      if (khi-2 >= ks) then
                         q(i,j,k,n) = eighth * (15*q(ihi+1,jlo-1,khi,n) - 10*q(ihi+1,jlo-1,khi-1,n) + 3*q(ihi+1,jlo-1,khi-2,n))
                      else
                         q(i,j,k,n) = half * (3*q(ihi+1,jlo-1,khi,n) - q(ihi+1,jlo-1,khi-1,n))
                      end if
                   end if
#endif

                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,2,n) .eq. HOEXTRAP .and. bc(2,2,n) .eq. HOEXTRAP) then

          if (hi(1) > ihi .and. hi(2) > jhi) then
             imin = max(lo(1),ihi+1)
             imax = hi(1)
             jmin = max(lo(2),jhi+1)
             jmax = hi(2)

             i = ihi+1
             j = jhi+1

             if (i.ge.imin .and. i.le.imax .and. j.ge.jmin .and. j.le.jmax) then

                do k = lo(3), hi(3)
                   
                   if (jhi-2 >= js) then
                      q(i,j,k,n) = half * eighth * (15*q(ihi+1,jhi,k,n) - 10*q(ihi+1,jhi-1,k,n) + 3*q(ihi+1,jhi-2,k,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(ihi+1,jhi,k,n) - q(ihi+1,jhi-1,k,n))
                   end if
                   
                   if (ihi-2 >= is) then
                      q(i,j,k,n) = q(ihi+1,jhi+1,k,n) + &
                           half * eighth * (15*q(ihi,jhi+1,k,n) - 10*q(ihi-1,jhi+1,k,n) + 3*q(ihi-2,jhi+1,k,n))
                   else
                      q(i,j,k,n) = q(ihi+1,jhi+1,k,n) + half * half * (3*q(ihi,jhi+1,k,n) - q(ihi-1,jhi+1,k,n))
                   end if
                   
#if (AMREX_SPACEDIM == 3)
                   if (k == klo-1 .and. bc(3,1,n) .eq. HOEXTRAP) then
                      if (klo+2 <= ke) then
                         q(i,j,k,n) = eighth * (15*q(ihi+1,jhi+1,klo,n) - 10*q(ihi+1,jhi+1,klo+1,n) + 3*q(ihi+1,jhi+1,klo+2,n))
                      else
                         q(i,j,k,n) = half * (3*q(ihi+1,jhi+1,klo,n) - q(ihi+1,jhi+1,klo+1,n))
                      end if
                   end if
                   
                   if (k == khi+1 .and. bc(3,2,n) .eq. HOEXTRAP) then
                      if (khi-2 >= ks) then
                         q(i,j,k,n) = eighth * (15*q(ihi+1,jhi+1,khi,n) - 10*q(ihi+1,jhi+1,khi-1,n) + 3*q(ihi+1,jhi+1,khi-2,n))
                      else
                         q(i,j,k,n) = half * (3*q(ihi+1,jhi+1,khi,n) - q(ihi+1,jhi+1,khi-1,n))
                      end if
                   end if
#endif
                   
                end do
             end if
          end if
       end if
#endif

#if AMREX_SPACEDIM == 3
       !
       ! Next correct the i-k edges
       !

       if (bc(1,1,n) .eq. HOEXTRAP .and. bc(3,1,n) .eq. HOEXTRAP) then

          if (lo(1) < ilo .and. lo(3) < klo) then
             imin = lo(1)
             imax = min(hi(1),ilo-1)
             kmin = lo(3)
             kmax = min(hi(3),klo-1)

             i = ilo-1
             k = klo-1

             if (i.ge.imin .and. i.le.imax .and. k.ge.kmin .and. k.le.kmax) then
             
                do j = lo(2), hi(2)

                   if (klo+2 <= ke) then
                      q(i,j,k,n) = half * eighth * (15*q(ilo-1,j,klo,n) - 10*q(ilo-1,j,klo+1,n) + 3*q(ilo-1,j,klo+2,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(ilo-1,j,klo,n) - q(ilo-1,j,klo+1,n))
                   end if
                   
                   if (ilo+2 <= ie) then
                      q(i,j,k,n) = q(ilo-1,j,klo-1,n) + &
                           half * eighth * (15*q(ilo,j,klo-1,n) - 10*q(ilo+1,j,klo-1,n) + 3*q(ilo+2,j,klo-1,n))
                   else
                      q(i,j,k,n) = q(ilo-1,j,klo-1,n) + half * half * (3*q(ilo,j,klo-1,n) - q(ilo+1,j,klo-1,n))
                   end if

                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,1,n) .eq. HOEXTRAP .and. bc(3,2,n) .eq. HOEXTRAP) then

          if (lo(1) < ilo .and. hi(3) > khi) then
             imin = lo(1)
             imax = min(hi(1),ilo-1)
             kmin = max(lo(3),khi+1)
             kmax = hi(3)

             i = ilo-1
             k = khi+1

             if (i.ge.imin .and. i.le.imax .and. k.ge.kmin .and. k.le.kmax) then

                do j = lo(2), hi(2)

                   if (khi-2 >= ks) then
                      q(i,j,k,n) = half * eighth * (15*q(ilo-1,j,khi,n) - 10*q(ilo-1,j,khi-1,n) + 3*q(ilo-1,j,khi-2,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(ilo-1,j,khi,n) - q(ilo-1,j,khi-1,n))
                   end if
                   
                   if (ilo+2 <= ie) then
                      q(i,j,k,n) = q(ilo-1,j,khi+1,n) + &
                           half * eighth * (15*q(ilo,j,khi+1,n) - 10*q(ilo+1,j,khi+1,n) + 3*q(ilo+2,j,khi+1,n))
                   else
                      q(i,j,k,n) = q(ilo-1,j,khi+1,n) + half * half * (3*q(ilo,j,khi+1,n) - q(ilo+1,j,khi+1,n))
                   end if
                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,2,n) .eq. HOEXTRAP .and. bc(3,1,n) .eq. HOEXTRAP) then

          if (hi(1) > ihi .and. lo(3) < klo) then
             imin = max(lo(1),ihi+1)
             imax = hi(1)
             kmin = lo(3)
             kmax = min(hi(3),klo-1)

             i = ihi+1
             k = klo-1

             if (i.ge.imin .and. i.le.imax .and. k.ge.kmin .and. k.le.kmax) then

                do j = lo(2), hi(2)

                   if (klo+2 <= ke) then
                      q(i,j,k,n) = half * eighth * (15*q(ihi+1,j,klo,n) - 10*q(ihi+1,j,klo+1,n) + 3*q(ihi+1,j,klo+2,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(ihi+1,j,klo,n) - q(ihi+1,j,klo+1,n))
                   end if
                   
                   if (ihi-2 >= is) then
                      q(i,j,k,n) = q(ihi+1,j,klo-1,n) + &
                           half * eighth * (15*q(ihi,j,klo-1,n) - 10*q(ihi-1,j,klo-1,n) + 3*q(ihi-2,j,klo-1,n))
                   else
                      q(i,j,k,n) = q(ihi+1,j,klo-1,n) + half * half * (3*q(ihi,j,klo-1,n) - q(ihi-1,j,klo-1,n))
                   end if
                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,2,n) .eq. HOEXTRAP .and. bc(3,2,n) .eq. HOEXTRAP) then
          
          if (hi(1) > ihi .and. hi(3) > khi) then
             imin = max(lo(1),ihi+1)
             imax = hi(1)
             kmin = max(lo(3),khi+1)
             kmax = hi(3)

             i = ihi+1
             k = khi+1

             if (i.ge.imin .and. i.le.imax .and. k.ge.kmin .and. k.le.kmax) then

                do j = lo(2), hi(2)

                   if (khi-2 >= ks) then
                      q(i,j,k,n) = half * eighth * (15*q(ihi+1,j,khi,n) - 10*q(ihi+1,j,khi-1,n) + 3*q(ihi+1,j,khi-2,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(ihi+1,j,khi,n) - q(ihi+1,j,khi-1,n))
                   end if
                   
                   if (ihi-2 >= is) then
                      q(i,j,k,n) = q(ihi+1,j,khi+1,n) + &
                           half * eighth * (15*q(ihi,j,khi+1,n) - 10*q(ihi-1,j,khi+1,n) + 3*q(ihi-2,j,khi+1,n))
                   else
                      q(i,j,k,n) = q(ihi+1,j,khi+1,n) + half * half * (3*q(ihi,j,khi+1,n) - q(ihi-1,j,khi+1,n))
                   end if
                   
                end do
             end if
          end if
       end if

       !
       ! Next correct the j-k edges
       !

       if (bc(2,1,n) .eq. HOEXTRAP .and. bc(3,1,n) .eq. HOEXTRAP) then

          if (lo(2) < jlo .and. lo(3) < klo) then
             jmin = lo(2)
             jmax = min(hi(2),jlo-1)
             kmin = lo(3)
             kmax = min(hi(3),klo-1)

             j = jlo-1
             k = klo-1

             if (j.ge.jmin .and. j.le.jmax .and. k.ge.kmin .and. k.le.kmax) then

                do i = lo(1), hi(1)

                   if (klo+2 <= ke) then
                      q(i,j,k,n) = half * eighth * (15*q(i,jlo-1,klo,n) - 10*q(i,jlo-1,klo+1,n) + 3*q(i,jlo-1,klo+2,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(i,jlo-1,klo,n) - q(i,jlo-1,klo+1,n))
                   end if
                   
                   if (jlo+2 <= je) then
                      q(i,j,k,n) = q(i,jlo-1,klo-1,n) + &
                           half * eighth * (15*q(i,jlo,klo-1,n) - 10*q(i,jlo+1,klo-1,n) + 3*q(i,jlo+2,klo-1,n))
                   else
                      q(i,j,k,n) = q(i,jlo-1,klo-1,n) + half * half * (3*q(i,jlo,klo-1,n) - q(i,jlo+1,klo-1,n))
                   end if

                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(2,1,n) .eq. HOEXTRAP .and. bc(3,2,n) .eq. HOEXTRAP) then

          if (lo(2) < jlo .and. hi(3) > khi) then
             jmin = lo(2)
             jmax = min(hi(2),jlo-1)
             kmin = max(lo(3),khi+1)
             kmax = hi(3)

             j = jlo-1
             k = khi+1

             if (j.ge.jmin .and. j.le.jmax .and. k.ge.kmin .and. k.le.kmax) then

                do i = lo(1), hi(1)

                   if (khi-2 >= ks) then
                      q(i,j,k,n) = half * eighth * (15*q(i,jlo-1,khi,n) - 10*q(i,jlo-1,khi-1,n) + 3*q(i,jlo-1,khi-2,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(i,jlo-1,khi,n) - q(i,jlo-1,khi-1,n))
                   end if
                   
                   if (jlo+2 <= je) then
                      q(i,j,k,n) = q(i,jlo-1,khi+1,n) + &
                           half * eighth * (15*q(i,jlo,khi+1,n) - 10*q(i,jlo+1,khi+1,n) + 3*q(i,jlo+2,khi+1,n))
                   else
                      q(i,j,k,n) = q(i,jlo-1,khi+1,n) + half * half * (3*q(i,jlo,khi+1,n) - q(i,jlo+1,khi+1,n))
                   end if
                   
                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(2,2,n) .eq. HOEXTRAP .and. bc(3,1,n) .eq. HOEXTRAP) then

          if (hi(2) > jhi .and. lo(3) < klo) then
             jmin = max(lo(2),jhi+1)
             jmax = hi(2)
             kmin = lo(3)
             kmax = min(hi(3),klo-1)

             j = jhi+1
             k = klo-1

             if (j.ge.jmin .and. j.le.jmax .and. k.ge.kmin .and. k.le.kmax) then

                do i = lo(1), hi(1)

                   if (klo+2 <= ke) then
                      q(i,j,k,n) = half * eighth * (15*q(i,jhi+1,klo,n) - 10*q(i,jhi+1,klo+1,n) + 3*q(i,jhi+1,klo+2,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(i,jhi+1,klo,n) - q(i,jhi+1,klo+1,n))
                   end if
                   
                   if (jhi-2 >= js) then
                      q(i,j,k,n) = q(i,jhi+1,klo-1,n) + &
                           half * eighth * (15*q(i,jhi,klo-1,n) - 10*q(i,jhi-1,klo-1,n) + 3*q(i,jhi-2,klo-1,n))
                   else
                      q(i,j,k,n) = q(i,jhi+1,klo-1,n) + half * half * (3*q(i,jhi,klo-1,n) - q(i,jhi-1,klo-1,n))
                   end if
                   
                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(2,2,n) .eq. HOEXTRAP .and. bc(3,2,n) .eq. HOEXTRAP) then

          if (hi(2) > jhi .and. hi(3) > khi) then
             jmin = max(lo(2),jhi+1)
             jmax = hi(2)
             kmin = max(lo(3),khi+1)
             kmax = hi(3)

             j = jhi+1
             k = khi+1

             if (j.ge.jmin .and. j.le.jmax .and. k.ge.kmin .and. k.le.kmax) then

                do i = lo(1), hi(1)

                   if (khi-2 >= ks) then
                      q(i,j,k,n) = half * eighth * (15*q(i,jhi+1,khi,n) - 10*q(i,jhi+1,khi-1,n) + 3*q(i,jhi+1,khi-2,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(i,jhi+1,khi,n) - q(i,jhi+1,khi-1,n))
                   end if
                   
                   if (jhi-2 >= js) then
                      q(i,j,k,n) = q(i,jhi+1,khi+1,n) + &
                           half * eighth * (15*q(i,jhi,khi+1,n) - 10*q(i,jhi-1,khi+1,n) + 3*q(i,jhi-2,khi+1,n))
                   else
                      q(i,j,k,n) = q(i,jhi+1,khi+1,n) + half * half * (3*q(i,jhi,khi+1,n) - q(i,jhi-1,khi+1,n))
                   end if
                   
                end do
             end if
          end if
       end if
#endif

    end do

  end subroutine my_filccn

  subroutine my_amrex_hoextraptocc (q, qlo, qhi, domlo, domhi, dx, xlo) &
       bind(c,name='my_amrex_hoextraptocc')
    integer, intent(in) :: qlo(3), qhi(3), domlo(*), domhi(*)
    real(amrex_real), intent(inout) :: q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    real(amrex_real), intent(in) :: dx(*), xlo(*)

#if (AMREX_SPACEDIM == 3)
    call my_hoextraptocc(q,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3),domlo,domhi,dx,xlo)
#elif (AMREX_SPACEDIM == 2)
    call my_hoextraptocc(q,qlo(1),qlo(2),qhi(1),qhi(2),domlo,domhi,dx,xlo)
#endif
  end subroutine my_amrex_hoextraptocc

end module my_amrex_filcc_module
