module face_velocity_multifab

contains

subroutine build_face_multifab(facemfab)
    use my_amr_module, only : verbose
    use amr_data_module
    implicit none
 
    integer :: idim
    logical :: nodal(3)
    type(amrex_multifab) :: facemfab(amrex_spacedim,0:amrex_max_level)
    integer :: lev, nlevs

    nlevs = amrex_get_finest_level()

    do lev=0,nlevs
       do idim = 1, amrex_spacedim
          nodal = .false.
          nodal(idim) = .true.
          call amrex_multifab_build(facemfab(idim,lev), phi_new(lev)%ba, phi_new(lev)%dm, 1, 0, nodal)
       end do
    end do

end subroutine build_face_multifab

subroutine destroy_face_multifab(facemfab)
    use my_amr_module, only : verbose
    use amr_data_module
    implicit none
 
    integer :: idim
    logical :: nodal(3)
    type(amrex_multifab) :: facemfab(amrex_spacedim,0:amrex_max_level)
    integer :: lev, nlevs

    nlevs = amrex_get_finest_level()

    do lev=0,nlevs
       do idim = 1, amrex_spacedim
          call amrex_multifab_destroy(facemfab(idim,lev))
       end do
    end do

end subroutine destroy_face_multifab


subroutine get_face_velocity_multifab(phi,facemfab)

    use my_amr_module, only : verbose
    use amr_data_module
    use face_velocity_module
    use fillpatch_module, only : fillpatch
    implicit none
  
    real(amrex_real)  :: time

    integer, parameter :: ngrow = 3
    integer :: idim
    logical :: nodal(3)
    type(amrex_multifab) :: phiborder,Uborder
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx, tbx, bxnodal
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pin,pout,pux,puy,puz,pfx,pfy,pfz, &
         pf, pfab, Uin, Gvolptr, Lvolptr
    type(amrex_fab) :: facevel(amrex_spacedim)
    type(amrex_multifab) :: facemfab(amrex_spacedim,0:amrex_max_level),phi(0:)
    integer :: lev, nlevs, check

    time=0.0d0

    nlevs = amrex_get_finest_level()

    do lev=0,nlevs

    call amrex_multifab_build(phiborder, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, ngrow)
    call amrex_multifab_build(Uborder, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, ngrow)

    call fillpatch(lev, time, phiborder,phi,phi)

    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.true.)
    do while(mfi%next())
       bx = mfi%tilebox()

       pin  => phiborder%dataptr(mfi)
       Uin  => Uborder%dataptr(mfi)

       do idim = 1, amrex_spacedim
          tbx = bx
          call tbx%nodalize(idim)
          call tbx%grow(1)
          call facevel(idim)%resize(tbx,1)
       end do

       pux => facevel(1)%dataptr()
       puy => facevel(2)%dataptr()
       puz => facevel(3)%dataptr()

       Gvolptr => Gvolmfab(lev)%dataptr(mfi)
       Lvolptr => Lvolmfab(lev)%dataptr(mfi)

       call get_face_velocity(pin(:,:,:,1),pin(:,:,:,2),pin(:,:,:,3),pin(:,:,:,4),&
                              pin(:,:,:,6),pin(:,:,:,7),pin(:,:,:,8),pin(:,:,:,9),pin(:,:,:,11),&
                              lbound(pin),ubound(pin), &
                              pux, lbound(pux), ubound(pux), &
                              puy, lbound(puy), ubound(puy), &
                              puz, lbound(puz), ubound(puz), &
                              Gvolptr, Lvolptr, lbound(Gvolptr), ubound(Gvolptr),&
                              amrex_geom(lev)%dx, amrex_problo,0.0d0)

         do idim = 1, amrex_spacedim
             pf => facemfab(idim,lev)%dataptr(mfi)
             pfab => facevel(idim)%dataptr()
             tbx = mfi%nodaltilebox(idim)
	     pf       (tbx%lo(1):tbx%hi(1),tbx%lo(2):tbx%hi(2),tbx%lo(3):tbx%hi(3),:) = &
             pfab(tbx%lo(1):tbx%hi(1),tbx%lo(2):tbx%hi(2),tbx%lo(3):tbx%hi(3),:)
          end do

     enddo

     call amrex_mfiter_destroy(mfi)

     do idim = 1, amrex_spacedim
       call amrex_fab_destroy(facevel(idim))
     end do
    !$omp end parallel

    call amrex_multifab_destroy(phiborder)
    call amrex_multifab_destroy(Uborder)

  enddo


end subroutine get_face_velocity_multifab

end module face_velocity_multifab
