module plotfile_module

  use amrex_amr_module
  use my_amr_module, only : plot_file, phi_new, t_new, stepno
  
  use amr_data_module

  implicit none

  private

  public :: writeplotfile, write_plotfile_multi2

contains

  subroutine writeplotfile ()

    integer :: nlevs
    character(len=127) :: name
    character(len=16)  :: current_step
    type(amrex_string) :: varname(1)
    
    if      (stepno(0) .lt. 1000000) then
       write(current_step,fmt='(i5.5)') stepno(0)
    else if (stepno(0) .lt. 10000000) then
       write(current_step,fmt='(i6.6)') stepno(0)
    else if (stepno(0) .lt. 100000000) then
       write(current_step,fmt='(i7.7)') stepno(0)
    else if (stepno(0) .lt. 1000000000) then
       write(current_step,fmt='(i8.8)') stepno(0)
    else
       write(current_step,fmt='(i15.15)') stepno(0)
    end if
    name = trim(plot_file) // current_step

    nlevs = amrex_get_numlevels()

    call amrex_string_build(varname(1), "phi")

    call amrex_write_plotfile(name, nlevs, phi_new, varname, amrex_geom, &
         t_new(0), stepno, amrex_ref_ratio)

  end subroutine writeplotfile

subroutine write_plotfile_multi2 ()
    type(amrex_multifab) :: plotmf(0:amrex_max_level)
    type(amrex_string), allocatable :: varname(:)
    integer, dimension(0:amrex_max_level) :: steps, rr
    integer :: ilev, nc
    character(len=127) :: name
    character(len=16)  :: current_step
    integer :: nlevs, icomp, num_comp

    character(len=10), allocatable :: varnamedummy(:)
       
    nc = 11

    allocate(varname(nc))
    allocate(varnamedummy(nc))

    do icomp=1, nc
      write(varnamedummy(icomp),'("phi",I2.2)')icomp
    enddo

    do icomp=1, nc
       call amrex_string_build(varname(icomp), varnamedummy(icomp))
    enddo

    nlevs=amrex_get_finest_level()

      do ilev = 0, nlevs
          call amrex_multifab_build(plotmf(ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, nc, 1)
   	do icomp=1,nc
          call plotmf(ilev) % copy(phi_new(ilev), icomp, icomp, 1, 0)
	enddo
    end do

    if      (stepno(0) .lt. 1000000) then
       write(current_step,fmt='(i5.5)') stepno(0)
    else if (stepno(0) .lt. 10000000) then
       write(current_step,fmt='(i6.6)') stepno(0)
    else if (stepno(0) .lt. 100000000) then
       write(current_step,fmt='(i7.7)') stepno(0)
    else if (stepno(0) .lt. 1000000000) then
       write(current_step,fmt='(i8.8)') stepno(0)
    else
       write(current_step,fmt='(i15.15)') stepno(0)
    end if

    name = trim(plot_file) // current_step

	
    call amrex_write_plotfile(name, nlevs+1, plotmf, varname, amrex_geom, t_new(0), stepno, amrex_ref_ratio)
    !call amrex_write_plotfile(name, nlevs+1, phi_new, varname(1), amrex_geom, t_new(0), stepno, amrex_ref_ratio)

    ! let's not rely on finalizer, which is feature not all compilers support properly.
    do ilev = 0, nlevs
       call amrex_multifab_destroy(plotmf(ilev))
    end do
  end subroutine write_plotfile_multi2


end module plotfile_module
