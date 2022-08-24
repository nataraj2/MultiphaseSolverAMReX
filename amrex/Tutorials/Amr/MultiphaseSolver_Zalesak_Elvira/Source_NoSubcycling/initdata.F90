module initdata_module

  use amrex_amr_module

  use my_amr_module, only : restart, plot_int
  use plotfile_module, only : writeplotfile, write_plotfile_multi2
  use averagedown_module, only : averagedown
  use restart_module

  implicit none

  private

  public :: initdata

contains

  subroutine initdata ()
    if (len_trim(restart) .eq. 0) then

       call amrex_init_from_scratch(0.0_amrex_real)

       call averagedown()
       !if (plot_int .gt. 0) call write_plotfile_multi2
       !if (plot_int .gt. 0) call writeplotfile
    else
       call readcheckpointfile ()
       !call amrex_abort("init from checkpoint not implemented yet")
    end if
  end subroutine initdata

end module initdata_module
