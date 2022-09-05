module solve_helmholtz

 use amrex_amr_module
 use amr_data_module
 use amrex_fi_mpi

  !private
  !public :: solve

contains  

  subroutine solve ()
       call solve_abeclaplacian()
  end subroutine solve
 
  subroutine solve_abeclaplacian ()
    type(amrex_multigrid) :: multigrid
    integer :: ilev, idim
    integer(c_long) :: npts
    real(amrex_real) :: err, avg1, avg2, offset
    logical :: nodal(3)
    type(amrex_multifab) :: null
    integer :: nlevs,do_solve, total_levels

    integer :: myrank, ierr
    type(amrex_boxarray) :: batemp
    type(amrex_distromap) :: dmtemp
    logical :: periodic(3)
    integer :: nc,ng
 
    nlevs=amrex_get_finest_level()

    do_solve=1

    if (composite_solve.and.do_solve.eq.1) then

	!call amrex_geometry_set_periodic([.true., .true., .true.])
			
       call amrex_abeclaplacian_build(abeclap, amrex_geom(0:nlevs), phi_new(0:nlevs)%ba, phi_new(0:nlevs)%dm, &
            metric_term=.false., agglomeration=agglomeration, consolidation=consolidation)
	 
       call abeclap % set_maxorder(linop_maxorder)

       ! This is set up to have homogeneous Neumann BC
       call abeclap % set_domain_bc([amrex_lo_neumann, amrex_lo_neumann, amrex_lo_neumann], &
            &                       [amrex_lo_neumann, amrex_lo_neumann, amrex_lo_neumann])

       !call abeclap % set_domain_bc([amrex_lo_periodic, amrex_lo_periodic, amrex_lo_periodic], &
       !     &                       [amrex_lo_periodic, amrex_lo_periodic, amrex_lo_periodic])

	!do ilev=0,nlevs
	!call Pmfab(ilev)%fill_boundary(amrex_geom(ilev))
	!call rhs(ilev)%fill_boundary(amrex_geom(ilev))
	!enddo

       do ilev = 0, nlevs
          ! for problem with pure homogeneous Neumann BC, we could pass an empty multifab
          !call abeclap % set_level_bc(ilev, null)
          call abeclap % set_level_bc(ilev, Pmfab(ilev))
       end do

       call abeclap % set_scalars(ascalar, bscalar)

        do ilev = 0, nlevs
          call abeclap % set_acoeffs(ilev, acoef(ilev))
          call abeclap % set_bcoeffs(ilev, beta(:,ilev))
       end do

       call amrex_multigrid_build(multigrid, abeclap)
       call multigrid % set_verbose(verbose)
       call multigrid % set_bottom_verbose(cg_verbose)
       call multigrid % set_max_iter(max_iter)
       call multigrid % set_max_fmg_iter(max_fmg_iter)


       err = multigrid % solve(Pmfab, rhs, 1.e-6_amrex_real, 0.0_amrex_real)


       call multigrid%get_grad_solution(gradPmfab)

       call amrex_abeclaplacian_destroy(abeclap)
       call amrex_multigrid_destroy(multigrid)

    end if

  end subroutine solve_abeclaplacian

end module solve_helmholtz
