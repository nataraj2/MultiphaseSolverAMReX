
module bc_module

  use amrex_base_module

  implicit none

  ! periodic bc.  See amrex_bc_types_module for a list of bc types.
  integer, allocatable :: lo_bc(:,:)  ! the second dimension is the
  integer, allocatable :: hi_bc(:,:)  ! number of components

end module bc_module
