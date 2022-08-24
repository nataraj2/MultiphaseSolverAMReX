AMREX_INSTALL_DIR = tmp_install_dir
DIM = 3
USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
USE_ACC = FALSE
COMP = gnu
DEBUG = FALSE
USE_PARTICLES = TRUE
USE_FORTRAN_INTERFACE = TRUE
USE_LINEAR_SOLVERS = TRUE
USE_HYPRE = FALSE
USE_PETSC = FALSE
USE_EB = FALSE
AMREX_XSDK = FALSE
ALLOW_DIFFERENT_COMP = FALSE
USE_SENSEI_INSITU = FALSE
USE_OMP_OFFLOAD = FALSE
TINY_PROFILE = FALSE
USE_COMPILE_PIC = FALSE


AMREX_HOME := $(shell pwd)

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

ifeq ($(USE_CUDA),TRUE)
  ALLOW_DIFFERENT_COMP = TRUE
endif

Pdirs := Base AmrCore Amr Boundary
ifeq ($(USE_PARTICLES),TRUE)
    Pdirs += Particle
   ifeq ($(USE_FORTRAN_INTERFACE),TRUE)
     Pdirs += F_Interfaces/Particle
   endif
endif
ifeq ($(USE_FORTRAN_INTERFACE),TRUE)
    Pdirs += F_Interfaces/Base F_Interfaces/Octree F_Interfaces/AmrCore
endif
ifeq ($(USE_LINEAR_SOLVERS),TRUE)
   Pdirs += LinearSolvers/MLMG LinearSolvers/C_CellMG
   ifeq ($(USE_FORTRAN_INTERFACE),TRUE)
     Pdirs += F_Interfaces/LinearSolvers
   endif
endif
ifeq ($(USE_EB),TRUE)
   Pdirs += EB
endif
ifeq ($(USE_HYPRE),TRUE)
   ifeq ($(USE_LINEAR_SOLVERS),TRUE)
      Pdirs += Extern/HYPRE
   endif
endif
ifeq ($(USE_PETSC),TRUE)
   ifeq ($(USE_LINEAR_SOLVERS),TRUE)
      Pdirs += Extern/PETSc
   endif
endif
ifeq ($(USE_SENSEI_INSITU),TRUE)
	Pdirs += Extern/SENSEI
endif
Ppack := $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)
include $(Ppack)

all: $(amrexlib)
	@echo SUCCESS

.PHONY: distclean install uninstall

distclean: realclean
	$(SILENT) $(RM) GNUmakefile

install: install_lib install_headers install_fortran_modules install_pkg_config

uninstall:
	@echo Uninstalling...
	$(SILENT) $(RM) -r $(AMREX_INSTALL_DIR)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
