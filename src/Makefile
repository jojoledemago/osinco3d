#===============================================================================
# Fortran Compiler and Flags
#===============================================================================
FC = gfortran

# Adjust include path if fftw3_f03.mod is elsewhere:
FFTW_INC = /usr/include
FFTW_LIB = /usr/lib

ifdef debug
  FFLAGS = -Wall -Wextra -std=f2008 -O0 -g -cpp -MMD -MP -J$(OBJDIR) \
           -fall-intrinsics -fPIC -fcheck=all -fbacktrace -I$(FFTW_INC) -L$(FFTW_LIB)
  DBGMSG = "*** DEBUG mode enabled ***"
else
  FFLAGS = -Wall -Wextra -std=f2008 -O3 -cpp -MMD -MP -J$(OBJDIR) \
           -fall-intrinsics -fPIC -I$(FFTW_INC) -L$(FFTW_LIB)
  DBGMSG = "--- Release mode (optimized) ---"
endif

LDFLAGS = -L$(FFTW_LIB)
LIBS    = -lfftw3 -lfftw3_threads -lm

#===============================================================================
# Project Structure
#===============================================================================
SRCDIR     ?= ./
OBJDIR     ?= $(SRCDIR)/build
EXEDIR     ?= $(SRCDIR)/../bin
EXECUTABLE = $(EXEDIR)/osinco3d.app

#===============================================================================
# Source Files
#===============================================================================
SOURCES = $(SRCDIR)/functions.f90 \
          $(SRCDIR)/derivation.f90 \
          $(SRCDIR)/noise_module.f90 \
          $(SRCDIR)/poisson.f90 \
          $(SRCDIR)/poisson_multigrid.f90 \
          $(SRCDIR)/initialization.f90 \
          $(SRCDIR)/IOfunctions.f90 \
          $(SRCDIR)/utils.f90 \
          $(SRCDIR)/visualization.f90 \
          $(SRCDIR)/differential_operators.f90 \
          $(SRCDIR)/initial_conditions.f90 \
          $(SRCDIR)/les_turbulence.f90 \
          $(SRCDIR)/integration.f90 \
          $(SRCDIR)/osinco3d_main.f90

#===============================================================================
# Build Rules
#===============================================================================
OBJECTS = $(patsubst $(SRCDIR)/%.f90, $(OBJDIR)/%.o, $(SOURCES))
DEPS    = $(OBJECTS:.o=.d)

#===============================================================================
# Default Target
#===============================================================================
all: directories $(EXECUTABLE)
	@echo $(DBGMSG)
	@echo "Compilation completed."

directories:
	@mkdir -p $(OBJDIR) $(EXEDIR)

$(EXECUTABLE): $(OBJECTS)
	$(FC) $(OBJECTS) -o $@ $(LDFLAGS) $(LIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c -o $@ $<

-include $(DEPS)

#===============================================================================
# Utility Targets
#===============================================================================
clean:
	rm -rf $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/*.d
	@echo "Object and module files cleaned."

mrproper: clean
	rm -rf $(EXECUTABLE)
	@echo "Full clean: executable removed."

.PHONY: all clean mrproper directories

