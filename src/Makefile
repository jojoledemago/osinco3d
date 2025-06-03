#===============================================================================
# Fortran Compiler and Flags
#===============================================================================
FC      = gfortran
FFLAGS  = -Wall -Wextra -std=f2008 -O3 -cpp -MMD -MP -J$(OBJDIR) -fall-intrinsics -fPIC -I/usr/include
LDFLAGS =
LIBS    = -lfftw3 -lm

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
          $(SRCDIR)/noise_mod.f90 \
          $(SRCDIR)/boundary_conditions.f90 \
          $(SRCDIR)/poisson.f90 \
          $(SRCDIR)/initialization.f90 \
          $(SRCDIR)/IOfunctions.f90 \
          $(SRCDIR)/utils.f90 \
          $(SRCDIR)/visualization.f90 \
          $(SRCDIR)/initial_conditions.f90 \
          $(SRCDIR)/differential_operators.f90 \
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
	@echo "Compilation completed."

# Create directories if they don't exist
directories:
	@mkdir -p $(OBJDIR) $(EXEDIR)

# Link object files into the final executable
$(EXECUTABLE): $(OBJECTS)
	$(FC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBS)

# Compile .f90 files into .o with module output in OBJDIR
$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c -o $@ $<

# Include dependency files for smart rebuilding
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

