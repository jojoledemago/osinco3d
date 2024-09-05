# Compilateur Fortran
FC = gfortran

# Options de compilation
FFLAGS = -Wall -Wextra -std=f2008 -O3 -cpp -MMD -MP -J$(OBJDIR) -fall-intrinsics

# Dossiers de sortie
SRCDIR = ./
OBJDIR = $(SRCDIR)/build
EXEDIR = $(SRCDIR)/../bin

# Noms des fichiers sources
SOURCES = $(SRCDIR)/functions.f90 \
          $(SRCDIR)/derivation.f90 \
          $(SRCDIR)/boundary_conditions.f90 \
          $(SRCDIR)/poisson.f90 \
          $(SRCDIR)/initialization.f90 \
          $(SRCDIR)/IOfunctions.f90 \
          $(SRCDIR)/utils.f90 \
          $(SRCDIR)/visualization.f90 \
          $(SRCDIR)/initial_conditions.f90 \
          $(SRCDIR)/differential_operators.f90 \
          $(SRCDIR)/integration.f90 \
          $(SRCDIR)/osinco3d.f90
# Nom de l'executable
EXECUTABLE = $(EXEDIR)/osinco3d.app

# Liste des fichiers objets à générer
OBJECTS = $(patsubst $(SRCDIR)/%.f90,$(OBJDIR)/%.o,$(SOURCES))
DEPS = $(OBJECTS:.o=.d)

# Règle par défaut pour compiler l'exécutable
all: directories $(EXECUTABLE)
	@echo "Compilation completed"

# Règle pour créer les répertoires de sortie
directories:
	@mkdir -p $(OBJDIR) $(EXEDIR)

# Règle pour créer l'exécutable
$(EXECUTABLE): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $^

# Règle pour créer les fichiers objets
$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c -o $@ $<

# Inclusion des fichiers de dépendances
-include $(DEPS)

# Règle pour nettoyer les fichiers temporaires et l'exécutable
clean:
	rm -rf $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/*.d 
