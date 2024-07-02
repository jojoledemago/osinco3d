# osinco3d

This project is a Fortran90 Computational Fluid Dynamics (CFD) code designed to solve the 3D Navier-Stokes equations for incompressible flow using a high-order finite difference method.

## Table of Contents 
- [Installation](#installation)
- [Usage](#usage)
- [Parameters](#parameters)
- [Contributing](#contributing)

## Installation

Clone the repository:
```sh
git clone git@github.com:jojoledemago/osinco3d.git
cd src/
make
```
Use make `clean` to clean compilation directory

## Usage

Got to the `bin` directory:
```sh
cd bin
```

Edit the paramters file: `parameters.o3d`

```sh
./osinco3d.app
```

## Paramètres

### FlowParam

```fortran
&FlowParam
typesim = 2  ! Type de simulation (0=lecture d'un fichier fields.bin, 1=vortex convecté, 2=vortex de Taylor-Green, 3=jet plan, 4=couche de mélange)
l0 = 1.       ! Longueur caractéristique de l'écoulement
u0 = 1.       ! Vitesse caractéristique de l'écoulement
re = 1000.    ! Nombre de Reynolds
ratio = 1.00  ! Ratio (à préciser)
/End
```


## Contributing

I would like to parallelize the code using the MPI library. Therefore, I'm looking for someone to help me with this task.

