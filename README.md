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
## Parameters

### FlowParam
- `typesim`: type of simulation
    - 0: read fields.bin file
    - 1: Convected vortex
    - 2: Taylor-Green vortex
    - 3: Planar Jet
    - 4: Mixing Layer
- `l0`: caracterise size of the flow
- `u0`: caracterise velocity of the flow
- `re`: Reynolds number of the flow $Re = \frac{u_0 \cdot l_0}{\nu}$

### Domain
- `nx`, `ny`, `nz`: number of discretisation points
- `xlx`, `yly`, `zlz`: sizes of the domain
- `x0`, `y0`, `z0`: origin of the thre direction

### BoundaryConditions
- `nbcx1`, `nbcxn`: BC in $x(1)$ and $x(nx)$
- `nbcy1`, `nbcyn`: BC in $y(1)$ and $y(nx)$
- `nbcz1`, `nbczn`: BC in $z(1)$ and $z(nx)$
    - 0: Periodic
    - 1: Free Slip
    - 2: Dirichlet

### AdvanceTime
- `itscheme`: Numerical scheme for time integration
    - 1: Euler
    - 2: Two order Adams-Bashforth
    - 3: Three order Adams-Bashforth
- `cfl`: Courant-Friedrichs-Lewy (CFL) number, a criterion for numerical stability, defined as $\text{CFL} = \frac{u \cdot \Delta t}{\Delta x}$
- `irestart`: 1 for a restart simulation with `typesim` = 0 else 0
- `itstart`: first time index 
- `itstop`: last time index

### Scalar
- `nscr`: 1 for scalar equation resolution 
- `sc`: Schmidt number defined as $Sc = \frac{\nu}{D}$


## Contributing

I would like to parallelize the code using the MPI library. Therefore, I'm looking for someone to help me with this task.

