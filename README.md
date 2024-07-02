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
    . 0: read fields.bin file
    . 1: Convected vortex
    . 2: Taylor-Green vortex
    . 3: Planar Jet
    . 4: Mixing Layer
- `l0`: caracterise size of the flow
- `u0`: caracterise velocity of the flow
- `re`: Reynolds number of the flow [re = \frac{u_0 \cdot l_0}{\nu}]

### Domain
- `nx`, `ny`, `nz`: number of discretisation points
- `xlx`, `yly`, `zlz`: sizes of the domain
- `x0`, `y0`, `z0`: origin of the thre direction

### BoundaryConditions
- `nbcx1`, `nbcxn`: BC in x(1) and x(nx)
- `nbcy1`, `nbcyn`: BC in y(1) and y(nx)
- `nbcz1`, `nbczn`: BC in z(1) and z(nx)

## Contributing

I would like to parallelize the code using the MPI library. Therefore, I'm looking for someone to help me with this task.

