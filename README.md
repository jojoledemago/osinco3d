# osinco3d

This project is a Fortran90 Computational Fluid Dynamics (CFD) code designed to solve the 3D Navier-Stokes equations for incompressible flow using a high-order finite difference method.

## Table of Contents 
- [Summary](#summary)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [Parameters](#parameters)

## Summary

- **Prediction-correction method**: Utilizes the Chorin method.
- **Time integration schemes**: Euler, Adams-Bashforth 2nd order (AB2), or Adams-Bashforth 3rd order (AB3).
- **Spatial discretization schemes**: Finite differences.
  - **First derivatives**: 6th order.
  - **Second derivatives**: 4th order.
- **Poisson equation solver**: Successive Over-Relaxation (SOR) method with a 2nd order finite differences scheme.
- **Boundary conditions**: Periodic or free-slip conditions in all directions and Dirichlet conditions in the x-direction only.
- **Visualization**:
  - **2D**: Using Gnuplot.
  - **3D**: Using ParaView (XDMF format).
- **Programming language**: Fortran90.
- **Configuration**: Single parameter file `parameters.o3d`.
- **Executable**: The executable file is named `osinco3d.app`.

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
## Contributing

I would like to parallelize the code using the MPI library and implement a non-reflective outflow condition at $x=L_x$. Therefore, I'm looking for someone to help me with these tasks.

## Parameters
The `parameters.o3d` file contains various parameters that control the simulation. Here's a description of each parameter group:

### FlowParam
- `typesim`: type of simulation
    - 0: read data from `fields.bin` file
    - 1: Convected vortex
    - 2: Taylor-Green vortex
    - 3: Planar Jet
    - 4: Co-planar Jet
    - 5: Mixing Layer
- `l0`: characteristic size of the flow field
- `u0`: characteristic velocity of the flow field
- `re`: Reynolds number of the flow, defined as $Re = \frac{u_0 \cdot l_0}{\nu}$ (nu is the kinematic viscosity)
- `ratio`: ratio between the maximum and minimum velocity in the domain. This parameter depends on the simulation type `typesim`

### Domain
- `nx`, `ny`, `nz`: number of discretisation points  in each direction ($x$, $y$, $z$)
- `xlx`, `yly`, `zlz`: sizes of the computational domain in each direction
- `x0`, `y0`, `z0`: origin of the computational domain in each direction

### BoundaryConditions
- `nbcx1`, `nbcxn`: BC for x-direction at faces 1 $x(1)$ and nx $x(nx)$
- `nbcy1`, `nbcyn`: BC for y-direction at faces 1 $y(1)$ and ny $y(ny)$
- `nbcz1`, `nbczn`: BC for z-direction at faces 1 $z(1)$ and nz $z(nz)$
    - 0: Periodic
    - 1: Free Slip
    - 2: Dirichlet, only for x-direction
- `sim2d`: Set 1 to execute a 2D simulation

### AdvanceTime
- `itscheme`: time integration scheme
    - 1: Euler
    - 2: second-order Adams-Bashforth
    - 3: third-order Adams-Bashforth
- `cfl`: Courant-Friedrichs-Lewy (CFL) number, a criterion for numerical stability, defined as $\text{CFL} = \frac{u_0 \cdot \Delta t}{\Delta x}$
- `irestart`
    - 1: restart simulation from `fields.bin` file (only valid for `typesim = 0`)
    - 0: start a new simulation
- `itstart`: Index of the first time step
- `itstop`: Index of the last time step

### PoissonEq
- `omega`: relaxation coefficient for the Successive Over-Relaxation (SOR) method. The theoretical optimum is: $\omega = \frac{2}{1+\sin(\pi \cdot L_x/n_x)}$
- `eps`: convergence criterion for the iterative solver
- `kmax`: maximum number of iterations allowed for the solver

### Scalar (optional)
- `nscr`: flag to enable scalar equation resolution (set to 1)
- `sc`: Schmidt number, defined as $Sc = \frac{\nu}{D}$ (D is the scalar diffusivity)

### VisuParameters
- `nfre`: frequency of saving data for visualization (XDMF format)
- `nsve`: frequency of saving data for post-processing
- `xpro`, `ypro`, `zpro`: coordinates of the cut plane for 2D visualization (Gnuplot format)

### InitInflow
- `iin`: inlet boundary condition type (0: classical, 1: turbulent)
- `inflow_noise`: turbulence intensity at the inlet boundary (0 to 1), representing a fraction of the characteristic velocity `u0`
- `ici`: initial condition type (0: classical, 1: oscillation, 2: noise, 3: both)
- `init_noise`: turbulence intensity for the initial condition (0 to 1), representing a fraction of the characteristic velocity `u0`
