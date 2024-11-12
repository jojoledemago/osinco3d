
# osinco3d

<table>
  <tr>
    <td>This project is a Fortran90 Computational Fluid Dynamics (CFD) code designed to solve the 3D Navier-Stokes equations for incompressible flow using a high-order finite difference method.</td>
    <td style="max-width: 100px;">
      <img src="https://github.com/jojoledemago/osinco3d/blob/main/examples/tgv_re2500_les/tgv_qcrit_vort_re2500_les.png" alt="TGV Re = 2500"/>
    </td>
  </tr>
</table>

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
- **Poisson equation solver**: Successive Over-Relaxation (SOR) method with a 2nd order finite difference scheme.
- **Boundary conditions**: Periodic or free-slip conditions in all directions, and Dirichlet conditions in the x-direction only.
- **Large Eddy Simulation**: Standard Smagorinsky model.
- **Visualization**:
  - **2D**: Using Gnuplot.
  - **3D**: Using ParaView (XDMF format).
- **Programming language**: Fortran90.
- **Configuration**: Single parameter file `parameters.o3d`.
- **Executable**: The executable file is named `osinco3d.app`.

## Installation

To install the project, clone the repository and compile the source code:

```sh
git clone git@github.com:jojoledemago/osinco3d.git
cd src/
make
```

To clean the compilation directory, use:

```sh
make clean
```

## Usage

Go to the `bin` directory:

```sh
cd bin
```

Edit the `parameters.o3d` file to set the desired simulation parameters.

Run the executable:

```sh
./osinco3d.app
```

## Contributing

I am looking to parallelize the code using the MPI library and to implement a non-reflective outflow condition at \(x = L_x\). If you're interested in contributing, please feel free to help with these tasks.

## Parameters

The `parameters.o3d` file contains various parameters that control the simulation. Below is a description of each parameter group:

### FlowParam
- `typesim`: Type of simulation
    - 0: Read data from `fields.bin` file
    - 1: Convected vortex
    - 2: Taylor-Green vortex
    - 3: Planar Jet
    - 4: Co-planar Jet
    - 5: Mixing Layer
- `l0`: Characteristic size of the flow field
- `u0`: Characteristic velocity of the flow field
- `re`: Reynolds number, defined as \(Re = ␌rac{u_0 \cdot l_0}{
u}\) (where \(
u\) is the kinematic viscosity)
- `ratio`: Ratio between the maximum and minimum velocity in the domain. This parameter depends on the simulation type `typesim`.

### Domain
- `nx`, `ny`, `nz`: Number of discretization points in each direction (\(x\), \(y\), \(z\))
- `xlx`, `yly`, `zlz`: Sizes of the computational domain in each direction
- `x0`, `y0`, `z0`: Origin of the computational domain in each direction

### BoundaryConditions
- `nbcx1`, `nbcxn`: BC for x-direction at faces 1 (\(x(1)\)) and \(nx\) (\(x(nx)\))
- `nbcy1`, `nbcyn`: BC for y-direction at faces 1 (\(y(1)\)) and \(ny\) (\(y(ny)\))
- `nbcz1`, `nbczn`: BC for z-direction at faces 1 (\(z(1)\)) and \(nz\) (\(z(nz)\))
    - 0: Periodic
    - 1: Free Slip
    - 2: Dirichlet (only for x-direction)
- `sim2d`: Set to 1 to execute a 2D simulation.

### AdvanceTime
- `itscheme`: Time integration scheme
    - 1: Euler
    - 2: Second-order Adams-Bashforth
    - 3: Third-order Adams-Bashforth
- `cfl`: Courant-Friedrichs-Lewy (CFL) number, a criterion for numerical stability, defined as \(	ext{CFL} = ␌rac{u_0 \cdot \Delta t}{\Delta x}\)
- `irestart`:
    - 1: Restart simulation from `fields.bin` file (only valid for `typesim = 0`)
    - 0: Start a new simulation
- `itstart`: Index of the first time step
- `itstop`: Index of the last time step

### PoissonEq
- `omega`: Relaxation coefficient for the Successive Over-Relaxation (SOR) method. The theoretical optimum is: \(\omega = ␌rac{2}{1+\sin(\pi \cdot L_x/n_x)}\)
- `idyn`: If 1, the relaxation coefficient \(\omega\) is calculated dynamically at each iteration.
- `eps`: Convergence criterion for the iterative solver
- `kmax`: Maximum number of iterations allowed for the solver

### Scalar (optional)
- `nscr`: Flag to enable scalar equation resolution (set to 1)
- `sc`: Schmidt number, defined as \(Sc = ␌rac{
u}{D}\) (where \(D\) is the scalar diffusivity)

### VisuParameters
- `nfre`: Frequency of saving data for visualization (XDMF format)
- `nsve`: Frequency of saving data for post-processing
- `initstat`: Time after which statistics are collected
- `xpro`, `ypro`, `zpro`: Coordinates of the cut plane for 2D visualization (Gnuplot format)

### InitInflow
- `iin`: Inlet boundary condition type
    - 0: Classical
    - 1: Turbulent
- `inflow_noise`: Turbulence intensity at the inlet boundary (0 to 1), representing a fraction of the characteristic velocity `u0`
- `ici`: Initial condition type
    - 0: Classical
    - 1: Oscillation
    - 2: Noise
    - 3: Both
- `init_noise`: Turbulence intensity for the initial condition (0 to 1), representing a fraction of the characteristic velocity `u0`

### LES (Large Eddy Simulation)
- `iles`: If 1, Large Eddy Simulation is enabled
- `cs`: Smagorinsky constant to evaluate \(
u_{t}\)
