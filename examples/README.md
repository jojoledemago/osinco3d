# Osinco3D - Beginner Mini-Guide

This quick guide explains how to compile and run a **Taylor-Green Vortex** simulation with Osinco3D.

---

## 1. Clone the repository

```bash
git clone https://github.com/jojoledemago/osinco3d.git
cd osinco3d/src/
```

---

## 2. Compile the code

```bash
make
```

- To clean and rebuild:

```bash
make clean
make
```

- The executable will be generated in the `bin/` folder under the name `osinco3d.app`.

---

## 3. Configure the Taylor-Green simulation

1. Go to the `bin` folder:

```bash
cd ../bin
```

2. Edit the `parameters.o3d` file:

```fortran
typesim = 2  ! Taylor-Green Vortex
```

- Also check the grid size (`nx`, `ny`, `nz`) and the time step (`dt`) according to your needs.

---

## 4. Run the simulation

```bash
./osinco3d.app
```

- The results will be saved in `fields.bin`.

---

## 5. Visualize the results

- **2D**: Use Gnuplot or Matplotlib.  
- **3D**: Use ParaView to open `.xdmf` or `.vtu` files.

---

## 6. Quick tip

For a quick test with minimal computation:

```fortran
nx = 32
ny = 32
nz = 32
dt = 0.01
```

This allows you to quickly observe the evolution of the Taylor-Green vortex.

---
