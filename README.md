# Stommel Gyre in Oceananigans  
🌊 Stommel gyre model setup in Oceananigans.jl  

The Stommel gyre is a basic model of ocean circulation in a rectangular basin, forced by a wind stress that varies in the north-south direction. This model demonstrates the effects of the Coriolis force and western intensification in ocean currents, making it a fundamental example in geophysical fluid dynamics.

## Theory

The Stommel gyre is governed by the barotropic vorticity equation:

$$\frac{\partial \zeta}{\partial t} + J(\psi, \zeta) = \frac{\beta \psi_x}{H} - \frac{\zeta}{R} + \nabla^2 \zeta + \frac{\nabla \times \tau}{\rho H}$$

where:
- $\zeta = \nabla^2 \psi$: The vorticity.
- $\psi$: The streamfunction, related to the velocity field as $u = -\partial \psi / \partial y$ and $v = \partial \psi / \partial x$.
- $J(\psi, \zeta)$: The Jacobian representing nonlinear advection of vorticity.
- $\beta$: The rate of change of the Coriolis parameter with latitude ($\beta$-effect).
- $R$: A linear drag coefficient.
- $\nabla^2 \zeta$: Viscous dissipation.
- $\nabla \times \tau$: The curl of the wind stress, which provides the vorticity forcing.
- $\rho$: The density of seawater.
- $H$: The depth of the ocean.

The simulation illustrates the formation of a wind-driven gyre with western boundary intensification due to the combined effects of the $\beta$-effect and friction.

## Setting Up the Julia Environment

To run the provided scripts and simulations, you'll need to set up a Julia environment based on the `Project.toml` file in the repository. Follow these steps:


1. **Install Julia**:  
   Download and install Julia from the [official website](https://julialang.org/downloads/).

2. **Clone the repository**:  
   Open a terminal and run:
   ```bash
   git clone https://github.com/alpsjur/Stommel-gyre.git
   cd Stommel-gyre
   ```

3. **Activate the Julia environment**:  
   Launch Julia in the repository directory and activate the project:
   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   ```
   This will install all required dependencies listed in `Project.toml`.

## Running scripts  

   To run a simulation, use Julia to execute the scripts. For example, to run the `Stommel_gyre.jl` script:
   ```bash
   julia --project=. Stommel_gyre.jl
   ```
   This will start the simulation and generate output files based on the script's configuration.


## Creating Animations

You can create animations of the simulation results using the `animate_results.jl` script. This script processes saved data, generates visualizations of fields (e.g., speed and surface elevation), and produces an MP4 animation.

---
Feel free to contribute or raise issues in the repository for suggestions or improvements!
