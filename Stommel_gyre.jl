using Oceananigans             # Ocean model 🌊
using Oceananigans.Units         
using Printf                   # formatting text
using CairoMakie               # plotting

Lx = 1000kilometers
Ly = 1000kilometers 
Lz =    4kilometers

Nx = 128
Ny = 128


# Define the grid
grid = RectilinearGrid(size=(Nx, Ny, 1),  
                       x=(0, Lx),     
                       y=(0, Ly), 
                       z=(-Lz,0),    
                       topology=(Bounded, Bounded, Bounded)
                       )


# Define the wind stress forcing 
τ₀ = 0.1       # Maximum wind stress [Nm⁻²]
ρ = 1025       # Density of seawater [kgm⁻³]
u_surface_stress(x, y, t) = -τ₀ * sin(π * y / Ly) / ρ  
u_surface_bc  = FluxBoundaryCondition(u_surface_stress)

# plot forcing
x, y, z = nodes(grid, (Center(), Center(), Center()))
τ = u_surface_stress.(1, y, 1)

fig = Figure()
ax = Axis(fig[1, 1], ylabel = "y [km]", xlabel = "Wind stress [Nm⁻²]")
lines!(ax, τ*ρ, y/1000)
save("figures/surface_forcing.png", fig)

# Define linear bottom drag
R = 1/60days
u_bottom_drag(x, y, t, u) = -R*u
v_bottom_drag(x, y, t, v) = -R*v

u_bottom_bc = FluxBoundaryCondition(u_bottom_drag, field_dependencies=:u)
v_bottom_bc = FluxBoundaryCondition(v_bottom_drag, field_dependencies=:v)


#
u_bcs = FieldBoundaryConditions(top = u_surface_bc,
                                bottom = u_bottom_bc)
v_bcs = FieldBoundaryConditions(bottom = v_bottom_bc)


# Define the Coriolis parameter, varying with latitude (simple beta-plane approximation)
f₀ = 1e-4        # [s⁻¹]
β = 1e-11        # [m⁻¹s⁻¹]
coriolis = BetaPlane(f₀=f₀, β=β)


# Initialize the model
model = HydrostaticFreeSurfaceModel(; grid,
                          coriolis = coriolis,
                          boundary_conditions = (u=u_bcs, v=v_bcs),
                          )

                        
# set up simulation
Δt = 1hour
stop_time = 365days
simulation = Simulation(model, Δt=Δt, stop_time=stop_time)  # Δt is the time step in seconds, stop_time is the total simulation time in seconds

# logging simulation progress
start_time = time_ns()
progress(sim) = @printf(
    "i: %d, sim time: % 15s, max(|u|): %.3f ms⁻¹, max(|v|): %.3f ms⁻¹, wall time: %s",
    sim.model.clock.iteration,
    prettytime(sim.model.clock.time),
    maximum(abs, u), 
    maximum(abs, v),
    prettytime(wall_time)
)

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))


# Add output writers for saving simulation output
save_interval = 1day

u, v, w = model.velocities
ζ = ∂x(v) - ∂y(u)
s = sqrt(u^2 + v^2)

simulation.output_writers[:timeseries] = JLD2OutputWriter(
    model, (; u, v, ζ, s),
    schedule = AveragedTimeInterval(save_interval),
    filename = "stommel_gyre_output.jld2",
    overwrite_existing = true
)
                                                        

# Run the simulation
run!(simulation)


