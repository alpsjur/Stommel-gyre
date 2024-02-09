using Oceananigans             # Ocean model 🌊
using Oceananigans.Units
using Printf                   # formatting text
using CairoMakie               # plotting

Lx = 2000kilometers
Ly = 2000kilometers
Lz =    4kilometers

Nx = 128
Ny = 128


# Define the grid
grid = RectilinearGrid(
    #CPU();
    GPU();
    size=(Nx, Ny, 1),
    x=(0, Lx),
    y=(0, Ly),
    z=(-Lz,0),
    topology=(Bounded, Bounded, Bounded),
)


# Define the wind stress forcing
τ₀ = 0.1       # Maximum wind stress [Nm⁻²]
ρ = 1025       # Density of seawater [kgm⁻³]
u_surface_stress(x, y, t) = -τ₀ * cos(π * y / Ly) / ρ
u_surface_bc  = FluxBoundaryCondition(u_surface_stress)

# plot forcing
x, y, z = nodes(grid, (Center(), Center(), Center()))
τ = u_surface_stress.(1, y, 1)

fig = Figure()
ax = Axis(fig[1, 1], ylabel = "y [km]", xlabel = "Wind stress [Nm⁻²]")
lines!(ax, τ*ρ, y/1000)
save("figures/surface_forcing.png", fig)

# Define linear bottom drag
r = 1/60days
u_bottom_drag(x, y, t, u) = -r*u
v_bottom_drag(x, y, t, v) = -r*v

u_bottom_bc = FluxBoundaryCondition(u_bottom_drag, field_dependencies=:u)
v_bottom_bc = FluxBoundaryCondition(v_bottom_drag, field_dependencies=:v)

# Define horizontal boundary condition
#horizontal_bc = ValueBoundaryCondition(0.0)  # No-slip boundary condition
horizontal_bc = FluxBoundaryCondition(0.0)   # Free-slip boundary condition

# Collect all boundary conditions
u_bcs = FieldBoundaryConditions(
    top = u_surface_bc,
    bottom = u_bottom_bc,
    north = horizontal_bc,
    south = horizontal_bc,
)
v_bcs = FieldBoundaryConditions(
    bottom = v_bottom_bc,
    east = horizontal_bc,
    west = horizontal_bc,
)


# Define the Coriolis parameter, varying with latitude (simple beta-plane approximation)
f₀ = 1e-4        # [s⁻¹]
β = 2e-11        # [m⁻¹s⁻¹]
coriolis = BetaPlane(f₀=f₀, β=β)


# Initialize the model
model = HydrostaticFreeSurfaceModel(; grid,
                          coriolis = coriolis,
                          boundary_conditions = (u=u_bcs, v=v_bcs),
                          momentum_advection = nothing,
                          #closure = ScalarDiffusivity(ν=2e-4, κ=2e-4),
                          )


# set up simulation
Δt = 20minutes
stop_time = 2*365days
simulation = Simulation(model, Δt=Δt, stop_time=stop_time)  # Δt is the time step in seconds, stop_time is the total simulation time in seconds

# logging simulation progress
start_time = time_ns()
progress(sim) = @printf(
    "i: %d, sim time: % 15s, max(|u|): %.3f ms⁻¹, max(|v|): %.3f ms⁻¹, wall time: %s\n",
    sim.model.clock.iteration,
    prettytime(sim.model.clock.time),
    maximum(abs, u),
    maximum(abs, v),
    prettytime(1e-9 * (time_ns() - start_time))
)

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))


# Add output writers for saving simulation output
save_interval = 1day
filename = "linear_stommel_gyre_output"

u, v, w = model.velocities
η = model.free_surface.η
ζ = ∂x(v) - ∂y(u)
s = sqrt(u^2 + v^2)

# Writer for JLD2 file format
simulation.output_writers[:JLD2] = JLD2OutputWriter(
    model, (; u, v, η, ζ, s
    ),
    schedule = AveragedTimeInterval(save_interval),
    filename = filename,
    overwrite_existing = true,
    with_halos = true,                     # for computation of derivatives at boundaries. Also error for η if this is left out?
)


# Info for NetCDFOutputWriter 
outputs = Dict(
    "u" => u, 
    "v" => v,  
    #"η" => η,     # Oceananigans throws error for writing η. Related to problem with JLD2 when halos not included?
    "ζ" => ζ,
    "s" => s,
    )

output_attributes = Dict(
    "ζ"  => Dict("long_name" => "Relative vorticity", "units" => "1/s"),
    "s"  => Dict("long_name" => "Speed", "units" => "m/s"),
)

# Remove netCDF file if it already exists 
if isfile(filename*".nc")
    rm(filename*".nc")
end

# Writer for netCDF file format
simulation.output_writers[:netCDF] = NetCDFOutputWriter(
    model, outputs,
    output_attributes=output_attributes,
    schedule = AveragedTimeInterval(save_interval),
    filename = filename,
    overwrite_existing = true,
    with_halos = true,                     # for computation of derivatives at boundaries. Also error for η if this is left out?
)



# Run the simulation
run!(simulation)
