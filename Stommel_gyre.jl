using Oceananigans             # Ocean model 🌊
using Oceananigans.Units
using Printf                   # formatting text
using CairoMakie               # plotting
using CUDA                     # for running on GPU

const Lx = 2000kilometers
const Ly = 2000kilometers
const Lz =    4kilometers

const Nx = 128
const Ny = 128


# Run on GPU (wow, fast!) if available. Else run on CPU
if CUDA.functional()
    architecture = GPU()
    @info "Running on GPU"
else
    architecture = CPU()
    @info "Running on CPU"
end

# Define the grid
grid = RectilinearGrid(
    architecture;
    size=(Nx, Ny, 1),
    x=(0, Lx),
    y=(0, Ly),
    z=(-Lz,0),
    topology=(Bounded, Bounded, Bounded),
)

# Boundary conditions and fluxes
"""
From the documentation of Oceananigans:
The sign convention is such that a positive flux represents the flux of a quantity in the 
positive direction. For example, a positive vertical flux implies a quantity is fluxed 
upwards, in the +z direction.

Due to this convention, a positive flux applied to the top boundary specifies that a 
quantity is fluxed upwards across the top boundary and thus out of the domain. As a result, 
a positive flux applied to a top boundary leads to a reduction of that quantity in the 
interior of the domain; for example, a positive, upwards flux of heat at the top of the 
domain acts to cool the interior of the domain. Conversely, a positive flux applied to the 
bottom boundary leads to an increase of the quantity in the interior of the domain. The 
same logic holds for east, west, north, and south boundaries.
"""

# Define the wind stress forcing
const τ₀ = 0.1       # Maximum wind stress [Nm⁻²]
const ρ = 1025       # Density of seawater [kgm⁻³]
const T = 10days     # Timecale for initial increasing surface forcing
u_surface_stress(x, y, t) = -τ₀ * cos(π * y*2 / (Ly)) / ρ * tanh(t/T) # minus sign comes from the sign convention described above

u_surface_bc  = FluxBoundaryCondition(u_surface_stress)   
# plot forcing
figpath = "figures/"

#check if directory figpath exists. Creates it if it does not exist 
if !isdir(figpath)
    mkdir(figpath)
end

x, y, z = nodes(grid, (Center(), Center(), Center()))
τ = u_surface_stress.(1, y, T*5)

fig = Figure()
ax = Axis(fig[1, 1], ylabel = "y [km]", xlabel = "Wind stress [Nm⁻²]")
lines!(ax, -τ*ρ, y/1000)              
save(figpath*"surface_forcing.png", fig)

# Define linear bottom drag
const R = Lz/(60days)
u_bottom_stress(x, y, t, u) = -R*u
v_bottom_stress(x, y, t, v) = -R*v

u_bottom_bc = FluxBoundaryCondition(u_bottom_stress, field_dependencies=:u)
v_bottom_bc = FluxBoundaryCondition(v_bottom_stress, field_dependencies=:v)

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
Δt = 20minutes           # time step 
stop_time = 2*365days      # total simulation time 
simulation = Simulation(model, Δt=Δt, stop_time=stop_time) 

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

simulation.callbacks[:progress] = Callback(progress, IterationInterval(1day/Δt))


# Add output writers for saving simulation output
save_interval = 1day

datapath = "data/"
filename = "linear_double_stommel_gyre_output"
filepath = datapath*filename

#check if directory datapath exists. Creates it if it does not exist 
if !isdir(datapath)
    mkdir(datapath)
end

u, v, w = model.velocities
η = Average(model.free_surface.η, dims=3)
ζ = Field(∂x(v) - ∂y(u))
s = Field(sqrt(u^2 + v^2))

# Writer for JLD2 file format
simulation.output_writers[:JLD2] = JLD2OutputWriter(
    model, (; u, v, η, ζ, s
    ),
    schedule = AveragedTimeInterval(save_interval),
    filename = filepath,
    overwrite_existing = true,
    with_halos = true,                     # for computation of derivatives at boundaries. Also error for η if this is left out?
)


# Info for NetCDFOutputWriter 
outputs = Dict(
    "u" => u, 
    "v" => v,  
    "eta" => η,     # Oceananigans throws error for writing η. Related to problem with JLD2 when halos not included?
    "zeta" => ζ,
    "s" => s,
    )

output_attributes = Dict(
    "zeta"  => Dict("long_name" => "Relative vorticity", "units" => "1/s"),
    "s"  => Dict("long_name" => "Speed", "units" => "m/s"),
    "eta"  => Dict("long_name" => "Sea surface height", "units" => "m"),
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
    filename = filepath,
    overwrite_existing = true,
    with_halos = true,                     # for computation of derivatives at boundaries. Also error for η if this is left out?
)



# Run the simulation
run!(simulation)
