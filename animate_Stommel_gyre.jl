using Oceananigans
using Oceananigans.Units
using Printf
using CairoMakie

filename = "stommel_gyre_output.jld2"
figpath = "figures/"

#check if directory figpath exists. Creates it if it does not exist 
if !isdir(figpath)
    mkdir(figpath)
end

s_timeseries = FieldTimeSeries(filename, "s")
u_timeseries = FieldTimeSeries(filename, "u")
v_timeseries = FieldTimeSeries(filename, "v")
η_timeseries = FieldTimeSeries(filename, "η")

times = s_timeseries.times

xf, yc, zc = nodes(s_timeseries[1])
xc, yc, zf = nodes(η_timeseries[1])

# Initialize logging for the animation creation process
@info "Making an animation from saved data..."

# Create an observable integer for indexing through time series data
n = Observable(1)

# Define a title for the animation with the time variable dynamically updated
title = @lift @sprintf("%s", prettytime(times[$n]))

# Extract the interior data for the fields at the current time step, dynamically updated
# +1 means that initial time-step is skipped. Needed for contour plot
sₙ = @lift interior(s_timeseries[$n+1], :, :, 1)
uₙ = @lift interior(u_timeseries[$n+1], :, :, 1)
vₙ = @lift interior(v_timeseries[$n+1], :, :, 1)
ηₙ = @lift interior(η_timeseries[$n+1], :, :, 1)


# Set limits for the velocity color scale
slim = maximum(interior(s_timeseries))
ηlim = maximum(abs, interior(η_timeseries))

# Define common axis keywords for both plots
axis_kwargs = (xlabel = "x [km]",
               ylabel = "y [km]",
               aspect = 1,
               titlesize = 20)

# Create a figure object for the animation
fig = Figure(size = (1200, 500))

# Create axes for speed
ax_s = Axis(fig[2, 1]; title = "speed [m/s]", axis_kwargs...)
ax_η = Axis(fig[2, 3]; title = "η [m]", axis_kwargs...)

# Add a title
fig[1, :] = Label(fig, title, fontsize=24, tellwidth=false)

# Create a heatmap for the speed
hm_s = heatmap!(ax_s, xf*1e-3, yc*1e-3, sₙ; colorrange = (0,slim), colormap = :speed)
Colorbar(fig[2, 2], hm_s)

# Create a heatmap for the speed
hm_η = heatmap!(ax_η, xc*1e-3, yc*1e-3, ηₙ; colorrange = (-ηlim,ηlim), colormap = :balance)
c_η = contour!(ax_η, xc*1e-3, yc*1e-3, ηₙ; color = "black")
Colorbar(fig[2, 4], hm_η)

# Define the frame range for the animation
frames = 1:length(times)-1

# Record the animation, updating the figure for each time step
record(fig, figpath*"animated_Stommel_gyre.mp4", frames, framerate=8) do i
    msg = string("Plotting frame ", i, " of ", frames[end])
    print(msg * " \r")  # Log progress without creating a new line for each frame
    n[] = i             # Update the observable to the current frame index
end

