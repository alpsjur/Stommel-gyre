using Oceananigans             
using Oceananigans.Units         
using Printf                   
using CairoMakie   

filename = "stommel_gyre_output.jld2"

s_timeseries = FieldTimeSeries(filename, "s")
u_timeseries = FieldTimeSeries(filename, "u")
v_timeseries = FieldTimeSeries(filename, "v")

times = s_timeseries.times

xc, yc, zc = nodes(s_timeseries[1]) 


# Initialize logging for the animation creation process
@info "Making an animation from saved data..."

# Create an observable integer for indexing through time series data
n = Observable(1)

# Define a title for the animation with the time variable dynamically updated
title = @lift @sprintf("%s", prettytime(times[$n]))

# Extract the interior data for the fields at the current time step, dynamically updated
sₙ = @lift interior(s_timeseries[$n], :, :, 1)
uₙ = @lift interior(u_timeseries[$n], :, :, 1)
vₙ = @lift interior(v_timeseries[$n], :, :, 1)


# Set limits for the velocity color scale
slim = maximum(interior(s_timeseries))

# Define common axis keywords for both plots
axis_kwargs = (xlabel = "x [km]",
               ylabel = "y [km]",
               aspect = 1,
               titlesize = 20)

# Create a figure object for the animation
fig = Figure()

# Create axes for speed
ax_s = Axis(fig[2, 1]; title = "speed [m/s]", axis_kwargs...)
#ax_u = Axis(fig[3, 1]; title = "streamlines", axis_kwargs...)

# Add a title 
fig[1, :] = Label(fig, title, fontsize=24, tellwidth=false)

# Create a heatmap for the speed
hm_s = heatmap!(ax_s, xc*1e-3, yc*1e-3, sₙ; colorrange = (0,slim), colormap = :speed)
Colorbar(fig[2, 2], hm_s)

# Define the frame range for the animation
frames = 1:length(times)

# Record the animation, updating the figure for each time step
record(fig, "figures/animated_Stommel_gyre.mp4", frames, framerate=8) do i
    msg = string("Plotting frame ", i, " of ", frames[end])
    print(msg * " \r")  # Log progress without creating a new line for each frame
    n[] = i             # Update the observable to the current frame index
end

"""
Testing streamplot, so far not working

fig = Figure()
ax_u = Axis(fig[1, 1]; title = "streamlines", axis_kwargs...)
# create a streamplot for the velocity field
field(i, j) = Point2f(u_timeseries[Int8(i), Int8(j), 1, end], v_timeseries[Int8(i), Int8(j), 1, end])
sp_u = streamplot!(ax_u, field, 1..129, 1..128)
"""