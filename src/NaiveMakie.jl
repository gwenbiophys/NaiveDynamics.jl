#module NaiveMakie
#using NaiveDynamics

#using GLMakie
#using Colors
#import NaiveDynamics.record_video

# this extension could not have been written without Molly
    # providing a relevant example for this rendering unaware tool
export 
    record_video

function record_video(output_path::AbstractString,
                                    simLog,
                                    collector::Collector,
                                    w::Float64 = 1.05, #window scale factor
                                    frameinterval::Integer = 1
                                    framerate::Integer = 25,
                                    color = :purple,
                                    markersize = 0.05,
                                    line_width = 2.0,
                                    transparency = true,
                                    show_boundary::Bool = true,
                                    boundary_linewidth = 2.0,
                                    boundary_color = :black,
                                    kwargs...
                                    )
                                
    fig = Figure()
    positions::Vector{Vector{MVector{3, Float32}}} = simLog
    #for each in eachindex(simLog)
       # push!(positions, simLog[each].position)
   # end
    positionsToPlot = Observable(Point3f.(positions[1]))

    
    ax = Axis3(fig[1, 1], aspect = :data, title = "graphiti")
    scatter!(ax, positionsToPlot; color=color, markersize=markersize, transparency=transparency,
                markerspace=:data, kwargs...)
    
    xlims!(ax, collector.min_xDim*w, collector.max_xDim*w)
    ylims!(ax, collector.min_yDim*w, collector.max_yDim*w)
    zlims!(ax, collector.min_zDim*w, collector.max_zDim*w)

    frames = 1:length(positions)
    #if frames[frame_i] % frame_i == 0

    record(fig, output_path, frames; framerate) do frame_i

            coord = Point3f.(positions[frame_i])
            positionsToPlot[] = coord

    end

end



#end