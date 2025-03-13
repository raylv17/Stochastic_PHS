begin
    using Revise, Pkg
    Pkg.activate(".")
    Pkg.develop(path="./PH_Project")
    ##########
    using PH_Project
    using CairoMakie, GLMakie
    using StaticArrays
    using BenchmarkTools
    import Statistics: mean
    import LinearAlgebra: norm
    ############
    # include("./input_params.jl")
    # properties = Dict(:λ => 2, :A => 5, :B => 0.3, :hamiltonian => 0.0, :dt => dt)
end



# seed = 42
# properties = Dict(:λ => 2, :A => 5, :B => 0.3, :hamiltonian => 0.0, :dt => 0.01, :sigma => 0.2)

# begin
#     model = initialize(10,11,5,ode_step, properties; seed)
#     for i in Agents.allagents(model)
#         println(i.id, i.pos, i.vel)
#     end
# end

# begin
#     data = []
#     mean_data = []
#     for i in 1:Int(T/dt)
#         push!(data, model.hamiltonian)
#         push!(mean_data, mean(data))
#         Agents.step!(model)
#         if mod(i,200) == 0
#             print("$i ")
#         end
#     end
#     println()
#     println(model.hamiltonian)
#     println(maximum(data))
#     println(mean(data))
# end

# @time begin
#     dt = 0.01
#     seed = 42
#     # println(Int(T/dt))
#     # properties = Dict(:λ => 2, :A => 5, :B => 0.3, :hamiltonian => 0.0, :dt => dt, :sigma => 0.1)
#     model = initialize(32,11,5, stochastic_ode_step, properties; seed)
#     Agents.abmvideo("stochastic_ode_cross2iem.mp4", model; 
#     framerate = 30,
#     frames = 500,
#     dt=10
#     )
# end


# @time begin
#     CairoMakie.activate!()
#     dt = 0.01
#     seed = 42
#     model = initialize(32,11,5, stochastic_ode_step, properties; seed)
#     mdata = [:hamiltonian] # max hamiltonian @ step 323
#     adata = [:pos]
#     adf, mdf = Agents.run!(model,10000;mdata,adata)
#     # CairoMakie.activate!()
#     a = CairoMakie.lines(mdf[:,:hamiltonian])
#     a
# end
# data
# mean_data = []
# a = CairoMakie.lines(mdf[:,:hamiltonian])
# a
# for i in 1:length(data)
#     push!(mean_data, mean(data[1:i]))
# end
# lines!(mean_data)
# CairoMakie.current_figure()
# # adf, mdf = Agents.run!(model,5000;mdata,adata)

# CairoMakie.lines(mdf[1:3000,:hamiltonian])

# using GLMakie
# GLMakie.activate!()
# fig = Figure()
# begin
#     # T = 10
#     # dt = 0.01
#     # println(Int(T/dt))
#     model = initialize(20,8,2, stochastic_ode_step, properties; seed=42)
#     # agent_color =
#     # plotkwargs = (;
#     # agent_color
#     # )
#     fig, abmobs = Agents.abmexploration(model)
#     fig
#     # Agents.add_interaction!(ax)
# end

# # model = initialize(2,11,5, stochastic_ode_step, properties; seed)
# # Agents.allagents(model)
# # step!(model);
# # Agents.allagents(model)
# model = initialize(0,8,2, stochastic_ode_step, properties; seed=nothing)
# for i in 1:2
#     Agents.add_agent!(model;
#         pos = Agents.normalize_position(SVector{2}([i,i]),model),
#         # pos = [rand()*x_len, rand()*y_len],  # Initial position
#         vel = [0.0,0.0], # Initial velocity
#         # uᵢ = mod(i,4) == 0 ? [1,0] : mod(i,4) == 1 ? [-1,0] : mod(i,4) == 2 ? [0, 1] : [0, -1], # custom
#         # uᵢ = mod(i,3) == 0 ? [1,0] : mod(i,3) == 1 ? [-1,0] : [1,1], # custom
#         # uᵢ = mod(i,2) == 0 ? [1,0] : [-1,0] # counter_flow
#         # uᵢ = mod(i,2) == 0 ? [0,1] : [1,0] # cross
#         uᵢ = [1.0,0.0]
#         # uᵢ = [0,0]
#     )
# end

# using Debugger
# Agents.step!(model)