push!(LOAD_PATH, "./try1.jl")
using Revise
include("./try1.jl")
using .ph_ped

begin
    T = 10
    dt = 0.01
    println(Int(T/dt))
    properties = Dict(:λ => 2, :A => 5, :B => 0.5, :hamiltonian => 0.0, :dt => dt)
    model = initialize(10,11,5,42; properties)
    for i in Agents.allagents(model)
        println(i.id, i.pos, i.vel)
    end
end

begin
    data = []
    mean_data = []
    for i in 1:Int(T/dt)
        push!(data, model.hamiltonian)
        push!(mean_data, mean(data))
        Agents.step!(model)
        if mod(i,200) == 0
            print("$i ")
        end
    end
    println()
    println(model.hamiltonian)
    println(maximum(data))
    println(mean(data))
end

begin
    using GLMakie
    GLMakie.activate!()
    T = 10
    dt = 0.001
    println(Int(T/dt))
    properties = Dict(:λ => 2, :A => 5, :B => 0.5, :hamiltonian => 0.0, :dt => dt)
    model = initialize(20,10,10,42; properties)
    fig = Figure()
    # ax = CairoMakie.Axis(f[1, 1], xlabel = "time_steps", ylabel = "H",
        # title = "")
    fig, ax, abmobs = Agents.abmplot(model, add_controls = true)
    fig
    # Agents.add_interaction!(ax)
end

begin
    dt = 0.01
    println(Int(T/dt))
    properties = Dict(:λ => 2, :A => 5, :B => 0.5, :hamiltonian => 0.0, :dt => dt)
    model = initialize(5,11,5,42; properties)
    Agents.abmvideo("ph_ped_walk.mp4", model; 
    framerate = 60,
    frames = 1000,
    # dt=10
    )
end


begin
    dt = 0.01
    properties = Dict(:λ => 2, :A => 5, :B => 0.5, :hamiltonian => 0.0, :dt => dt)
    model = initialize(10,11,5,42; properties)
    mdata = [:hamiltonian]
    mdf = Agents.run!(model,1000;mdata)
end


mdf