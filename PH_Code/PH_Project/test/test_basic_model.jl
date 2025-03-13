cd("/home/rafay/Unsynced/Programming_Stuff/julia_works/projects/ph_trys/PH_Pedestrian_Dynamics_org/PH_Project/test/")
using Revise, Pkg
Pkg.activate("../")

##########
using PH_Project
using Agents, Random, DifferentialEquations
using CairoMakie
##########

begin
    dt = 0.01
    properties = Dict(:λ => 2, :A => 5, :B => 0.3, :hamiltonian => 0.0, :dt => dt)
    model = initialize(20,11,5,leapfrog_step, properties; seed = 42)
end


step!(model)

begin
    mdata = [:hamiltonian]
    adf, mdf = Agents.run!(model,2000;mdata)
    # CairoMakie.activate!()
    CairoMakie.lines(mdf[:,:hamiltonian])
end

begin
    # println(Int(T/dt))
    model = initialize(20,11,5,leapfrog_step, properties; seed = 42)
    abmvideo("simple_mode1.mp4", model;
    frames = 200,
    framerate = 30,
    dt = 10)
end
