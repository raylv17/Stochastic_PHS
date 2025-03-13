cd("/home/rafay/Unsynced/Programming_Stuff/julia_works/projects/ph_trys/PH_Pedestrian_Dynamics_org/PH_Project/test/")
using Revise, Pkg
Pkg.activate("../")

##########
using PH_Project
using Agents, Random, DifferentialEquations
using CairoMakie
function f(du,u,p,t)
    du[1] = -u[1]
end

# u0 = [10.0]
# tspan = (0.0, Inf)
# prob = ODEProblem(f, u0, tspan)
# integrator = init(prob, Tsit5())
# integrator
# step!(integrator)
# integrator

@agent struct Pedestrian(ContinuousAgent{2,Float64})
    uᵢ::Vector{Float64} # desired_velocity : custom property
    # old_pos::Vector{Float64}
    # old_vel::Vector{Float64}
end

# plot(sol)
x_len, y_len = 11,5
space = ContinuousSpace((x_len, y_len); periodic = true)
properties = Dict(:λ => 2, :A => 5, :B => 0.3, :hamiltonian => 0.0, :dt => 0.01)
rng = Random.seed!(42)
# num_solver = euler_step

model = StandardABM(
    Pedestrian,
    space;
    container = Vector,
    agent_step! = (agent, model) -> ode_step!(agent::Pedestrian, model),
    properties,
    rng,
    scheduler = Schedulers.Randomly()
)

number_of_peds = 3
for i in 1:number_of_peds
    add_agent!(model;
        # pos = [rand()*x_len, rand()*y_len],  # Initial position
        pos = [i,0],
        vel = [0,0], # Initial velocity
        # old_pos = [0,0],
        # old_vel = [0,0],
        # uᵢ = i-1 < (number_of_peds ÷ 2) ? [1,0] : [-1,0]
        uᵢ = [1,0]
    )
end



allagents(model)

step!(model)

function ph_ped_model(du,u,p,t)
    agent, model = p
    vel, pos = u[1:2], u[3:4]
    dist = (a1,a2) -> euclidean_distance(a1,a2,model)
    potU = (pos1,pos2) -> - ((pos1 - pos2)/dist(pos1,pos2))*model.A*exp(-dist(pos1,pos2)/B)
    du[1:2] = model.λ(agent.uᵢ - vel) - sum( [(i.id != agent.id) ? potU(pos, i.pos) : [0,0] for i in allagents(model)] )
    du[3:4] = vel
end

function ode_step!(agent, model)
    
    u0 = [agent.pos; agent.vel]
    println(u0)
    tspan = (0.0, Inf)
    p = (agent, model)
    prob = ODEProblem{false}(ph_ped_model, u0, tspan, p)
    
    # Use the model's integrator
    integrator = DifferentialEquations.init(prob, Euler(), dt=model.dt)
    
    # Solve and update agent state
    step!(integrator, model.dt)
    pos = integrator.u[1:2]
    agent.vel = integrator.u[3:4]
    pos = normalize_position(pos,model)
    move_agent!(agent,pos,model)
end




u0 = [model[1].vel; model[1].pos]
println(u0)
tspan = (0.0, Inf)
p = (model[1], model)
prob = ODEProblem{false}(ph_ped_model, u0, tspan, p)

# Use the model's integrator
integrator = init(prob, Euler(), dt=model.dt)

# u0 = Vector{Float64}()
# for agent in allagents(model)
#     push!(u0, agent.pos[1])
#     push!(u0, agent.pos[2])
#     push!(u0, agent.vel[1]) #Initial velocity x
#     push!(u0, agent.vel[2]) #Initial velocity y
# end

# prob = ODEProblem(f!,u0,(0.0,10.0), model)
# integrator = init(prob)
# step!(integrator)


# for i in allagents(model)
#     println(i.pos[1])
# end
