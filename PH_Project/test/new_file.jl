using Pkg, Revise
using Agents, Random, CairoMakie, DifferentialEquations


### Agents Side
@agent struct Pedestrian(ContinuousAgent{2,Float64})
    uᵢ::Vector{Float64} # desired_velocity : custom property
    # old_pos::Vector{Float64}
    # old_vel::Vector{Float64}
end


# Example move!
function simple_move(agent::Pedestrian, model)
    p = agent.pos + model.dt.*[1,0]
    v = agent.vel + model.dt.*rand(2)
    return v, p
end

function agent_step!(agent::Pedestrian, model, num_solver)
    # agent.vel, agent.pos = euler_step(agent,model)
    # agent.vel, agent.pos = leapfrog_step(agent,model)
    agent.vel, p = num_solver(agent,model)
    println(agent.vel, p)
    # move_agent!(agent,model,dt)
    # print("$(agent.id):: $p")
    p = normalize_position(p,model)
    # println(", $p")
    move_agent!(agent, p, model)
end


function initialize(number_of_peds, x_len, y_len, num_solver)
    space = ContinuousSpace((x_len, y_len); periodic = true)
    properties = Dict(:λ => 2, :A => 5, :B => 0.3, :dt => 0.01)
    rng = Random.seed!(42)
    # num_solver = euler_step

    model = StandardABM(
        Pedestrian,
        space;
        container = Vector,
        agent_step! = (agent, model) -> agent_step!(agent::Pedestrian, model, num_solver),
        properties,
        rng,
        scheduler = Schedulers.Randomly()
    )

    for i in 1:number_of_peds
        add_agent!(model;
            # pos = [rand()*x_len, rand()*y_len],  # Initial position
            pos = [i,0],
            vel = [0,0], # Initial velocity
            # uᵢ = i-1 < (number_of_peds ÷ 2) ? [1,0] : [-1,0]
            uᵢ = [1,0]
        )
    end

    return model
end

model = initialize(4,11,5, simple_move)
step!(model)
allagents(model)
### DiffEq Side

function ph_ped_model!(du,u,p,t)
    agent, model = p
    vel, pos = u[1:2], u[3:4]
    dist = (a1,a2) -> euclidean_distance(a1,a2,model)
    potU = (pos1,pos2) -> - ((pos1 - pos2)/dist(pos1,pos2))*model.A*exp(-dist(pos1,pos2)/model.B)
    du[1:2] = model.λ*(agent.uᵢ - vel) - sum( [(i.id != agent.id) ? potU(pos, i.pos) : [0,0] for i in allagents(model)] )
    du[3:4] = vel
end

function ode_step(agent, model)
    
    u0 = [agent.vel[1], agent.vel[2], agent.pos[1], agent.pos[2]]
    # u0 = [agent.vel, agent.pos]
    tspan = (0.0, Inf)
    p = (agent, model)
    prob = ODEProblem{true}(ph_ped_model!, u0, tspan, p)
    
    # Use the model's integrator
    integrator = DifferentialEquations.init(prob, Euler(), dt=model.dt)
    
    # Solve and update agent state
    step!(integrator, model.dt)
    vel = integrator.u[1:2]
    pos = integrator.u[3:4]
    # pos = normalize_position(pos,model)
    return vel, SVector{2}(pos)
end

model = initialize(4,11,5, ode_step)
allagents(model)
step!(model)
allagents(model)

# agent = model[1]
# u0 = Vector{Float64}()
# push!(u0,agent.vel[1])
# push!(u0,agent.vel[2])
# push!(u0,agent.pos[1])
# push!(u0,agent.pos[2])
# println(u0)
# tspan = (0.0, Inf)
# p = [agent, model]
# prob = ODEProblem{true}(ph_ped_model!, u0, tspan, p)

# # Use the model's integrator
# integrator = DifferentialEquations.init(prob, Euler(), dt=model.dt)

# # Solve and update agent state
# step!(integrator, model.dt)