# using Revise, Pkg
# cd("./PH_Project/")
# Pkg.activate(".")
using Agents, Random, DifferentialEquations, DiffEqPhysics
using LinearAlgebra, Statistics, Distributions
## Define Agent
@agent struct Pedestrian(ContinuousAgent{2,Float64})
    uᵢ::Vector{Float64} # desired_velocity : custom property
    # old_pos::Vector{Float64}
    # old_vel::Vector{Float64}
end

## helper functions: movement
# function del_potU(a1, a2, model)
#     A = abmproperties(model)[:A]
#     B = abmproperties(model)[:B]
#     dist = euclidean_distance(a1, a2, model)
#     result = (-(a1.pos-a2.pos)./dist)*A*exp(-dist/B)
#     # println(":: $result")
#     return result
# end

# function sum_del_potU(a1::Pedestrian, model)
#     # checking id is not necessary, since the distance between agents with same id (the same agent) would be zero!
#     # which won't matter in the sum!
#     return sum([i.id != a1.id ? del_potU(a1, i, model) : [0,0] for i in allagents(model)])
# end

# function acc(a::Pedestrian, model)
#     λ = abmproperties(model)[:λ]
#     return λ.*(a.uᵢ - a.vel) - sum_del_potU(a, model)
# end

##

function dU(pos1, pos2, model)
    A = abmproperties(model)[:A]
    B = abmproperties(model)[:B]
    dist = euclidean_distance(pos1, pos2, model)
    result = ((get_direction(pos1,pos2,model))/dist)*A*exp(-dist/B)
    # if pos1 == model[1].pos
    #     println("$pos1 - $pos2 -- dist:$dist")
    # end
    # print(result)
    return result
end

function acc(pos::SVector{2,Float64}, vel::SVector{2,Float64}, agent, model)
    λ = abmproperties(model)[:λ]
    val = λ.*(agent.uᵢ - vel) - sum(dU(pos,i.pos,model) for i in allagents(model) if i.id != agent.id)
    if agent.id == 1
        # println("acc: $val \nterm1: $(λ.*(agent.uᵢ - vel))\nterm2: $(sum([agent.id != i.id ? del_potU2(pos,i.pos,model) : [0 ,0] for i in allagents(model)]))")
    end
    return val
end

## num methods
function euler_step(agent::Pedestrian, model)
    accold1 = acc(agent.pos, agent.vel, agent, model)
    # accold2 = acc(agent,model)
    # println("acc1:$accold1, acc2:$accold2")
    # println("$(agent.id) -- acc: $accold")
    old_vel = agent.vel
    v = agent.vel + model.dt.*accold1
    p = agent.pos + model.dt.*old_vel
    return v, p
end

function leapfrog_step(agent::Pedestrian, model)
    acc_old = acc(agent.pos, agent.vel, agent, model)
    p = agent.pos + model.dt.*agent.vel + ((model.dt^2)/2).*acc_old
    p = normalize_position(p, model)
    acc_new = acc(p, agent.vel, agent, model)
    v = agent.vel + ((model.dt)/(2 + model.λ*model.dt)).*(acc_old + acc_new)
    return v, p
end 

### DiffEq Model

## to be used in ode_step
# function ph_model!(du,u,p,t)
#     agent, model = p
#     vel, pos = u[1:2], u[3:4]
#     dist = (a1,a2) -> euclidean_distance(a1,a2,model)
#     potU = (pos1,pos2) -> - (get_direction(pos1,pos2,model)/dist(pos1,pos2))*model.A*exp(-dist(pos1,pos2)/model.B)
#     du[1:2] = model.λ*(agent.uᵢ - vel) + sum( [(i.id != agent.id) ? potU(pos, i.pos) : [0,0] for i in allagents(model)] )
#     du[3:4] = vel
# end

## to be used in dynamic_ode_step as f1 i.e. p(x)
# function ph_p!(dv,v,u,p,t)
#     agent, model = p
#     vel = v[1:2]
#     pos = u[1:2]
#     dist = (a1,a2) -> euclidean_distance(a1,a2,model)
#     potU = (pos1,pos2) -> - (get_direction(pos1,pos2,model)/dist(pos1,pos2))*model.A*exp(-dist(pos1,pos2)/model.B)
#     dv[1:2] = model.λ*(agent.uᵢ - vel) + sum( [(i.id != agent.id) ? potU(pos, i.pos) : [0,0] for i in allagents(model)] )
# end
##############
function ph_p!(dv,v,u,p,t)
    agent, model = p
    vel = v[1:2]
    pos = u[1:2]
    dv[1:2] = acc(agent.pos, agent.vel, agent, model)
end


function ph_q!(du,v,u,p,t)
    agent, model = p
    vel = v[1:2]
    du[1:2] = vel
end

function dynamic_ode_step(agent, model)
    v0 = [agent.vel[1], agent.vel[2]]
    u0 = [agent.pos[1], agent.pos[2]]
    tspan = (0.0, Inf)
    p = (agent, model)
    prob = DynamicalODEProblem(ph_p!,ph_q!,v0,u0, tspan, p)
    
    # Use the model's integrator
    integrator = DifferentialEquations.init(prob, McAte5(), dt=model.dt)
    
    # Solve and update agent state
    step!(integrator, model.dt)
    vel = integrator.u[1:2]
    pos = integrator.u[3:4]
    # pos = normalize_position(pos,model)
    return vel, SVector{2}(pos)
end

##########

function ph_model!(du,u,p,t)
    agent, model = p
    vel, pos = u[1:2], u[3:4]
    # dist = (a1,a2) -> euclidean_distance(a1,a2,model)
    # potU = (pos1,pos2) -> - (get_direction(pos1,pos2,model)/dist(pos1,pos2))*model.A*exp(-dist(pos1,pos2)/model.B)
    du[1:2] = acc(agent.pos, agent.vel, agent, model)
    du[3:4] = vel
end

function ode_step(agent, model)
    u0 = [agent.vel[1], agent.vel[2], agent.pos[1], agent.pos[2]]
    tspan = (0.0, Inf)
    p = (agent, model)
    prob = ODEProblem{true}(ph_model!, u0, tspan, p)
    
    # Use the model's integrator
    integrator = DifferentialEquations.init(prob, RadauIIA3(), dt=model.dt)
    
    # Solve and update agent state
    step!(integrator, model.dt)
    vel = integrator.u[1:2]
    pos = integrator.u[3:4]
    # pos = normalize_position(pos,model)
    return vel, SVector{2}(pos)
end


function ode_step_exeuler(agent, model)
    u0 = [agent.vel[1], agent.vel[2], agent.pos[1], agent.pos[2]]
    tspan = (0.0, Inf)
    # print(" $(std(mdf[:,:hamiltonian])), $(std(mdf[:,:dH]))")
    p = (agent, model)
    prob = ODEProblem{true}(ph_model!, u0, tspan, p)
    
    # Use the model's integrator
    integrator = init(prob, Euler(), dt=model.dt)
    
    # Solve and update agent state
    step!(integrator, model.dt)
    vel = integrator.u[1:2]
    pos = integrator.u[3:4]
    # pos = normalize_position(pos,model)
    return vel, SVector{2}(pos)
end

function ode_step_imeuler(agent, model)
    u0 = [agent.vel[1], agent.vel[2], agent.pos[1], agent.pos[2]]
    tspan = (0.0, Inf)
    p = (agent, model)
    prob = ODEProblem{true}(ph_model!, u0, tspan, p)
    
    # Use the model's integrator
    integrator = DifferentialEquations.init(prob, ImplicitEuler(), dt=model.dt)
    
    # Solve and update agent state
    step!(integrator, model.dt)
    vel = integrator.u[1:2]
    pos = integrator.u[3:4]
    # pos = normalize_position(pos,model)
    return vel, SVector{2}(pos)
end

#######################

function ham(p,q,params)
    agent, model = params
    vel, pos = p,q
    A = model.A
    B = model.B
    kinetic_energy = 0.5.*sum([transpose(vel)*vel for a in allagents(model)])
    potential_energy = 0.5.*sum([sum_potU(a,allagents(model),model) for a in allagents(model)])
    return kinetic_energy + potential_energy
end

function hamiltonian_ode_step(agent, model)
    v0 = [agent.vel[1], agent.vel[2]]
    u0 = [agent.pos[1], agent.pos[2]]
    tspan = (0.0, Inf)
    p = (agent, model)
    prob = HamiltonianProblem(ham,v0,u0, tspan, p)
    
    # Use the model's integrator
    integrator = init(prob, VerletLeapfrog(), dt=model.dt)
    
    # Solve and update agent state
    step!(integrator, model.dt)
    vel = integrator.u[1:2]
    pos = integrator.u[3:4]
    # pos = normalize_position(pos,model)
    return vel, SVector{2}(pos)
end


## SDE
## explicit/implicit euler_maruyama method
## q -> explicit euler
## p -> implicit euler maruyama
function ph_stoc_drift!(du,u,p,t)
    agent, model = p
    vel, pos = u[1:2], u[3:4]
    # dist = (a1,a2) -> euclidean_distance(a1,a2,model)
    # potU = (pos1,pos2) -> - (get_direction(pos1,pos2,model)/dist(pos1,pos2))*model.A*exp(-dist(pos1,pos2)/model.B)
    # du[1:2] = model.λ*(agent.uᵢ - vel) + sum( [(i.id != agent.id) ? potU(pos, i.pos) : [0,0] for i in allagents(model)] )
    # du[3:4] = vel
    du[1:2] = acc(agent.pos, agent.vel, agent, model)
    du[3:4] = vel
end

function ph_stoc_diffu!(du,u,p,t)
    agent, model = p
    vel, pos = u[1:2], u[3:4]
    σ = model.sigma
    # du[1:2] = σ.*vel
    du[1] = σ
    du[2] = σ
    du[3:4] = [0.0,0.0]
end

function stochastic_ode_step(agent, model)
    u0 = [agent.vel[1], agent.vel[2], agent.pos[1], agent.pos[2]]
    tspan = (0.0, Inf)
    p = (agent, model) 
    W = WienerProcess(0.0, zeros(2), zeros(2))
    prob = SDEProblem(ph_stoc_drift!,ph_stoc_diffu!,u0, tspan, p, noise=W ,noise_rate_prototype=zeros(4, 2))
    
    # Use the model's integrator
    # integrator = DifferentialEquations.init(prob, ImplicitEM(), dt=model.dt)
    # integrator = DifferentialEquations.init(prob, EulerHeun(), dt=model.dt)
    integrator = DifferentialEquations.init(prob, LambaEulerHeun(), dt=model.dt)
    
    # Solve and update agent state
    step!(integrator)
    vel = integrator.u[1:2]
    pos = integrator.u[3:4]
    wiener_increment = integrator.W.dW
    # println(wiener_increment)
    # pos = normalize_position(pos,model)
    return vel, SVector{2}(pos), wiener_increment
end

## Dummy Methods
function simple_move(agent::Pedestrian, model)
    p = agent.pos + model.dt.*[1,0]
    v = agent.vel + model.dt.*rand(2)
    return v, p
end

# ABM functions
function agent_step!(agent::Pedestrian, model, num_solver)
    # agent.vel, agent.pos = euler_step(agent,model)
    # agent.vel, agent.pos = leapfrog_step(agent,model)
    agent.vel, p = num_solver(agent,model)
    # move_agent!(agent,model,dt)
    # print("$(agent.id):: $p")
    p = normalize_position(p,model)
    # println(", $p")
    move_agent!(agent, p, model)
    model.hamiltonian = calc_hamiltonian(model)
end


function model_step!(model,num_solver)
    ## simple_move
    random_list = []
    for agent in allagents(model)
        v, p, rest... = num_solver(agent,model)
        if length(rest) > 0
            push!(random_list, rest[1])
        end
        p = normalize_position(p,model)
        agent.vel = v
        agent.pos = p
    end
    # println(random_list, length(random_list))
    model.hamiltonian = calc_hamiltonian(model)
    model.dH = calc_dH(model)
    model.no_disp_H = calc_no_disp_hamiltonian(model)
    # model.dist12 = euclidean_distance(model[model.a1], model[model.a2], model)
    model.alignment = calc_alignment(model)
    model.stoch_dH = calc_stoch_dH(model, random_list)
end


## helper functions: calculation of hamiltonian
function potU(a1::Pedestrian,a2::Pedestrian, model)
    A = model.A
    B = model.B
    dist = euclidean_distance(a1, a2, model)
    return A*B*exp(-dist/B)
end

function sum_potU(a1::Pedestrian, agents, model)
    return sum([i.id != a1.id ? potU(a1, i, model) : 0 for i in agents])
end

function calc_hamiltonian(model)
    kinetic_energy = 0.5.*sum([transpose(a.vel)*a.vel for a in allagents(model)])
    potential_energy = 0.5.*sum([sum_potU(a,allagents(model),model) for a in allagents(model)])
    return kinetic_energy + potential_energy
end

function calc_no_disp_hamiltonian(model)
    return 0.5*sum(transpose(i.uᵢ) * i.uᵢ for i in allagents(model))
end

function calc_dH(model)
    # return sum([model.λ*norm(i.vel)*norm(i.uᵢ - i.vel) for i in allagents(model)])
    return sum([model.λ*transpose(i.vel)*(i.uᵢ - i.vel) for i in allagents(model)])
end

# ddU
function ddU(a1, a2, model)
    x = get_direction(a1.pos,a2.pos,model) # x
    dist = euclidean_distance(a1,a2,model)
    val_ddU = -model.A*( 
        (1/dist)*(I - (x*transpose(x)/dist^2)) - (1/(model.B*dist^2))*(x*transpose(x)/dist^2)
        ) * exp(-dist/model.B)
    return tr(val_ddU)
end
# trace of the hessian of H(z(t))
function calc_trace_ddH(model)
    total_ddU = 0
    for i in allagents(model)
        for j in allagents(model)
            if i.id != j.id
                total_ddU += ddU(i,j,model)
            end
        end
    end
    return (model.sigma^2/4)*total_ddU + (model.sigma^2/2)*nagents(model)
end

function calc_stoch_dH(model, random_list)
    if length(random_list) == 0
        # dW = [rand(Normal(0, sqrt(model.dt)),2) for _ in 1:nagents(model)] 
        return 0
    else
        dW = random_list # from sdesolver
    end
    # println(dW)
    drift_H = calc_dH(model) + calc_trace_ddH(model)
    # diffu_H = model.sigma * sum(transpose(i.vel) * dW for i in allagents(model))
    diffu_H = model.sigma * sum(transpose(i.vel) * dW[i.id] for i in allagents(model))
    return (drift_H + diffu_H)
end

function calc_alignment(model)
    # return mean([norm(i.uᵢ)/norm(i.vel) for i in allagents(model)])
    return mean([transpose(i.vel/norm(i.vel))*(i.uᵢ) for i in allagents(model)])
end
   


# Add agents with random initial positions/velocities
"""
(number_of_peds = 32, x_len = 11, y_len, num_solver; properties, seed)
"""
function initialize(number_of_peds::Int64 = 32, x_len::Real = 11, y_len::Real = 5, num_solver=euler_step, 
    properties::Union{Nothing,Dict} = nothing; seed::Union{Nothing,Int64} = nothing)
    if isnothing(properties)
        properties = Dict(:λ => 2, :A => 5, :B => 0.3, :hamiltonian => 0.0, :dH => 0.0, :no_disp_H => 0.0, :dt => 0.01, :sigma => 0.1, :alignment => 0.0, :stoch_dH => 0.0)
    end
    rng = Xoshiro(seed)
    space = ContinuousSpace((x_len, y_len); periodic = true, spacing=0.25)
    model = StandardABM(
        Pedestrian,
        space;
        container = Vector,
        # agent_step! = (agent, model) -> agent_step!(agent::Pedestrian, model, num_solver),
        model_step! = (model) -> model_step!(model, num_solver),
        properties,
        rng,
        scheduler = Schedulers.fastest
    )
    for i in 1:number_of_peds
        add_agent!(model;
            # pos = Agents.normalize_position(SVector{2}([i,i]),model),
            pos = [rand(rng)*x_len, rand(rng)*y_len],  # Initial position
            vel = [0.0,0.0], # Initial velocity
            # uᵢ = mod(i,4) == 0 ? [1,0] : mod(i,4) == 1 ? [-1,0] : mod(i,4) == 2 ? [0, 1] : [0, -1], # custom
            # uᵢ = mod(i,3) == 0 ? [1,0] : mod(i,3) == 1 ? [-1,0] : [1,1], # custom
            # uᵢ = mod(i,2) == 0 ? [1,0] : [-1,0] # counter_flow
            # uᵢ = mod(i,2) == 0 ? [0,1] : [1,0] # cross
            # uᵢ = [1.0,0.0] # uni flow
            uᵢ = [0,0]
        )
    end
    return model
end
# println("hello")

