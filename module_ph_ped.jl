module ph_ped
using Agents, Random
using CairoMakie
import LinearAlgebra: norm
import Statistics: mean

begin
export Agents, Random, CairoMakie
export norm, mean
export Pedestrian
export del_potU, sum_del_potU, acc
export potU, sum_potU, calc_hamiltonian
export euler_step
export agent_step!, initialize
end

## Define Agent
@agent struct Pedestrian(ContinuousAgent{2,Float64})
    uᵢ::Vector{Float64} # desired_velocity : custom property
    old_pos::Vector{Float64}
    old_vel::Vector{Float64}
end

## helper functions: movement
function del_potU(a1, a2, model, is_old_pos)
    A = abmproperties(model)[:A]
    B = abmproperties(model)[:B]
    # dist = euclidean_distance(a1, a2, model)

    if is_old_pos == true
        use_a1pos = a1.old_pos
        use_a2pos = a2.old_pos
        use_dist = euclidean_distance(a1.old_pos, a2.old_pos, model)
    else
        use_a1pos = a1.pos
        use_a2pos = a2.pos
        use_dist = euclidean_distance(a1.pos, a2.pos, model)
    end

    return (-(use_a1pos-use_a2pos)./use_dist)*A*exp(-use_dist/B)
end

function sum_del_potU(a1::Pedestrian, model, is_old_pos)
    return sum([i.id != a1.id ? del_potU(a1, i, model, is_old_pos) : [0,0] for i in allagents(model)])
end

function acc(a::Pedestrian, model, is_old_pos::Bool = true, is_old_vel::Bool = true)
    if is_old_vel == true
        use_vel = a.old_vel
    else
        use_vel = a.vel
    end
    λ = abmproperties(model)[:λ]
    return λ.*(a.uᵢ - use_vel) - sum_del_potU(a, model, is_old_pos)
end

## helper functions: calculation of hamiltonian
function potU(a1::Pedestrian,a2::Pedestrian, model)
    A = model.A
    B = abmproperties(model)[:B]
    dist = euclidean_distance(a1, a2, model)
    return A*B*exp(-dist/B)
end

function sum_potU(a1::Pedestrian, agents, model)
    return sum([i.id != a1.id ? potU(a1, i, model) : 0 for i in agents])
end

function calc_hamiltonian(model)
    kinetic_energy = 0.5.*sum([transpose(a.vel)*a.vel for a in allagents(model)])
    potential_energy = 0.5.*sum(([sum_potU(a,allagents(model),model) for a in allagents(model)]))
    return kinetic_energy + potential_energy
end

## num methods
function euler_step(agent::Pedestrian, model)
    agent.old_pos = agent.pos
    agent.old_vel = agent.vel
    agent.vel = agent.vel + model.dt*acc(agent, model) # may lead to out of extent errors
    agent.pos = agent.pos + model.dt*agent.old_vel
    return agent.vel, agent.pos
end

function leapfrom_step(agent::Pedestrian, model)
    agent.old_pos = agent.pos
    agent.old_vel = agent.vel
    agent.pos = agent.old_pos + model.dt*agent.old_vel + ((model.dt^2)/2)*acc(agent,model)
    agent.vel = agent.old_vel + ((model.dt)/(2 + model.λ*model.dt))*(acc(agent, model) + acc(agent,model,false,true))
    return agent.vel, agent.pos
end 

# ABM functions
function agent_step!(agent::Pedestrian, model)
    # agent.vel, agent.pos = euler_step(agent,model)
    agent.vel, agent.pos = leapfrom_step(agent,model)
    # move_agent!(agent,model,dt)
    agent.pos = normalize_position(agent.pos,model)
    model.hamiltonian = calc_hamiltonian(model)
end
# Add agents with random initial positions/velocities
function initialize(number_of_peds = 4, x_len = 11, y_len = 5, seed=42; properties = nothing)
    if isnothing(properties)
        properties = Dict(:λ => 2, :A => 5, :B => 0.5, :hamiltonian => 0.0, :dt => dt)
    end
    space = ContinuousSpace((x_len, y_len); periodic = true)
    model = StandardABM(
        Pedestrian,
        space;
        container = Vector,
        agent_step!,
        properties,
        rng = Random.seed!(seed),
        scheduler = Schedulers.Randomly()
    )
    for i in 1:number_of_peds
        add_agent!(model;
            pos = [rand()*x_len, rand()*y_len],  # Initial position
            vel = [0,0], # Initial velocity
            old_pos = [0,0],
            old_vel = [0,0],
            uᵢ = i-1 < (number_of_peds ÷ 2) ? [1,0] : [-1,0]
        )
    end
    return model
end
println("hello")
end

### MAIN ### 

# using .ph_ped

# begin
#     model = initialize(6,11,5,42)
#     for i in Agents.allagents(model)
#         println(i.id, i.pos, i.vel)
#     end

#     data = []
#     mean_data = []
#     for i in 1:5000
#         push!(data, model.hamiltonian)
#         push!(mean_data, mean(data))
#         Agents.step!(model)
#         if mod(i,100) == 0
#             print("$i ")
#         end
#     end
#     println()
#     println(model.hamiltonian)
#     println(maximum(data))
#     println(mean(data))
# end

# begin
#     f = CairoMakie.Figure()
#     ax = CairoMakie.Axis(f[1, 1], xlabel = "time_steps", ylabel = "H",
#         title = "")
#     CairoMakie.lines!(ax,data)
#     CairoMakie.lines!(ax,mean_data)
#     f
# end