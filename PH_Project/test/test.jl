using Revise, Pkg
cd("PH_Project/src/test")
Pkg.activate("../../")
##########
using PH_Project
using Agents, Random, DifferentialEquations



x_len, y_len = 11,5
space = ContinuousSpace((x_len, y_len); periodic = true)
properties = Dict(:λ => 2, :A => 5, :B => 0.3, :hamiltonian => 0.0, :dt => 0.01)
rng = Random.seed!(42)
num_solver = euler_step

model = StandardABM(
    Pedestrian,
    space;
    container = Vector,
    agent_step! = (agent, model) -> agent_step!(num_solver, agent::Pedestrian, model),
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

model[3]

allagents(model)
step!(model)
allagents(model)

