using DifferentialEquations

function f1!(du, u, p, t)
    du[1] = p[1] * u[1] - p[2] * u[1] * u[2]
end

function f2!(du, u, p, t)
    du[2] = -p[3] * u[2] + p[4] * u[1] * u[2]
end

function g!(du, u, p, t)
    du[2] = p[5] * u[2]
end

u0 = [1.0, 1.0]
tspan = (0.0, 10.0)
p = [1.5, 1.0, 3.0, 1.0, 0.1]  # Include noise parameter


prob = SplitODEProblem(ODEFunction(f1!), SDEFunction(f2!, g!), u0, tspan, p)
alg = CompositeAlgorithm((Euler(), ImplicitEM()))
sol = solve(prob, alg, dt=0.01)


function drift1!(du,u,p,t)
    vel, pos = u[1:2], u[3:4]
    du[1:2] = 2*([1,0] - vel) + 0.1*pos
    du[3:4] = vel
end

function diffu1!(du,u,p,t)
    vel, pos = u[1:2], u[3:4]
    σ = p[1]
    du[1:2] = [σ,σ]
    du[3:4] = [0,0] 
end

u0 = Vector{Float64}([0,0,0,0])
tspan = (0.0, 10.0)
p = [0.1]
W = WienerProcess(0.0, 0.0, 0.0)
prob = SDEProblem(drift1!,diffu1!, u0, tspan, p, noise = W)
sol = solve(prob, ImplicitEM(), dt=0.01)

integrator = DifferentialEquations.init(prob, EM(), dt=0.01)
step!(integrator)
sol
using Plots
plot(sol)
CairoMakie.activate!()
vel1 = [i[3] for i in sol.u]
vel2 = [i[4] for i in sol.u]
CairoMakie.scatter(vel1,vel2)