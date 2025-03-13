

using DiffEqPhysics, DifferentialEquations, Plots
function H(ℒ, θ, params, t)
    g = params[1] #gravitational acceleration
    m = params[2] #mass
    l = params[3] #length

    return ℒ^2/(2*m*l^2) + m*g*l*(1-cos(θ))
end

g = 9.81
m = 2.
l = 1.
params = [g,m,l]
θ₀ = 1.
ℒ₀ = 1.

prob2 = HamiltonianProblem(H, ℒ₀, θ₀, (0., 100.), params)
sol2 = solve(prob2, SofSpa10(), dt = .05);
#==========================================================#


#Plotting
#==========================================================#
# lines(sol1, vars=1, tspan=(0,20), label="ODE Method")
lines(sol2.t, sol2[2,:], label="Hamiltonian Method")
#They produce the same solution!

sol2