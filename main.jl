
include("./main/activate_package.jl")
# include("./examples/initialize_run.jl")
GLMakie.activate!()
begin
    seed = 100
    properties = Dict(
        :λ => 2, 
        :A => 5, 
        :B => 0.3, 
        :dt => 0.01, 
        :sigma => 0.0,
        :hamiltonian => 0.0, 
        :dH => 0.0,
        :no_disp_H => 0.0
    )
    number_of_peds = 32
    x_len = 11
    y_len = 5
    model = initialize(number_of_peds,x_len, y_len, leapfrog_step, properties; seed)
    # Agents.allagents(model)



    # GLMakie.activate!()
    params = Dict(
        :λ => 0:0.1:5,
        :A => 0:0.1:10,
        :B => 0:0.1:2,
        :sigma => 0:0.1:2
        # :solar_change => -0.1:0.01:0.1,
    )
    agent_color(a) = mode(a.id,2) == [1,0] ? "#bf2642" : "#2b2b33"
    agent_size(a) = 100
    # agent_marker = :utriangle
    mdata = [:hamiltonian, :dH]
    plotkwargs = (;
        agent_color, agent_size)
    fig, abmobs = Agents.abmexploration(model; dt=1:0.2:5, params, mdata)
    fig
end