
include("./activate_package/activate_package.jl")
# include("./examples/initialize_run.jl")

function model_init()
    begin
        seed = 42
        properties = Dict(
            :λ => 2, :A => 5, :B => 0.3, :dt => 0.01, :sigma => 0.1,
            :no_disp_H => 0.0, :hamiltonian => 0.0, :dH => 0.0, :stoch_dH => 0.0,
            :alignment => 0.0
        )
        number_of_peds = 0
        x_len, y_len = 11, 5
        num_solver = leapfrog_step
        model = initialize(number_of_peds, x_len, y_len, num_solver, properties; seed)
    end

    begin
        rng = Distributions.Xoshiro(seed)
        number_of_peds = 40
        for i in 1:number_of_peds
            Agents.add_agent!(model;
                # pos = Agents.normalize_position(SVector{2}([i,i]),model),
                pos = [rand(rng)*x_len, rand(rng)*y_len],  # Initial position
                vel = [0.0,0.0], # Initial velocity
                # uᵢ = mod(i,4) == 0 ? [1,0] : mod(i,4) == 1 ? [-1,0] : mod(i,4) == 2 ? [0, 1] : [0, -1], # custom
                # uᵢ = mod(i,3) == 0 ? [1,0] : mod(i,3) == 1 ? [-1,0] : [1,1], # custom
                uᵢ = mod(i,2) == 0 ? [1,0] : [-1,0] # counter_flow
                # uᵢ = mod(i,2) == 0 ? [0,1] : [1,0] # cross
                # uᵢ = [1.0,0.0] # uni flow
                # uᵢ = [0,0]
            )
        end
    end
    return model

    begin
        seed = 42
        properties = Dict(
            :λ => 2, :A => 5, :B => 0.3, :dt => 0.01, :sigma => 0.1,
            :no_disp_H => 0.0, :hamiltonian => 0.0, :dH => 0.0, :stoch_dH => 0.0,
            :alignment => 0.0
        )
        number_of_peds = 0
        x_len, y_len = 11, 5
        num_solver = leapfrog_step
        model = initialize(number_of_peds, x_len, y_len, num_solver, properties; seed)
    end

    begin
        rng = Distributions.Xoshiro(seed)
        number_of_peds = 40
        for i in 1:number_of_peds
            Agents.add_agent!(model;
                # pos = Agents.normalize_position(SVector{2}([i,i]),model),
                pos = [rand(rng)*x_len, rand(rng)*y_len],  # Initial position
                vel = [0.0,0.0], # Initial velocity
                # uᵢ = mod(i,4) == 0 ? [1,0] : mod(i,4) == 1 ? [-1,0] : mod(i,4) == 2 ? [0, 1] : [0, -1], # custom
                # uᵢ = mod(i,3) == 0 ? [1,0] : mod(i,3) == 1 ? [-1,0] : [1,1], # custom
                uᵢ = mod(i,2) == 0 ? [1,0] : [-1,0] # counter_flow
                # uᵢ = mod(i,2) == 0 ? [0,1] : [1,0] # cross
                # uᵢ = [1.0,0.0] # uni flow
                # uᵢ = [0,0]
            )
        end
    end
    return model
end

begin
    T = 100
    total_steps = T / model.dt
    mdata = [:hamiltonian, :dH, :no_disp_H, :alignment, :stoch_dH] # max hamiltonian @ step 323
    adata = [:pos, :vel]
    adf, mdf = Agents.run!(model,total_steps;mdata,adata)
end
total_steps - 100
begin
    CairoMakie.activate!()
    t = model.dt:model.dt:T
    h_data, h_star_data = mdf[2:end,:hamiltonian] , mdf[2:end,:no_disp_H]
    fig, ax = lines(t, h_data)
    lines!(ax, t, h_star_data, linestyle=:dash, label=L"H^*")
    ax.ylabel = "H"
    ax.xlabel = "time [s]"
    ax.title = "Hamiltonian for counter flow | n = $number_of_peds"
    axislegend()
    save("Images/example_plot.png", fig)
    fig
end

begin
    path_length = 100
    start_time = Int(T/model.dt) - path_length
    end_time = start_time + path_length
    fig = Figure()
    ax = Axis(fig[1, 1], limits = (0,x_len,0,y_len), xlabel = "x_pos", ylabel = "y_pos")
    for i in 1:number_of_peds
        if mod(i,2) == 0
            scatter!(ax, adf[adf.id .== i, :].pos[start_time:end_time], 
                    color = 0:path_length, colormap = :Reds, markersize=5)
        else
            scatter!(ax, adf[adf.id .== i, :].pos[start_time:end_time], 
                    color = 0:path_length, colormap = :Greens, markersize=5)
        end
    end
    save("Images/example_traj.png", fig)
    fig
end

begin
    GLMakie.activate!()
    model = model_init()
    Agents.step!(model)
    params = Dict(
        :λ => 0:0.1:5,
        :A => 0:0.1:10,
        :B => 0:0.1:2,
        :sigma => 0:0.05:1
        # :solar_change => -0.1:0.01:0.1,
    )
    mdata = [:hamiltonian, :dH, :alignment]
    fig, abmobs = Agents.abmexploration(model; dt=1:0.2:5, params, mdata)
    fig
end