begin
    seed = 42
    properties = Dict(
        :λ => 2, 
        :A => 5, 
        :B => 0.3, 
        :dt => 0.1, 
        :sigma => 0.2,
        :hamiltonian => 0.0,
        :dH => 0.0,
        :no_disp_H => 0.0
    )
    number_of_peds = 32
    x_len = 11
    y_len = 5     
end;

begin
    solver_list = [ode_step_imeuler, leapfrog_step, euler_step]
    # dt_list = [0.01, 0.05, 0.1, 0.16, 0.2]
    dt_list = 0.01:0.01:0.2
    n_list = map(x -> 20/x, dt_list)
    all_experiment_time = []
end;
for (solver, name) in zip(solver_list, ["im_euler", "leapfrog", "ex_euler"])
    all_timedata = []
    for (dt, n) in zip(dt_list, n_list)
        properties[:dt] = dt
        model = initialize(number_of_peds, x_len, y_len, solver, properties; seed);
        mdata = [:hamiltonian, :dH, :no_disp_H] # max hamiltonian @ step 323
        adata = [:pos]
        print("running $(String(Symbol(solver))) with $(properties[:dt]) dt | ")
        data = @timed Agents.run!(model,n;mdata,adata)
        timedata = data.time
        push!(all_timedata,timedata)
    end
    push!(all_experiment_time, all_timedata)
end
begin
    CairoMakie.activate!()
    fig = CairoMakie.Figure()
    ax = CairoMakie.Axis(fig[1, 1])
    ax.title = "Comparison of Computational Time | dt:[$(dt_list[1]),$(dt_list[end])]"
    ax.xlabel = "Stepsize [dt]"
    # ax.yticks = 0:11
    ax.ylabel = "Computational Time"
    for (timedata,name) in zip(all_experiment_time,["im_euler", "leapfrog", "ex_euler"])
        lines!(ax, dt_list, timedata, label=name)
    end
    CairoMakie.axislegend(position=:rt)
    save("./Images/time1.png", fig)
end
fig