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
        :no_disp_H => 0.0,
        :dist12 => 0.0,
        :a1 => 1,
        :a2 => 1
    )
    number_of_peds = 32
    x_len = 11
    y_len = 5     
end;
model = initialize(number_of_peds, x_len, y_len, leapfrog_step, properties; seed);


begin
    solver_list = [ode_step_exeuler, leapfrog_step, ode_step_imeuler]
    names = ["ex_euler", "leapfrog", "im_euler"]
    # dt_list = [0.01, 0.05, 0.1, 0.16, 0.2]
    dt_list = 0.01:0.01:0.2
    n_list = map(x -> 20/x, dt_list)
    all_experiment_data = []
end;
for (solver, name) in zip(solver_list, names)
    all_diffdata = []
    for (dt, n) in zip(dt_list, n_list)
        properties[:dt] = dt
        model = initialize(number_of_peds, x_len, y_len, solver, properties; seed);
        mdata = [:hamiltonian, :dH, :no_disp_H] # max hamiltonian @ step 323
        adata = [:pos]
        print("running $(String(Symbol(solver))) with $(properties[:dt]) dt | ")
        adf, mdf = @time Agents.run!(model,n;mdata,adata)
        # timedata = data.time
        diffdata = abs(mean(mdf[3:end,:dH] - diff(mdf[2:end,:hamiltonian])*(1/properties[:dt])))
        push!(all_diffdata,diffdata)
    end
    push!(all_experiment_data, all_diffdata)
end
begin
    CairoMakie.activate!()
    fig = CairoMakie.Figure()
    ax = CairoMakie.Axis(fig[1, 1])
    ax.title = "Comparison of Error [dH - (H(t) - H(t-1)] | dt: [$(dt_list[1]),$(dt_list[end])]"
    ax.xlabel = "Stepsize [dt]"
    ax.ylabel = "Error"
    for (diffdata,name) in zip(all_experiment_data, names)
        scatterlines!(ax, dt_list, diffdata, label=name)
    end
    CairoMakie.axislegend(position=:lt)
    save("./Images/leapfrog_compare.png", fig)
end
fig


model = initialize(number_of_peds, x_len, y_len, leapfrog_step, properties; seed)
mdata = [:hamiltonian, :dH, :no_disp_H] # max hamiltonian @ step 323
adata = [:pos]
# print("running $(String(Symbol(solver))) with $(properties[:dt]) dt | ")
adf, mdf = @time Agents.run!(model,10000;mdata,adata)


model.dt
diff(mdf[2:end,:hamiltonian])/model.dt
mdf[2:end, :dH]

mdf[:,:hamiltonian]
mdf[:,:dH]

mean(mdf[2:end, :dH] - diff(mdf[1:end, :hamiltonian])/model.dt)

# dx = [0,2,4,8,16]
# x = [1,2,4,8,16,32]
# diff(x[2:end])
# diff(x[1:end])[2:end]

# dx[2:end] - diff(x[2:end])