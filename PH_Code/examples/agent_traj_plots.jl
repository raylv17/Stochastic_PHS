begin
    seed = 42
    properties = Dict(
        :λ => 2, 
        :A => 5, 
        :B => 0.3, 
        :dt => 0.01, 
        :sigma => 0.2,
        :hamiltonian => 0.0,
        :dH => 0.0,
        :no_disp_H => 0.0,
        :dist12 => 0.0,
        :a1 => 1,
        :a2 => 1
    )
    number_of_peds = 100
    x_len = 11
    y_len = 5
    T = 15000
    name = "cross_$(number_of_peds)_$(properties[:A])"
end
function run_model(num_solver)
    model = initialize(number_of_peds, x_len, y_len, num_solver, properties; seed)
    @time begin
        mdata = [:hamiltonian, :dH, :no_disp_H] # max hamiltonian @ step 323
        adata = [:pos, :vel]
        adf, mdf = Agents.run!(model,T;mdata,adata)
        data = mdf[:,:dH];
        # fig, ax = CairoMakie.lines(data)
        # mean_data = []
        # for i in 1:length(data)
        #     push!(mean_data, mean(skipmissing(data[1:i])))
        # end
        # CairoMakie.lines!(mean_data)
        # fig
    end
return model, adf, mdf
end
model = initialize(number_of_peds, x_len, y_len, leapfrog_step, properties; seed)
@time begin
    mdata = [:hamiltonian, :dH, :no_disp_H] # max hamiltonian @ step 323
    adata = [:pos, :vel]
    adf, mdf = Agents.run!(model,T;mdata,adata)
    data = mdf[:,:dH];
    # fig, ax = CairoMakie.lines(data)
    # mean_data = []
    # for i in 1:length(data)
    #     push!(mean_data, mean(skipmissing(data[1:i])))
    # end
    # CairoMakie.lines!(mean_data)
    # fig

end;
begin
CairoMakie.activate!()
fig, ax = CairoMakie.lines(mdf[2:end,:hamiltonian])
CairoMakie.lines!(mdf[2:end,:no_disp_H], label = L"H^*", linestyle=:dash)
CairoMakie.axislegend()
ax.ylabel = "H"
ax.xlabel = "time-step [t]"
fig
save("./Images/H_$name.png",fig)
end
begin
    fig, ax = CairoMakie.lines(mdf[2:end,:dH])
    # ax.title = "Hamiltonian "
    ax.ylabel = "dH"
    # ax.xticks = 1:1000:4000
    ax.xlabel = "time-step [t]"
    fig
    save("./Images/dH_$name.png",fig)
end    
begin
    CairoMakie.activate!()
    x = 14900
    path_length = 100
    fig = Figure()
    ax = Axis(fig[1, 1], limits = (0,x_len,0,y_len), xlabel = "x_pos", ylabel = "y_pos")
    # ax = Axis(fig[1, 1], limits = (0,11,0,5), xlabel = "x_pos", ylabel = "y_pos",
        # title = "Agent Trajectories @ t = $(x+path_length) : Counter flow")
        # title = "Agent Trajectories :@t=$(x+path_length) | $name Flow")
    # scatter!(ax, adf[adf.id .== 5, :].pos[1:end], markersize=2 )
    for i in 1:number_of_peds
        if mod(i,4) == 0
            scatter!(ax, adf[adf.id .== i, :].pos[x:x+path_length], color = 0:path_length, colormap = :Reds, markersize=6)
            # scatter!(ax, adf[adf.id .== i, :].pos[2:end], markersize=2, color=:red)
            # scatter!(ax, adf[adf.id .== i, :].pos[1:1], markersize=10, color=:red)
        elseif mod(i,4) == 1
            scatter!(ax, adf[adf.id .== i, :].pos[x:x+path_length], color = 0:path_length, colormap = :Greens, markersize=6)
            # scatter!(ax, adf[adf.id .== i, :].pos[2:end], markersize=2, color=:green)
            # scatter!(ax, adf[adf.id .== i, :].pos[1:1], markersize=10, color=:green)
        elseif mod(i,4) == 2
            scatter!(ax, adf[adf.id .== i, :].pos[x:x+path_length], color = 0:path_length, colormap = :Blues, markersize=6)
        else
            scatter!(ax, adf[adf.id .== i, :].pos[x:x+path_length], color = 0:path_length, colormap = :Purples, markersize=6)
        end
        # scatter!(ax, adf[adf.id .== i, :].pos[x:x+path_length], color = 0:path_length, colormap = :Reds)
        # scatter!(ax, adf[adf.id .== i, :].pos[x:x+path_length], color = 0:path_length, colormap = :Reds, markersize=4)
        # scatter!(ax, adf[adf.id .== i, :].pos[1:end], markersize=2 )
    end
    save("./Images/$(name)flow_$(x+path_length).png", fig)
    fig
end

fig = Figure()
ax = Axis(fig[1, 1])
model1, adf1, mdf1 = run_model(stochastic_ode_step)
model2, adf2, mdf2 = run_model(leapfrog_step)
begin
    CairoMakie.activate!()
    # ax = Axis(fig[1, 1], limits = (0,11,0,5), xlabel = "x_pos", ylabel = "y_pos",
        # title = "Agent Trajectories @ t = $(x+path_length) : Counter flow")
        # title = "Agent Trajectories :@t=$(x+path_length) | $name Flow")
    # scatter!(ax, adf[adf.id .== 5, :].pos[1:end], markersize=2 )
    norm_speed1 = [(mean(adf1[adf1.time .== i, :].vel[:])[1]) for i in 0:T]
    lines!(ax, norm_speed1)
    norm_speed2 = [(mean(adf2[adf2.time .== i, :].vel[:])[1]) for i in 0:T]
    lines!(ax, norm_speed2)
    fig
end
