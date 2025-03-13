begin
    seed = 42
    properties = Dict(
        :λ => 2,
        :A => 5, 
        :B => 0.3, 
        :dt => 0.01, 
        :sigma => 0.0,
        :hamiltonian => 0.0, 
        :dH => 0.0,
        :no_disp_H => 0.0,
        :dist12 => 0.0,
        :a1 => 1,
        :a2 => 3
    )
    number_of_peds = 32
    x_len = 11
    y_len = 5
    timesteps = 10000
    runs = 5
    name = "cross_sigmaxy"
    # Agents.allagents(model)
end

begin

begin
CairoMakie.activate!()
fig1 = CairoMakie.Figure()
ax1 = Axis(fig1[1,1], xlabel="timesteps [dt]", ylabel="H(t)")
fig2 = CairoMakie.Figure()
ax2 = Axis(fig2[1,1], xlabel="timesteps [dt]", ylabel="dH(t)")
end
## Deterministic Plot
begin
    model = initialize(number_of_peds,x_len, y_len, leapfrog_step, properties; seed)
    model.sigma = 0
    adata = [:pos, :vel]
    mdata = [:hamiltonian, :dH, :no_disp_H]
    adf, mdf = @time Agents.run!(model,timesteps;mdata,adata)
    CairoMakie.lines!(ax1, mdf[2:end, :no_disp_H], linestyle=:dash, label=L"H^*")
    CairoMakie.lines!(ax1, mdf[2:end,:hamiltonian],color=:black, label="σ = 0")
    CairoMakie.lines!(ax2, mdf[2:end,:dH],color=:black, label="σ = 0")
end
    
## Sigma 0.2 # 3 runs
begin
for (sig, clr) in zip([0.05,0.1,0.2,0.5,1], [:red, :green, :blue, :orange, :magenta])
model.sigma = sig
all_dh = []
all_ham = []
    for i in 1:runs
        model = initialize(number_of_peds,x_len, y_len, stochastic_ode_step, properties; seed)
        adata = [:pos]
        mdata = [:hamiltonian, :dH]
        print("$i - $sig |")
        adf, mdf = @time Agents.run!(model,timesteps;mdata,adata)
        push!(all_ham, mdf[:,:hamiltonian])
        push!(all_dh, mdf[:,:dH])
        CairoMakie.lines!(ax1, mdf[2:end,:hamiltonian]  , color = clr, alpha=0.2)
        CairoMakie.lines!(ax2, mdf[2:end,:dH]           , color = clr, alpha=0.2)
        fig1
        fig2
    end
fig1
mean_list = []
for i in 1:length(all_ham[1])
    el = 0
    summ = 0
    for j in all_ham
        summ = summ + j[i]
    end
    push!(mean_list, summ/length(all_ham))
end
CairoMakie.lines!(ax1, mean_list[2:end], color = clr, label="σ = $sig")

mean_list = []
for i in 1:length(all_dh[1])
    el = 0
    summ = 0
    for j in all_dh
        summ = summ + j[i]
    end
    push!(mean_list, summ/length(all_ham))
end
CairoMakie.lines!(ax2, mean_list[2:end], color = clr, label="σ = $sig")
end
end
# CairoMakie.axislegend(fig1[1,2], position=:rb)
# CairoMakie.axislegend(fig2[1,2], position=:rt)
# fig1[1,2] = CairoMakie.Legend(fig1, ax1)
# fig2[1,2] = CairoMakie.Legend(fig2, ax2)
CairoMakie.axislegend(ax1)
CairoMakie.axislegend(ax2)
fig1
fig2
end
fig2

begin
    model = initialize(number_of_peds,x_len, y_len, leapfrog_step, properties; seed)
    model.sigma = 0
    adata = [:pos]
    mdata = [:hamiltonian, :dH, :no_disp_H]
    adf, mdf = @time Agents.run!(model,timesteps;mdata,adata)
    CairoMakie.lines!(ax1, mdf[2:end, :no_disp_H], linestyle=:dash, label=L"H^*")
    CairoMakie.lines!(ax1, mdf[2:end,:hamiltonian],color=:black, label="σ = 0")
    CairoMakie.lines!(ax2, mdf[2:end,:dH],color=:black, label="σ = 0")
end

save("Images/H_stochasic_$name.png", fig1)
save("Images/dH_stochasic_$name.png", fig2)
fig = Figure(size=[1000,200])
ax = Axis(fig[1,1])
x = range(1,4*π,100)
CairoMakie.lines(x,sin.(x))