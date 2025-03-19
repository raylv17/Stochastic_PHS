begin
    seed = 42
    properties = Dict(
        :λ => 2,
        :A => 5, 
        :B => 0.3, 
        :dt => 0.01, 
        :sigma => 0.1,
        :hamiltonian => 0.0, 
        :dH => 0.0,
        :no_disp_H => 0.0,
        :alignment => 0.0,
        :stoch_dH => 0.0
    )
    number_of_peds = 32
    x_len = 11
    y_len = 5
    timesteps = 10000
    runs = 3
    name = "counter"
    startat = 30
end

@time begin
    begin
        CairoMakie.activate!()
        fig1 = CairoMakie.Figure()
        ax1 = Axis(fig1[1,1], xlabel="timesteps [dt]", ylabel="H(t)")
        fig2 = CairoMakie.Figure()
        ax2 = Axis(fig2[1,1], xlabel="timesteps [dt]", ylabel="dH(t)")
        fig3 = CairoMakie.Figure()
        ax3 = Axis(fig3[1,1], xlabel="timesteps [dt]", ylabel=L"Alignment: $p^T u$")
    end
    ## Deterministic Plot
    # begin
    #     model = initialize(number_of_peds,x_len, y_len, leapfrog_step, properties; seed)
    #     model.sigma = 0
    #     adata = [:pos, :vel]
    #     mdata = [:hamiltonian, :dH, :no_disp_H, :alignment]
    #     adf, mdf = @time Agents.run!(model,timesteps;mdata,adata)
    #     ###################
    #     # std_H = std_range(mdf[:,:hamiltonian], 100)
    #     # std_dH = std_range(mdf[:,:dH], 100)
    #     # std_nH = std_range(mdf[:,:no_disp_H], 100)
    #     # CairoMakie.lines!(ax1, std_nH[1],std_nH[2], linestyle=:dash, label=L"H^*")
    #     # CairoMakie.lines!(ax1, std_H[1],std_H[2],color=:black, label="σ = 0")
    #     # CairoMakie.lines!(ax2, std_dH[1],std_dH[2],color=:black, label="σ = 0")
    #     # print(" $(std(mdf[:,:hamiltonian])), $(std(mdf[:,:dH])) \n")
    #     ###################
    #     CairoMakie.lines!(ax1, mdf[2:end, :no_disp_H], linestyle=:dash, label=L"H^*")
    #     CairoMakie.lines!(ax1, mdf[2:end,:hamiltonian],color=:black, label="σ = 0")
    #     CairoMakie.lines!(ax2, mdf[2:end,:dH],color=:black, label="σ = 0")
    #     CairoMakie.lines!(ax3, mdf[2:end,:alignment],color=:black, label="σ = 0")
    # end
        
    ## Sigma 0.2 # 3 runs
    begin
        # for (sig, clr) in zip([1.0], [:magenta, :orange, :blue, :green, :red])
        for (sig, clr) in zip([1.0,0.5,0.2,0.1,0.05], [:magenta, :orange, :blue, :green, :red])
            properties[:sigma] = sig
            all_dh = []
            all_ham = []
            all_aligns = []
            for i in 1:runs
                model = initialize(number_of_peds,x_len, y_len, stochastic_ode_step, properties; seed)
                adata = [:pos]
                mdata = [:hamiltonian, :dH, :alignment, :stoch_dH]
                print("$i - $sig |")
                adf, mdf = @time Agents.run!(model,timesteps;mdata,adata)
                ################
                # std_H = std_range(mdf[:,:hamiltonian], 100)
                # std_dH = std_range(mdf[:,:dH], 100)
                # push!(all_ham, std_H[2])
                # push!(all_dh, std_dH[2])
                # CairoMakie.lines!(ax1, std_H[1],  std_H[2] , color = clr, alpha=0.2)
                # CairoMakie.lines!(ax2, std_dH[1],std_dH[2] , color = clr, alpha=0.2)
                # print(" $(std(mdf[:,:hamiltonian])), $(std(mdf[:,:dH])) \n")
                ################
                push!(all_ham, mdf[:,:hamiltonian])
                push!(all_dh, mdf[:,:stoch_dH])
                # push!(all_dh, diff(mdf[2:end,:hamiltonian])) # for stochastic dH
                push!(all_aligns, mdf[:,:alignment])
                CairoMakie.lines!(ax1, mdf[2:end,:hamiltonian]  , color = clr, alpha=0.2)
                CairoMakie.lines!(ax2, mdf[1:end,:stoch_dH]           , color = clr, alpha=0.2) 
                # CairoMakie.lines!(ax2, diff(mdf[2:end,:hamiltonian]), color = clr, alpha=0.2) # for stochastic dH
                CairoMakie.lines!(ax3, startat:Agents.abmtime(model)+1, mdf[startat:end,:alignment]         , color = clr, alpha=0.2)
                fig1
                fig2
            end
            mean_list = []
            for i in 1:length(all_ham[1])
                el = 0
                summ = 0
                for j in all_ham
                    summ = summ + j[i]
                end
                push!(mean_list, summ/length(all_ham))
            end
            CairoMakie.lines!(ax1, mean_list[2:end],  color = clr, label="σ = $sig")
            # CairoMakie.lines!(ax1, std_H[1][2:end], mean_list[2:end],  color = clr, label="σ = $sig")

            mean_list = []
            for i in 1:length(all_dh[1])
                el = 0
                summ = 0
                for j in all_dh
                    summ = summ + j[i]
                end
                push!(mean_list, summ/length(all_ham))
            end
            CairoMakie.lines!(ax2, mean_list[1:end], color = clr, label="σ = $sig")

            mean_list = []
            for i in 1:length(all_aligns[1])
                el = 0
                summ = 0
                for j in all_aligns
                    summ = summ + j[i]
                end
                push!(mean_list, summ/length(all_ham))
            end
            CairoMakie.lines!(ax3, startat:Agents.abmtime(model)+1, mean_list[startat:end], color = clr, label="σ = $sig")
            # CairoMakie.lines!(ax2,  std_dH[1][2:end], mean_list[2:end], color = clr, label="σ = $sig")
        end
    end
    # CairoMakie.axislegend(fig1[1,2], position=:rb)
    # CairoMakie.axislegend(fig2[1,2], position=:rt)
    # fig1[1,2] = CairoMakie.Legend(fig1, ax1)
    # fig2[1,2] = CairoMakie.Legend(fig2, ax2)
    begin
        model = initialize(number_of_peds,x_len, y_len, leapfrog_step, properties; seed)
        model.sigma = 0
        adata = [:pos]
        mdata = [:hamiltonian, :dH, :no_disp_H, :alignment, :stoch_dH]
        print("1 - $(model.sigma) |")
        adf, mdf = @time Agents.run!(model,timesteps;mdata,adata)
        ###################
        # std_H = std_range(mdf[:,:hamiltonian], 100)
        # std_dH = std_range(mdf[:,:dH], 100)
        # std_nH = std_range(mdf[:,:no_disp_H], 100)
        # CairoMakie.lines!(ax1, std_nH[1],std_nH[2], linestyle=:dash, label=L"H^*")
        # CairoMakie.lines!(ax1, std_H[1],std_H[2],color=:black, label="σ = 0")
        # CairoMakie.lines!(ax2, std_dH[1],std_dH[2],color=:black, label="σ = 0")
        ###################
        CairoMakie.lines!(ax1, mdf[2:end,:hamiltonian],color=:black, label="σ = 0")
        CairoMakie.lines!(ax2, mdf[2:end,:dH],color=:black, label="σ = 0")
        CairoMakie.lines!(ax3, startat:Agents.abmtime(model)+1, mdf[startat:end,:alignment],color=:black, label="σ = 0")    
        CairoMakie.lines!(ax1, mdf[10:end, :no_disp_H], linestyle=:dash, label=L"H^*")
    end
    CairoMakie.axislegend(ax1, position=:rb) # H
    CairoMakie.axislegend(ax2, position=:rt) # dH
    CairoMakie.axislegend(ax3, position=:rb) # Alignment
    save("Images/stochastic_collective/H_stochasic_$name.png", fig1)
    save("Images/stochastic_collective/dH_stochasic_$name.png", fig2)
    save("Images/stochastic_collective/align_stochastic_$name.png", fig3)
    fig3
end
fig3


fig1
fig2
fig3


# ###############
model = initialize(number_of_peds,x_len, y_len, leapfrog_step, properties; seed)
begin
model.sigma = 0
model.dt = 0.01
adata = [:pos, :vel]
mdata = [:hamiltonian, :dH, :no_disp_H, :stoch_dH, :alignment]
adf, mdf = @time Agents.run!(model,100;mdata,adata)
end;


fig, ax = lines(mdf[20:end, :hamiltonian])
lines!(ax, mdf[1:end, :no_disp_H], linestyle=:dash, label=L"H^*" )
fig, ax = lines(mdf[20:end, :dH])

# lines!(ax, mdf[1:end, :alignment])
# lines!(ax, diff(mdf[2:end,:hamiltonian])/model.dt)
fig

mdf
Agents.abmtime(model)