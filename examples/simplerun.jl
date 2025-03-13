using LinearAlgebra
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
    name = "uni"
    # Agents.allagents(model)
end

CairoMakie.activate!()
model = initialize(number_of_peds,x_len, y_len, leapfrog_step, properties; seed)
mdata = [:hamiltonian, :dH, :dist12]
adata = [:pos]
adf, mdf = @time Agents.run!(model,10000;mdata,adata)

begin
fig = Figure()
ax = Axis(fig[1, 1], limits = (0,11,0,5), xlabel = "x_pos", ylabel = "y_pos")
# for i in Agents.allagents(model)
#     CairoMakie.text!(ax, text="$(i.id)", adf[adf.id .== i.id, :].pos[end])
# end
# # for i in Agents.allagents(model)
#     CairoMakie.text!(ax, text="$(i.id)", adf[adf.id .== i.id, :].pos[1])
# end

# for i in [i for i in Agents.allagents(model) if mod(i.id,2) != 0]
for i in [model[i] for i in 1:2:3]
    CairoMakie.scatter!(ax, adf[adf.id .== i.id, :].pos[9000:end], label="$(i.id)", markersize=2)
end
# CairoMakie.axislegend()
fig
end
# CairoMakie.lines(mdf[2:end,:dist12])
# CairoMakie.lines(diff(mdf[2:end,:dist12]))
begin
fig = Figure()

ax = Axis(fig[1, 1], xlabel = "timesteps [dt]", ylabel = "distance", title="Counter flow | Distances")
for agent in [i for i in Agents.allagents(model) if mod(i.id,2) == 1]
    p = 1
    if agent.id != p && agent.id != 21
        distance_y = [i[2] for i in adf[adf.id .== p, 2:end].pos] - [i[2] for i in adf[adf.id .== agent.id, 2:end].pos] 
        CairoMakie.lines!(distance_y, color =:red)
        # CairoMakie.lines!(diff(distance_y), color =:red)
    end
end
fig
# save("Images/distances_$name.png", fig)
end
# fig, ax = CairoMakie.lines(mdf[:, :hamiltonian])
# fig
[i[2] for i in adf[adf.id .==1, :].pos] - [i[2] for i in adf[adf.id .== 3, :].pos]

CairoMakie.lines(mdf[2:end,:dist12])
CairoMakie.lines(diff(mdf[2:end,:dist12]))
