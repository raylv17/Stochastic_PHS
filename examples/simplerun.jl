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
        :alignment => 0.0,
        :stoch_dH => 0.0
    )
    number_of_peds = 32
    x_len = 11
    y_len = 5
    name = "uni"
    # Agents.allagents(model)
end