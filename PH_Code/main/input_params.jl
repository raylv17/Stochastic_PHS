
begin
    number_of_peds = 20
    x_len = 11
    y_len = 5
    seed = 42 # or nothing


    properties::Dict = Dict(
        :λ => 2, 
        :A => 5, 
        :B => 0.3, 
        :hamiltonian => 0.0, 
        :dt => 0.01,
        :sigma => 0.1
    );
end
