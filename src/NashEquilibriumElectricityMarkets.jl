module NashEquilibriumElectricityMarkets

    using JuMP
    using BilevelJuMP
    using DataFrames
    using CSV

    include("structures.jl")
    include("input_data.jl")
    include("constraints.jl")

end 
