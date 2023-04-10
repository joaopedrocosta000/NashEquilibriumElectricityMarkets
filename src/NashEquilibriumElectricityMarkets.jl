module NashEquilibriumElectricityMarkets

    using JuMP
    using BilevelJuMP
    using DataFrames
    using CSV
    using Statistics
    using Gurobi

    include("structures.jl")
    include("input_data.jl")
    include("variables.jl")
    include("constraints.jl")
    include("objective_functions.jl")
    include("problems.jl")
end 
