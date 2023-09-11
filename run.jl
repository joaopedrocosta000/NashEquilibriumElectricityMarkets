using JuMP, BilevelJuMP, CSV, DataFrames, Gurobi
import Pkg
Pkg.activate(".")
Pkg.instantiate()
include("src/NashEquilibriumElectricityMarkets.jl")

#Case 1 - Small test system
#T            = 24
#path         = joinpath(pwd(), "data_ThreeBus")
#price_makers = ["GENCO1","GENCO2","GENCO3"]

#Case 2 - BR system
T             = 24
path          = joinpath(pwd(), "data_CaseBR")
price_makers  = ["GENCO 22","GENCO 43","GENCO 49","GENCO 64"]

#Nash parameters
big_N         = 1e7
price_bid_cap = 0.4
time_limit    = 60 * 60
owner         = 1
price         = "zonal"

#Transforming input data into system object
system = NashEquilibriumElectricityMarkets.create_input_data(path, price_makers; T = 24);

#Run Audited Costs study
output = NashEquilibriumElectricityMarkets.audited_costs(system, path; T = 24);
#Run Competitive Equilibrium study
output = NashEquilibriumElectricityMarkets.competitive_equilibrium(system, path; T = 24);
#Run Nash Equilibrium study
output = NashEquilibriumElectricityMarkets.nash(system, path; T = 24, iteration_max = 100, count_revenue_max = 5, price = "zonal", 
                                                price_bid_cap = 0.4, time_limit = 60 * 60, big_N = 1e6);