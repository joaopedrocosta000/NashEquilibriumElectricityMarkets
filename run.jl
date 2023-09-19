import Pkg
Pkg.activate(".")
Pkg.instantiate()
include("src/NashEquilibriumElectricityMarkets.jl")

# If the packages below are not yet installed, you must install it once using the following command:
# Pkg.add("PACKAGE_NAME")
# After the installation, run the command line below to load the packages.
using JuMP, BilevelJuMP, CSV, DataFrames, Gurobi

#Case 1 - Small test system (3 buses, 3 hydros, 3 thermals, 3 loads)
T            = 24                                 #Number of periods
path         = joinpath(pwd(), "data_ThreeBus")   #Case study folder name, where all input data is stored
price_makers = ["GENCO1","GENCO2","GENCO3"]

#Case 2 - BR system
#T             = 24
#path          = joinpath(pwd(), "data_CaseBR")
#price_makers  = ["GENCO 22","GENCO 43","GENCO 49","GENCO 64"]

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