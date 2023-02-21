"Create price bid constraints for generators (hydro and thermal) of specific owner."
function add_price_bid_constraints!(model::Ml, η::Float64, maximum_thermal_uvc::Float64, T::Int64,
                                        hydro::Vector{HydroGenerator},
                                        thermal::Vector{ThermalGenerator}) where {Ml}
    
    Igc, Jgc = length(thermal), length(hydro)
    @constraint(model, [t = 1:T, i = 1:Igc], 0 ≤ model[:λt][t, thermal[i].number] ≤ (1 + η) * thermal[i].uvc)
    @constraint(model, [t = 1:T, j = 1:Jgc], 0 ≤ model[:λh][t, hydro[j].number] ≤ (1 + η) * maximum_thermal_uvc)
end

"Create quantity bid constraints for generators (hydro and thermal) of specific owner."
function add_quantity_bid_constraints!(model::Ml, T::Int64, 
                                            hydro::Vector{HydroGenerator},
                                            thermal::Vector{ThermalGenerator}) where {Ml}

    Igc, Jgc = length(thermal), length(hydro)            
    @constraint(model, [t = 1:T, i = 1:Igc], thermal[i].g_min[t] ≤ model[:μt][t, thermal[i].number] ≤ thermal[i].g_max[t])
    @constraint(model, [t = 1:T, j = 1:Jgc], hydro[j].g_min ≤ model[:μh][t, hydro[j].number] ≤ hydro[j].g_max)
end

"Create ramp constraints for generators (hydro and thermal) of specific/all owners."
function add_ramp_constraints!(model::Ml, T::Int64,
                                hydro::Vector{HydroGenerator}, 
                                thermal::Vector{ThermalGenerator}) where {Ml}

    I, J = length(thermal), length(hydro)          
    @constraint(model, [t = 2:T, i = 1:I], model[:p_grid][t, thermal[i].number] - model[:p_grid][t - 1, thermal[i].number] ≤  thermal[i].ramp_up)       
    @constraint(model, [t = 2:T, i = 1:I], model[:p_grid][t - 1, thermal[i].number] - model[:p_grid][t, thermal[i].number] ≤  thermal[i].ramp_down)         
    @constraint(model, [t = 2:T, j = 1:J], model[:g_grid][t, hydro[j].number] - model[:g_grid][t - 1, hydro[j].number] ≤  hydro[j].ramp_up)       
    @constraint(model, [t = 2:T, j = 1:J], model[:g_grid][t - 1, hydro[j].number] - model[:g_grid][t, hydro[j].number] ≤  hydro[j].ramp_down)                 
end

"Create grid generation constraints for generators (hydro and thermal) of specific owner based on quantity bids."
function add_generation_bid_constraints!(model::Ml, T::Int64,
                                                hydro::Vector{HydroGenerator}, 
                                                thermal::Vector{ThermalGenerator};
                                                grid::Bool = true) where {Ml}

    I, J = length(thermal), length(hydro)
    
    if grid
        @constraint(model, [t = 1:T, i = 1:I], thermal[i].g_min ≤ model[:p_grid][t, thermal[i].number] ≤ model[:μt][t, thermal[i].number])
        @constraint(model, [t = 1:T, j = 1:J], hydro[j].g_min ≤ model[:g_grid][t, hydro[j].number] ≤ model[:μh][t, hydro[j].number])
    else
        @constraint(model, [t = 1:T, i = 1:I], thermal[i].g_min ≤ model[:p_market][t, thermal[i].number] ≤ model[:μt][t, thermal[i].number])
        @constraint(model, [t = 1:T, j = 1:J], hydro[j].g_min ≤ model[:g_market][t, hydro[j].number] ≤ model[:μh][t, hydro[j].number])
    end
end

"Create grid generation constraints for generators (hydro and thermal) of all but specific owner based on quantity bids."
function add_generation_bid_constraints!(model::Ml, T::Int64,
                                                hydro::Vector{HydroGenerator}, 
                                                thermal::Vector{ThermalGenerator},
                                                QBidt_EQ::Matrix{Float64}, QBidh_EQ::Matrix{Float64};
                                                grid::Bool = true) where {Ml}
                                            
    I, J = length(thermal), length(hydro)     

    if grid
        @constraint(model, [t = 1:T, i = 1:I], thermal[i].g_min ≤ model[:p_grid][t, thermal[i].number] ≤ QBidt_EQ[t, thermal[i].number])
        @constraint(model, [t = 1:T, j = 1:J], hydro[j].g_min ≤ model[:g_grid][t, hydro[j].number] ≤ QBidh_EQ[t, hydro[j].number])   
    else
        @constraint(model, [t = 1:T, i = 1:I], thermal[i].g_min ≤ model[:p_market][t, thermal[i].number] ≤ QBidt_EQ[t, thermal[i].number])
        @constraint(model, [t = 1:T, j = 1:J], hydro[j].g_min ≤ model[:g_market][t, hydro[j].number] ≤ QBidh_EQ[t, hydro[j].number])   
    end
end

"Create grid generation constraints for generators (hydro and thermal) of all but specific owner based on generation capacity."
function add_generation_capacity_constraints!(model::Ml, T::Int64,
                                                        hydro::Vector{HydroGenerator}, 
                                                        thermal::Vector{ThermalGenerator};
                                                        grid::Bool = true) where {Ml}

    I, J = length(thermal), length(hydro)     

    if grid
        @constraint(model, [t = 1:T, i = 1:I], thermal[i].g_min ≤ model[:p_grid][t, thermal[i].number] ≤ thermal[i].g_max)
        @constraint(model, [t = 1:T, j = 1:J], hydro[j].g_min ≤ model[:g_grid][t, hydro[j].number] ≤ hydro[j].g_max)
    else
        @constraint(model, [t = 1:T, i = 1:I], thermal[i].g_min ≤ model[:p_market][t, thermal[i].number] ≤ thermal[i].g_max)
        @constraint(model, [t = 1:T, j = 1:J], hydro[j].g_min ≤ model[:g_market][t, hydro[j].number] ≤ hydro[j].g_max)
    end
end

function add_balance_constraints!(model::Ml, T::Int64, system::Dict; grid::Bool = true) where {Ml}
                                                       
    load    = system["load"]     
    thermal = system["thermal"]  
    hydro   = system["hydro"]    

    if grid
        bus  = system["bus"]     
        line = system["line"]    

        B = length(bus)
        L = length(line)

        @constraint(model, KCL_grid[t = 1:T, b = 1:B], sum(model[:p_grid][t, i] for i ∈ findall(i -> i == b, getfield.(thermal, :bus)); init = 0.0) 
                                                     + sum(model[:g_grid][t, j] for j ∈ findall(j -> j == b, getfield.(hydro, :bus)); init = 0.0)
                                                     + sum(model[:f_grid][t, l] for l ∈ findall(l -> l == b, getfield.(line, :to)); init = 0.0)
                                                     - sum(model[:f_grid][t, l] for l ∈ findall(l -> l == b, getfield.(line, :from)); init = 0.0)
                                                     == sum(load[d].value[t] for d ∈ findall(d -> d == b, getfield.(load, :bus)); init = 0.0)
                                                     - model[:δ_grid][t, b])

        @constraint(model, [t = 1:T, l = 1:L], line[l].f_min ≤ model[:f_grid][t, l] ≤ line[l].f_max)

        @constraint(model, [t = 1:T, l = 1:L], model[:f_grid][t, l] == (model[:θ][t, line[l].from] - model[:θ][t, line[l].to])/line[l].reactance)
    else
        zone     = system["zone"]
        exchange = system["exchange"] 

        Z = length(zone)
        E = length(exchange)
        
        @constraint(model, KCL_market[t = 1:T, z = 1:Z], sum(model[:p_market][t, i] for i ∈ findall(i -> i == z, getfield.(thermal, :zone)); init = 0.0) 
                                                        + sum(model[:g_market][t, j] for j ∈ findall(j -> j == z, getfield.(hydro, :zone)); init = 0.0)
                                                        + sum(model[:f_market][t, e] for e ∈ findall(e -> e == z, getfield.(exchange, :to)); init = 0.0)
                                                        - sum(model[:f_market][t, e] for e ∈ findall(e -> e == z, getfield.(exchange, :from)); init = 0.0)
                                                        == sum(load[d].value[t] for d ∈ findall(d -> d == z, getfield.(load, :zone)); init = 0.0)
                                                        - model[:δ_market][t, z])

        @constraint(model, [t = 1:T, e = 1:E], exchange[e].f_min ≤ model[:f_market][t, e] ≤ exchange[e].f_max)

    end
end

function create_future_cost_function(model::Ml, T::Int64, hydro::Vector{HydroGenerator}; nash::Bool = true) where {Ml}

    Jgc = length(hydro)

    if nash

        return sum(hydro[j].γ * (model[:v][T, hydro[j].number] - hydro[j].v_initial) for j = 1:Jgc)
    else

        return sum(hydro[j].γ * (model[:v_grid][T, hydro[j].number] - hydro[j].v_initial) for j = 1:Jgc) + 
                sum(hydro[j].γ * (model[:v_market][T, hydro[j].number] - hydro[j].v_initial) for j = 1:Jgc)
    end
end



using Gurobi, JuMP, BilevelJuMP, CSV, DataFrames

include("src/structures.jl")
include("src/input_data.jl")

path = "data"
T = 24
system = create_input_data(path, T)

model = Model(Gurobi.Optimizer)

#bigN = 1e9
#nashmodel = BilevelModel(Gurobi.Optimizer, mode = BilevelJuMP.FortunyAmatMcCarlMode(primal_big_M = bigN, dual_big_M = bigN))
#set_optimizer_attribute(nashmodel, "NonConvex", 2)

#η = 0.5

I = length(system["thermal"])
J = length(system["hydro"])
L = length(system["line"])
E = length(system["exchange"])
B = length(system["bus"])
Z = length(system["zone"])

# thermal_idx = findall(i -> i == 3, getfield.(system["thermal"], :owner))
# hydro_idx   = findall(j -> j == 3, getfield.(system["hydro"], :owner))

# thermal = system["thermal"][thermal_idx]
# hydro   = system["hydro"][hydro_idx]

# @variable(model, λt[t = 1:T, thermal_idx])
# @variable(model, λh[t = 1:T, hydro_idx])
# @variable(model, μt[t = 1:T, thermal_idx])
# @variable(model, μh[t = 1:T, hydro_idx])
@variable(model, p_grid[t = 1:T, i = 1:I])
@variable(model, g_grid[t = 1:T, j = 1:J])
@variable(model, p_market[t = 1:T, i = 1:I])
@variable(model, g_market[t = 1:T, j = 1:J])
@variable(model, f_grid[t = 1:T, l = 1:L])
@variable(model, f_market[t = 1:T, e = 1:E])
@variable(model, θ[t = 1:T, b = 1:B])
@variable(model, δ_grid[t = 1:T, b = 1:B])
@variable(model, δ_market[t = 1:T, z = 1:Z])
@variable(model, v[1:T, j = 1:J])
@variable(model, v_grid[1:T, j = 1:J])
@variable(model, v_market[1:T, j = 1:J])

# @variable(Upper(nashmodel), λt[t = 1:T, thermal_idx])
# @variable(Upper(nashmodel), λh[t = 1:T, hydro_idx])
# @variable(Upper(nashmodel), μt[t = 1:T, thermal_idx])
# @variable(Upper(nashmodel), μh[t = 1:T, hydro_idx])
# @variable(Upper(nashmodel), p_grid[t = 1:T, thermal_idx])
# @variable(Upper(nashmodel), g_grid[t = 1:T, hydro_idx])

add_price_bid_constraints!(model, η, maximum_thermal_uvc, T, hydro, thermal)
add_price_bid_constraints!(Upper(nashmodel), η, maximum_thermal_uvc, T, hydro, thermal)

add_quantity_bid_constraints!(model, T, hydro, thermal)
add_quantity_bid_constraints!(Upper(nashmodel), T, hydro, thermal)

add_ramp_constraints!(model, T, hydro, thermal)
add_ramp_constraints!(Upper(nashmodel), T, hydro, thermal)

add_balance_constraints!(model, T, system; grid = false) 

