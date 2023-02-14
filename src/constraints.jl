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
function add_grid_generation_bid_constraints!(model::Ml, T::Int64,
                                                hydro::Vector{HydroGenerator}, 
                                                thermal::Vector{ThermalGenerator})

    I, J = length(thermal), length(hydro)     
    @constraint(model, [t = 1:T, i = 1:I], model[:p_grid][t, thermal[i].number] ≥ thermal[i].g_min)
    @constraint(model, [t = 1:T, i = 1:I], model[:p_grid][t, thermal[i].number] ≤ model[:μt][t, thermal[i].number])
    @constraint(model, [t = 1:T, j = 1:J], model[:g_grid][t, hydro[j].number] ≥ hydro[j].g_min)
    @constraint(model, [t = 1:T, j = 1:J], model[:g_grid][t, hydro[j].number] ≤ model[:μh][t, hydro[j].number])
end

"Create grid generation constraints for generators (hydro and thermal) of all but specific owner based on quantity bids."
function add_grid_generation_bid_constraints!(model::Ml, T::Int64,
                                                hydro::Vector{HydroGenerator}, 
                                                thermal::Vector{ThermalGenerator},
                                                QBidt_EQ::Matrix{Float64}, QBidh_EQ::Matrix{Float64} )
                                            
    I, J = length(thermal), length(hydro)     
    @constraint(model, [t = 1:T, i = 1:I], model[:p_grid][t, thermal[i].number] ≥ thermal[i].g_min)
    @constraint(model, [t = 1:T, i = 1:I], model[:p_grid][t, thermal[i].number] ≤ QBidt_EQ[t, thermal[i].number])
    @constraint(model, [t = 1:T, j = 1:J], model[:g_grid][t, hydro[j].number] ≥ hydro[j].g_min)
    @constraint(model, [t = 1:T, j = 1:J], model[:g_grid][t, hydro[j].number] ≤ QBidh_EQ[t, hydro[j].number])                               
end

"Create grid generation constraints for generators (hydro and thermal) of all but specific owner based on generation capacity."
function add_grid_generation_capacity_constraints!(model::Ml, T::Int64,
                                                        hydro::Vector{HydroGenerator}, 
                                                        thermal::Vector{ThermalGenerator})

    I, J = length(thermal), length(hydro)     
    @constraint(model, [t = 1:T, i = 1:I], model[:p_grid][t, thermal[i].number] ≥ thermal[i].g_min)
    @constraint(model, [t = 1:T, i = 1:I], model[:p_grid][t, thermal[i].number] ≤ thermal[i].g_max)
    @constraint(model, [t = 1:T, j = 1:J], model[:g_grid][t, hydro[j].number] ≥ hydro[j].g_min)
    @constraint(model, [t = 1:T, j = 1:J], model[:g_grid][t, hydro[j].number] ≤ hydro[j].g_max)
end


using Gurobi, JuMP, BilevelJuMP

model = Model(Gurobi.Optimizer)

bigN = 1e9
nashmodel = BilevelModel(Gurobi.Optimizer, mode = BilevelJuMP.FortunyAmatMcCarlMode(primal_big_M = bigN, dual_big_M = bigN))
set_optimizer_attribute(nashmodel, "NonConvex", 2)

η = 0.5
T = 24
maximum_thermal_uvc = 5298.

thermal_idx = findall(i -> i == 3, getfield.(system["thermal"], :owner))
hydro_idx   = findall(j -> j == 3, getfield.(system["hydro"], :owner))

thermal = system["thermal"][thermal_idx]
hydro   = system["hydro"][hydro_idx]

@variable(model, λt[t = 1:T, thermal_idx])
@variable(model, λh[t = 1:T, hydro_idx])
@variable(model, μt[t = 1:T, thermal_idx])
@variable(model, μh[t = 1:T, hydro_idx])
@variable(model, p_grid[t = 1:T, thermal_idx])
@variable(model, g_grid[t = 1:T, hydro_idx])

@variable(Upper(nashmodel), λt[t = 1:T, thermal_idx])
@variable(Upper(nashmodel), λh[t = 1:T, hydro_idx])
@variable(Upper(nashmodel), μt[t = 1:T, thermal_idx])
@variable(Upper(nashmodel), μh[t = 1:T, hydro_idx])
@variable(Upper(nashmodel), p_grid[t = 1:T, thermal_idx])
@variable(Upper(nashmodel), g_grid[t = 1:T, hydro_idx])

add_price_bid_constraints!(model, η, maximum_thermal_uvc, T, hydro, thermal)
add_price_bid_constraints!(Upper(nashmodel), η, maximum_thermal_uvc, T, hydro, thermal)

add_quantity_bid_constraints!(model, T, hydro, thermal)
add_quantity_bid_constraints!(Upper(nashmodel), T, hydro, thermal)

add_ramp_constraints!(model, T, hydro, thermal)
add_ramp_constraints!(Upper(nashmodel), T, hydro, thermal)
