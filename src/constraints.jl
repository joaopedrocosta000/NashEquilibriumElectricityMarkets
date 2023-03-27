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
                                                thermal::Vector{ThermalGenerator},
                                                problem::String) where {Ml}

    I, J = length(thermal), length(hydro)
    
    if problem == "grid"
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
                                                QBidt_EQ::Matrix{Float64}, QBidh_EQ::Matrix{Float64},
                                                problem::String) where {Ml}
                                            
    I, J = length(thermal), length(hydro)     

    if problem == "grid"
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
                                                        thermal::Vector{ThermalGenerator},
                                                        problem::String) where {Ml}

    I, J = length(thermal), length(hydro)     

    if problem == "grid"
        @constraint(model, [t = 1:T, i = 1:I], thermal[i].g_min ≤ model[:p_grid][t, thermal[i].number] ≤ thermal[i].g_max)
        @constraint(model, [t = 1:T, j = 1:J], hydro[j].g_min ≤ model[:g_grid][t, hydro[j].number] ≤ hydro[j].g_max)
    else
        @constraint(model, [t = 1:T, i = 1:I], thermal[i].g_min ≤ model[:p_market][t, thermal[i].number] ≤ thermal[i].g_max)
        @constraint(model, [t = 1:T, j = 1:J], hydro[j].g_min ≤ model[:g_market][t, hydro[j].number] ≤ hydro[j].g_max)
    end
end

"Create balance constraints for grid and market problems."
function add_balance_constraints!(model::Ml, T::Int64, system::Dict, problem::String) where {Ml}
                                                       
    load    = system["load"]     
    thermal = system["thermal"]  
    hydro   = system["hydro"]    

    if problem == "grid"
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

        @constraint(model, [t = 1:T, b = 1:B], 0 ≤ model[:δ_grid][t, b] ≤ load[findall(d -> d == b, getfield.(load, :bus))].value[t])
        
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

        @constraint(model, [t = 1:T, z = 1:Z], 0 ≤ model[:δ_market][t, z] ≤ sum(load[d].value[t] for d ∈ findall(d -> d == z, getfield.(load, :zone)); init = 0.0))

        @constraint(model, [t = 1:T, e = 1:E], exchange[e].f_min ≤ model[:f_market][t, e] ≤ exchange[e].f_max)

    end
end

"Create hydro constraints for GENCO, grid and market problems."
function add_hydro_constraints!(model::Ml, T::Int64,
                                    hydro::Vector{HydroGenerator}, maximum_travel_time::Int64, problem::String) where {Ml}

    J = length(hydro)     
    c = 1/(10^6/3600)

    if problem == "grid"

        @constraint(model, [t = 1:T, j = 1:J], hydro[j].q_min ≤ model[:q_grid][t, hydro[j].number] ≤ hydro[j].q_max)
        @constraint(model, [t = 1:T, j = 1:J], 0 ≤ model[:s_grid][t, hydro[j].number])
        @constraint(model, [t = 1:T, j = 1:J], hydro[j].v_min ≤ model[:v_grid][t, hydro[j].number] ≤ hydro[j].v_max)

        @constraint(model, [t = 1:T, j = 1:J], model[:g_grid][t, hydro[j].number] == hydro[j].ρ * model[:q_grid][t, hydro[j].number])
        @constraint(model, [j = 1:J], model[:v_grid][0, hydro[j].number] == hydro[j].v_initial)

        @constraint(model, [t = (-maximum_travel_time + 1):0, j = 1:J], model[:q_grid][t, hydro[j].number] == hydro[j].q_initial[t + maximum_travel_time])
        @constraint(model, [t = (-maximum_travel_time + 1):0, j = 1:J], model[:s_grid][t, hydro[j].number] == hydro[j].s_initial[t + maximum_travel_time])

        @constraint(model, [t = 1:T, j = 1:J], model[:v_grid][t, hydro[j].number] == model[:v_grid][t - 1, hydro[j].number] - 
                                                c * (model[:q_grid][t - 1, hydro[j].number] + model[:s_grid][t - 1, hydro[j].number] - hydro[j].inflow[t] -
                                                sum(model[:q_grid][t - hydro[j].cascade[m], m] + model[:s_grid][t - hydro[j].cascade[m], m] for m in keys(hydro[j].cascade))))

    elseif problem == "market"

        @constraint(model, [t = 1:T, j = 1:J], hydro[j].q_min ≤ model[:q_market][t, hydro[j].number] ≤ hydro[j].q_max)
        @constraint(model, [t = 1:T, j = 1:J], 0 ≤ model[:s_market][t, hydro[j].number])
        @constraint(model, [t = 1:T, j = 1:J], hydro[j].v_min ≤ model[:v_market][t, hydro[j].number] ≤ hydro[j].v_max)

        @constraint(model, [t = 1:T, j = 1:J], model[:g_market][t, hydro[j].number] == hydro[j].ρ * model[:q_market][t, hydro[j].number])
        @constraint(model, [j = 1:J], model[:v_market][0, hydro[j].number] == hydro[j].v_initial)

        @constraint(model, [t = (-maximum_travel_time + 1):0, j = 1:J], model[:q_market][t, hydro[j].number] == hydro[j].q_initial[t + maximum_travel_time])
        @constraint(model, [t = (-maximum_travel_time + 1):0, j = 1:J], model[:s_market][t, hydro[j].number] == hydro[j].s_initial[t + maximum_travel_time])

        @constraint(model, [t = 1:T, j = 1:J], model[:v_market][t, hydro[j].number] == model[:v_market][t - 1, hydro[j].number] - 
                                                c * (model[:q_market][t - 1, hydro[j].number] + model[:s_market][t - 1, hydro[j].number] - hydro[j].inflow[t] -
                                                sum(model[:q_market][t - hydro[j].cascade[m], m] + model[:s_market][t - hydro[j].cascade[m], m] for m in keys(hydro[j].cascade))))

    else
        
        @constraint(model, [t = 1:T, j = 1:J], hydro[j].q_min ≤ model[:q][t, hydro[j].number] ≤ hydro[j].q_max)
        @constraint(model, [t = 1:T, j = 1:J], 0 ≤ model[:s][t, hydro[j].number])
        @constraint(model, [t = 1:T, j = 1:J], hydro[j].v_min ≤ model[:v][t, hydro[j].number] ≤ hydro[j].v_max)

        @constraint(model, [t = 1:T, j = 1:J], model[:g_grid][t, hydro[j].number] == hydro[j].ρ * model[:q][t, hydro[j].number])
        @constraint(model, [j = 1:J], model[:v][0, hydro[j].number] == hydro[j].v_initial)

        @constraint(model, [t = (-maximum_travel_time + 1):0, j = 1:J], model[:q][t, hydro[j].number] == hydro[j].q_initial[t + maximum_travel_time])
        @constraint(model, [t = (-maximum_travel_time + 1):0, j = 1:J], model[:s][t, hydro[j].number] == hydro[j].s_initial[t + maximum_travel_time])

        @constraint(model, [t = 1:T, j = 1:J], model[:v][t, hydro[j].number] == model[:v][t - 1, hydro[j].number] - 
                                                c * (model[:q][t - 1, hydro[j].number] + model[:s][t - 1, hydro[j].number] - hydro[j].inflow[t] -
                                                sum(model[:q][t - hydro[j].cascade[m], m] + model[:s][t - hydro[j].cascade[m], m] for m in keys(hydro[j].cascade))))
    end
end



