"Create future cost function."
function create_future_cost_function(model::Ml, T::Int64, hydro::Vector{HydroGenerator}) where {Ml}
    J = length(hydro)

    return sum(hydro[j].γ * (model[:v_grid][T, hydro[j].number] - hydro[j].v_initial) for j = 1:J) + 
                sum(hydro[j].γ * (model[:v_market][T, hydro[j].number] - hydro[j].v_initial) for j = 1:J)
end

"Create cost function considering a specific owner among all GENCOs." 
function create_grid_market_cost_function!(model::Ml, T::Int64,
                                                hydro::Vector{HydroGenerator}, 
                                                thermal::Vector{ThermalGenerator},
                                                bus::Vector{Bus},
                                                zone::Vector{Zone},
                                                PBidt_EQ::Matrix{Float64}, PBidh_EQ::Matrix{Float64},
                                                genco::Int64) where {Ml}

    thermal_owner = thermal[findall(i -> i == genco, getfield.(thermal, :owner))]
    hydro_owner   = hydro[findall(i -> i == genco, getfield.(hydro, :owner))]

    thermal_system = thermal[findall(i -> i != genco, getfield.(thermal, :owner))]
    hydro_system   = hydro[findall(i -> i != genco, getfield.(hydro, :owner))]

    Igc  = length(thermal_owner)
    Jgc  = length(hydro_owner)
    Isys = length(thermal_system)
    Jsys = length(hydro_system)
    B, Z = length(bus), length(zone)  
    
    if Igc > 0
        thermal_grid_cost_owner    = sum(sum(Upper(model)[:λt][t, thermal_owner[i].number] * Lower(model)[:p_grid][t, thermal_owner[i].number] for i in 1:Igc) for t in 1:T)
        thermal_market_cost_owner  = sum(sum(Upper(model)[:λt][t, thermal_owner[i].number] * Lower(model)[:p_market][t, thermal_owner[i].number] for i in 1:Igc) for t in 1:T)
    else
        thermal_grid_cost_owner   = 0
        thermal_market_cost_owner = 0
    end

    if Jgc > 0
        hydro_grid_cost_owner    = sum(sum(Upper(model)[:λh][t, hydro_owner[j].number] * Lower(model)[:g_grid][t, hydro_owner[j].number] for j in 1:Jgc) for t in 1:T)
        hydro_market_cost_owner  = sum(sum(Upper(model)[:λh][t, hydro_owner[j].number] * Lower(model)[:g_market][t, hydro_owner[j].number] for j in 1:Jgc) for t in 1:T)
    else
        hydro_grid_cost_owner   = 0
        hydro_market_cost_owner = 0
    end

    grid_cost_owner    = sum(thermal_grid_cost_owner + hydro_grid_cost_owner for t in 1:T) 
    market_cost_owner  = sum(thermal_market_cost_owner + hydro_market_cost_owner for t in 1:T) 

    if Isys > 0
        thermal_grid_cost_system    = sum(sum(PBidt_EQ[t, thermal_system[i].number] * Lower(model)[:p_grid][t, thermal_system[i].number] for i in 1:Isys) for t in 1:T)
        thermal_market_cost_system  = sum(sum(PBidt_EQ[t, thermal_system[i].number] * Lower(model)[:p_market][t, thermal_system[i].number] for i in 1:Isys) for t in 1:T)
    else
        thermal_grid_cost_system   = 0
        thermal_market_cost_system = 0
    end

    if Jsys > 0
        hydro_grid_cost_system    = sum(sum(PBidh_EQ[t, hydro_system[j].number] * Lower(model)[:g_grid][t, hydro_system[j].number] for j in 1:Jsys) for t in 1:T)
        hydro_market_cost_system  = sum(sum(PBidh_EQ[t, hydro_system[j].number] * Lower(model)[:g_market][t, hydro_system[j].number] for j in 1:Jsys) for t in 1:T)
    else
        hydro_grid_cost_system   = 0
        hydro_market_cost_system = 0
    end

    grid_cost_system    = sum(thermal_grid_cost_system + hydro_grid_cost_system for t in 1:T) 
    market_cost_system  = sum(thermal_market_cost_system + hydro_market_cost_system for t in 1:T) 
    
    grid_deficit_cost = sum(sum(bus[b].deficit_cost * Lower(model)[:δ_grid][t, b] for b in 1:B) for t in 1:T)

    market_deficit_cost = sum(sum(zone[z].deficit_cost * Lower(model)[:δ_market][t, z] for z in 1:Z) for t in 1:T)
        
    return grid_cost_owner + market_cost_owner + grid_cost_system + market_cost_system + grid_deficit_cost + market_deficit_cost
end

"Create cost function for all GENCOs based on bids."
function create_grid_market_cost_function!(model::Ml, T::Int64,
                                                hydro::Vector{HydroGenerator}, 
                                                thermal::Vector{ThermalGenerator},
                                                bus::Vector{Bus},
                                                zone::Vector{Zone},
                                                PBidt_EQ::Matrix{Float64}, PBidh_EQ::Matrix{Float64};
                                                genco::Int64 = 99999) where {Ml}

    thermal_system = thermal[findall(i -> i != genco, getfield.(thermal, :owner))]
    hydro_system   = hydro[findall(i -> i != genco, getfield.(hydro, :owner))]

    Isys = length(thermal_system)
    Jsys = length(hydro_system)
    B, Z = length(bus), length(zone)  
    
    thermal_grid_cost_system    = sum(sum(PBidt_EQ[t, thermal_system[i].number] * model[:p_grid][t, thermal_system[i].number] for i in 1:Isys; init = 0.0) for t in 1:T)
    hydro_grid_cost_system    = sum(sum(PBidh_EQ[t, hydro_system[j].number] * model[:g_grid][t, hydro_system[j].number] for j in 1:Jsys; init = 0.0) for t in 1:T)

    thermal_market_cost_system  = sum(sum(PBidt_EQ[t, thermal_system[i].number] * model[:p_market][t, thermal_system[i].number] for i in 1:Isys; init = 0.0) for t in 1:T)
    hydro_market_cost_system  = sum(sum(PBidh_EQ[t, hydro_system[j].number] * model[:g_market][t, hydro_system[j].number] for j in 1:Jsys; init = 0.0) for t in 1:T)


    grid_cost_system    = sum(thermal_grid_cost_system + hydro_grid_cost_system for t in 1:T) 
    market_cost_system  = sum(thermal_market_cost_system + hydro_market_cost_system for t in 1:T) 
    
    grid_deficit_cost = sum(sum(bus[b].deficit_cost * model[:δ_grid][t, b] for b in 1:B) for t in 1:T)
    market_deficit_cost = sum(sum(zone[z].deficit_cost * model[:δ_market][t, z] for z in 1:Z) for t in 1:T)
        
    return grid_cost_system + market_cost_system + grid_deficit_cost + market_deficit_cost
end

"Create cost function for all GENCOs based on costs." 
function create_grid_market_cost_function!(model::Ml, T::Int64,
                                                hydro::Vector{HydroGenerator}, 
                                                thermal::Vector{ThermalGenerator},
                                                bus::Vector{Bus},
                                                zone::Vector{Zone}, study::String) where {Ml}
    
    I, J = length(thermal), length(hydro)
    B, Z = length(bus), length(zone)    

    if study == "audited_costs"
        grid_cost = sum(sum(thermal[i].uvc * model[:p_grid][t, thermal[i].number] for i in 1:I) + 
                        sum(bus[b].deficit_cost * model[:δ_grid][t, b] for b in 1:B) for t in 1:T)

        market_cost = sum(sum(thermal[i].uvc * model[:p_market][t, thermal[i].number] for i in 1:I) + 
                        sum(zone[z].deficit_cost * model[:δ_market][t, z] for z in 1:Z) for t in 1:T)
                
        return grid_cost + market_cost
        
    elseif study == "competitive"
        grid_cost = sum(sum(thermal[i].uvc * model[:p_grid][t, thermal[i].number] for i in 1:I) +
                        sum(hydro[j].water_value * model[:g_grid][t, hydro[j].number] for j in 1:J) +
                        sum(bus[b].deficit_cost * model[:δ_grid][t, b] for b in 1:B) for t in 1:T)
                        
        market_cost = sum(sum(thermal[i].uvc * model[:p_market][t, thermal[i].number] for i in 1:I) + 
                        sum(hydro[j].water_value * model[:g_market][t, hydro[j].number] for j in 1:J) +
                        sum(zone[z].deficit_cost * model[:δ_market][t, z] for z in 1:Z) for t in 1:T)
                
        return grid_cost + market_cost
    else
        @error("Study not registered. Only audited_costs or competitive allowed.")
    end
end

"Create revenue function for specific owner."
function create_revenue_function(model::Ml, T::Int64,
                                    hydro::Vector{HydroGenerator},
                                    thermal::Vector{ThermalGenerator},
                                    bus::Vector{Bus},
                                    zone::Vector{Zone},
                                    price::String) where {Ml}      
              
    Igc, Jgc = length(thermal), length(hydro)

    if price == "zonal"
        return sum(sum((Upper(model)[:πz][t, thermal[i].zone] - thermal[i].uvc) * Lower(model)[:p_grid][t, thermal[i].number] for i in 1:Igc; init = 0.0) + 
                      sum((Upper(model)[:πz][t, hydro[j].zone]) * Lower(model)[:g_grid][t, hydro[j].number] for j in 1:Jgc) for t in 1:T)
    else
        return sum(sum((Upper(model)[:πb][t, thermal[i].bus] - thermal[i].uvc) * Lower(model)[:p_grid][t, thermal[i].number] for i in 1:Igc; init = 0.0) + 
                      sum((Upper(model)[:πb][t, hydro[j].bus]) * Lower(model)[:g_grid][t, hydro[j].number] for j in 1:Jgc) for t in 1:T)
    end
end                                   