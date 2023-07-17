"Create price bid variables for a specific owner"
function add_price_bid_variables!(model::Ml, T::Int64, 
                                thermal::Vector{ThermalGenerator},
                                hydro::Vector{HydroGenerator}) where {Ml}
  
    thermal_idx = getfield.(thermal, :number)   
    hydro_idx   = getfield.(hydro, :number)            
                        
    @variable(model, λt[t = 1:T, thermal_idx])
    @variable(model, λh[t = 1:T, hydro_idx])
end
"Create quantity bid variables for a specific owner"
function add_quantity_bid_variables!(model::Ml, T::Int64, 
                                    thermal::Vector{ThermalGenerator},
                                    hydro::Vector{HydroGenerator}) where {Ml}

    thermal_idx = getfield.(thermal, :number)   
    hydro_idx   = getfield.(hydro, :number)    
    
    @variable(model, μt[t = 1:T, thermal_idx])
    @variable(model, μh[t = 1:T, hydro_idx])
end

function add_dispatch_variables!(model::Ml, T::Int64,
                                 thermal::Vector{ThermalGenerator},
                                 hydro::Vector{HydroGenerator},
                                 line::Vector{Line},
                                 exchange::Vector{Exchange},
                                 bus::Vector{Bus},
                                 zone::Vector{Zone}) where {Ml}
                            
    I, J = length(thermal), length(hydro)
    B, Z = length(bus), length(zone)
    L, E = length(line), length(exchange)

    @variable(model, p_grid[t = 1:T, i = 1:I])
    @variable(model, g_grid[t = 1:T, j = 1:J])

    @variable(model, p_market[t = 1:T, i = 1:I])
    @variable(model, g_market[t = 1:T, j = 1:J])

    @variable(model, f_grid[t = 1:T, l = 1:L])
    @variable(model, f_market[t = 1:T, e = 1:E])

    @variable(model, θ[t = 1:T, b = 1:B])

    @variable(model, δ_grid[t = 1:T, b = 1:B])
    @variable(model, δ_market[t = 1:T, z = 1:Z])
end

function add_hydro_variables!(model::Ml, T::Int64,
                              hydro::Vector{HydroGenerator},
                              maximum_travel_time::Int64) where {Ml}

    J = length(hydro)

    @variable(model, v_grid[0:T, j = 1:J])
    @variable(model, v_market[0:T, j = 1:J])
    @variable(model, q_grid[-maximum_travel_time + 1:T, j = 1:J])
    @variable(model, q_market[-maximum_travel_time + 1:T, j = 1:J])
    @variable(model, s_grid[-maximum_travel_time + 1:T, j = 1:J])
    @variable(model, s_market[-maximum_travel_time + 1:T, j = 1:J])
end

"Market and Grid Operators' clearing price as the dual variable balance constraints"
function add_shadow_price_variable!(model::Ml, T::Int64, 
                                        bus::Vector{Bus},
                                        zone::Vector{Zone}) where {Ml}

    B = length(bus)
    Z = length(zone)
                                        
    @variable(Upper(model), πb[t = 1:T, b = 1:B] >= 0, BilevelJuMP.DualOf.(Lower(model)[:KCL_grid][t, b]))
    @variable(Upper(model), πz[t = 1:T, z = 1:Z] >= 0, BilevelJuMP.DualOf.(Lower(model)[:KCL_market][t, z]))
end