function clearing(model::Ml, system::Dict, mode::String)

    thermal  = system["thermal"]
    hydro    = system["hydro"]
    line     = system["line"]
    exchange = system["exchange"]
    bus      = system["bus"]
    zone     = system["zone"]
    load     = system["load"]

    maximum_travel_time = get_maximum_travel_time(hydro)

    # Adding decision variables to the optimization model (audited costs or competitive equilibrium)
    add_dispatch_variables!(model, T, thermal, hydro, line, exchange, bus, zone)
    add_hydro_variables!(model, T, hydro, maximum_travel_time)

    # Adding constraints to the optimization model (audited costs or competitive equilibrium)
    add_ramp_constraints!(model, T, hydro, thermal)

    add_generation_capacity_constraints!(model, T, hydro, thermal, "market")
    add_generation_capacity_constraints!(model, T, hydro, thermal, "grid")

    add_balance_constraints!(model, T, system, "market")
    add_balance_constraints!(model, T, system, "grid")

    add_hydro_constraints!(model, T, hydro, maximum_travel_time, "market")
    add_hydro_constraints!(model, T, hydro, maximum_travel_time, "grid")

    # Controlled flows' limits =======> ADD LATER


    # Create objective function
    future_cost_function          = create_future_cost_function(model, T, hydro, mode)
    grid_and_market_cost_function = create_grid_market_cost_function(model, T, hydro, thermal, bus, zone, mode)

    @objective(model, Min, future_cost_function + grid_and_market_cost_function)

    optimize!(model)
    termination_status(model) # =======> ADD LATER

end

function clearing(model::Ml, system::Dict, mode::String, QBidt_EQ::Matrix{Float64}, QBidh_EQ::Matrix{Float64}
                                                         PBidt_EQ::Matrix{Float64}, PBidh_EQ::Matrix{Float64})
end


function audited_costs()

    # Cria modelo JumMP

    # Chama a função clearing para mode = "audited_costs"
end

function competitive_equilibrium()

    # Cria modelo JumMP

    # Chama a função clearing para mode = "competitive_equilibrium"
end

function nash()

    # Cria modelo JumMP

    # Chama a função clearing modificada

    # Roda carrossel
    
end