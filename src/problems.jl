function clearing!(model::Ml, system::Dict, mode::String; T::Int64 = 24) where {Ml}

    thermal  = system["thermal"]
    hydro    = system["hydro"]
    line     = system["line"]
    exchange = system["exchange"]
    bus      = system["bus"]
    zone     = system["zone"]
    load     = system["load"]

    @info("Getting maximum travel time")
    maximum_travel_time = NashEquilibriumElectricityMarkets.get_maximum_travel_time(hydro)

    # Adding decision variables to the optimization model (audited costs or competitive equilibrium)
    @info("Adding dispatch variables")
    NashEquilibriumElectricityMarkets.add_dispatch_variables!(model, T, thermal, hydro, line, exchange, bus, zone)
    @info("Adding hydro variables")
    NashEquilibriumElectricityMarkets.add_hydro_variables!(model, T, hydro, maximum_travel_time)

    # Adding constraints to the optimization model (audited costs or competitive equilibrium)
    @info("Adding ramp constraints")
    NashEquilibriumElectricityMarkets.add_ramp_constraints!(model, T, hydro, thermal)

    @info("Adding generation capacity constraints for mode = market")
    NashEquilibriumElectricityMarkets.add_generation_capacity_constraints!(model, T, hydro, thermal, "market")
    @info("Adding generation capacity constraints for mode = grid")
    NashEquilibriumElectricityMarkets.add_generation_capacity_constraints!(model, T, hydro, thermal, "grid")

    @info("Adding balance constraints for mode = market")
    add_balance_constraints!(model, T, system, "market")
    @info("Adding balance constraints for mode = grid")
    add_balance_constraints!(model, T, system, "grid")

    @info("Adding hydro constraints for mode = market")
    add_hydro_constraints!(model, T, hydro, maximum_travel_time, "market")
    @info("Adding hydro constraints for mode = grid")
    add_hydro_constraints!(model, T, hydro, maximum_travel_time, "grid")

    # Controlled flows' limits =======> ADD LATER


    # Create objective function
    future_cost_function          = create_future_cost_function(model, T, hydro, mode)
    grid_and_market_cost_function = create_grid_market_cost_function(model, T, hydro, thermal, bus, zone, mode)

    @objective(model, Min, future_cost_function + grid_and_market_cost_function)

end

#function clearing(model::Ml, system::Dict, mode::String, QBidt_EQ::Matrix{Float64}, QBidh_EQ::Matrix{Float64},
#                                                         PBidt_EQ::Matrix{Float64}, PBidh_EQ::Matrix{Float64}) where {Ml}
#end

function audited_costs(system::Dict; T::Int64 = 24)

    # Chama a função clearing para mode = "audited_costs"

    model = Model(Gurobi.Optimizer)
    clearing!(model, system, "audited_costs"; T = 24)

    optimize!(model)
    println("Number of solutions: ", result_count(model))
    
    if result_count(model) ≥ 1
        println("Audited Costs: optimal or time limit")
        return create_output_grid(system, model, T)
    else
        println("Audited Costs: infeasible or unbounded")
        return nothing
    end
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

function calc_revenue(system::Dict, output_grid::OutputGrid, output_market::OutputMarket, T::Int64)

    thermal  = system["thermal"]
    hydro    = system["hydro"]
    line     = system["line"]
    exchange = system["exchange"]
    bus      = system["bus"]
    zone     = system["zone"]
    load     = system["load"]

    owner_list_thermal = getfield.(thermal, :owner)
    owner_list_hydro   = getfield.(hydro, :owner)

    owner_list = unique(vcat(owner_list_thermal, owner_list_hydro))

    revenue_nodal  = Vector{Float64}(undef, length(owner_list))
    revenue_zonal  = Vector{Float64}(undef, length(owner_list))
    revenue_uplift = Vector{Float64}(undef, length(owner_list))

    for (k, owner) in enumerate(owner_list)

        revenue_nodal_k  = 0.0
        revenue_zonal_k  = 0.0
        revenue_uplift_k = 0.0

        if owner in owner_list_thermal

            thermal_k = findall(m -> m == owner, getfield.(thermal, :owner))

            revenue_nodal_k += sum(sum((output_grid.nodal_price[t, thermal[i].bus] - thermal[i].uvc) * output_grid.p[t, thermal[i].number] for t in 1:T) 
                                        for i in thermal_k)

            revenue_zonal_k += sum(sum((output_market.zonal_price[t, thermal[i].zone] - thermal[i].uvc) * output_grid.p[t, thermal[i].number] for t in 1:T) 
                                        for i in thermal_k)

            revenue_uplift_k += sum(sum(max((thermal[i].uvc - output_market.zonal_price[t, thermal[i].zone]) * output_grid.p[t, thermal[i].number], 0) for t in 1:T) 
                                        for i in thermal_k)

        end

        if owner in owner_list_hydro

            hydro_k = findall(m -> m == owner, getfield.(hydro, :owner))

            revenue_nodal_k += sum(sum((output_grid.nodal_price[t, hydro[j].bus]) * output_grid.p[t, hydro[j].number] for t in 1:T) 
                                        for j in hydro_k)

            revenue_zonal_k += sum(sum((output_market.zonal_price[t, hydro[j].zone]) * output_grid.p[t, hydro[j].number] for t in 1:T) 
                                        for j in hydro_k)

        end

        revenue_nodal[k]  = revenue_nodal_k
        revenue_zonal[k]  = revenue_zonal_k
        revenue_uplift[k] = revenue_uplift_k
    end

    revenue_df = DataFrame(owner = owner_list, revenue_nodal = revenue_nodal, revenue_zonal = revenue_zonal, revenue_uplift = revenue_uplift)

    return revenue_df
end

function create_output(system::Dict, model::Ml, T::Int64) where {Ml}

    nodal_price = getdual(KCL_grid)
    v_grid      = JuMP.value.(v_grid).data
    q_grid      = JuMP.value.(q_grid).data[end - T + 1:end, :]
    s_grid      = JuMP.value.(s_grid).data[end - T + 1:end, :]
    p_grid      = JuMP.value.(p_grid).data
    g_grid      = JuMP.value.(g_grid).data
    f_grid      = JuMP.value.(f_grid).data
    δ_grid      = JuMP.value.(δ_grid).data
    θ_grid      = JuMP.value.(θ_grid).data

    zonal_price = getdual(KCL_market)
    v_market    = JuMP.value.(v_market).data
    q_market    = JuMP.value.(q_market).data[end - T + 1:end, :]
    s_market    = JuMP.value.(s_market).data[end - T + 1:end, :]
    p_market    = JuMP.value.(p_market).data
    g_market    = JuMP.value.(g_market).data
    f_market    = JuMP.value.(f_market).data
    δ_market    = JuMP.value.(δ_market).data

    output_grid    = OutputGrid(nodal_price, v_grid, q_grid, s_grid, p_grid, g_grid, f_grid, δ_grid, θ_grid)
    output_market  = OutputMarket(zonal_price, v_market, q_market, s_market, p_market, g_market, f_market, δ_market)
    revenue        = calc_revenue(system, output_grid, output_market, T)

    return Output(output_grid, output_market, revenue)

end


