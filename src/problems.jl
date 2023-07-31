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
    NashEquilibriumElectricityMarkets.add_balance_constraints!(model, T, system, "market")
    @info("Adding balance constraints for mode = grid")
    NashEquilibriumElectricityMarkets.add_balance_constraints!(model, T, system, "grid")

    @info("Adding hydro constraints for mode = market")
    add_hydro_constraints!(model, T, hydro, maximum_travel_time, "market")
    @info("Adding hydro constraints for mode = grid")
    add_hydro_constraints!(model, T, hydro, maximum_travel_time, "grid")

    # Controlled flows' limits =======> ADD LATER


    # Create objective function
    @info("Creating objective function")
    future_cost_function          = NashEquilibriumElectricityMarkets.create_future_cost_function(model, T, hydro)
    grid_and_market_cost_function = NashEquilibriumElectricityMarkets.create_grid_market_cost_function(model, T, hydro, thermal, bus, zone, mode)

    @objective(model, Min, grid_and_market_cost_function - future_cost_function)
    # @objective(model, Min, grid_and_market_cost_function)
end

function clearing!(model::Ml, system::Dict, QBidt_EQ::Matrix{Float64}, QBidh_EQ::Matrix{Float64},
                                                        PBidt_EQ::Matrix{Float64}, PBidh_EQ::Matrix{Float64}, genco::Int64; T::Int64 = 24) where {Ml}

    thermal  = system["thermal"]
    hydro    = system["hydro"]
    line     = system["line"]
    exchange = system["exchange"]
    bus      = system["bus"]
    zone     = system["zone"]
    load     = system["load"]

    @info("Getting maximum travel time")
    maximum_travel_time = NashEquilibriumElectricityMarkets.get_maximum_travel_time(hydro)

    # Adding decision variables to the bilevel optimization model
    @info("Adding dispatch variables")
    NashEquilibriumElectricityMarkets.add_dispatch_variables!(Lower(model), T, thermal, hydro, line, exchange, bus, zone)
    @info("Adding hydro variables")
    NashEquilibriumElectricityMarkets.add_hydro_variables!(Lower(model), T, hydro, maximum_travel_time)

    # Adding constraints to the bilevel optimization model
    @info("Adding ramp constraints")
    NashEquilibriumElectricityMarkets.add_ramp_constraints!(Lower(model), T, hydro, thermal)

    @info("Adding generation bid constraints for mode = market")
    NashEquilibriumElectricityMarkets.add_generation_bid_constraints!(model, T, hydro, thermal, QBidt_EQ, QBidh_EQ, "market", genco)

    @info("Adding generation bid constraints for mode = grid")
    NashEquilibriumElectricityMarkets.add_generation_bid_constraints!(model, T, hydro, thermal, QBidt_EQ, QBidh_EQ, "grid", genco)

    @info("Adding balance constraints for mode = market")
    NashEquilibriumElectricityMarkets.add_balance_constraints!(Lower(model), T, system, "market")
    @info("Adding balance constraints for mode = grid")
    NashEquilibriumElectricityMarkets.add_balance_constraints!(Lower(model), T, system, "grid")

    @info("Adding hydro constraints for mode = market")
    add_hydro_constraints!(Lower(model), T, hydro, maximum_travel_time, "market")
    @info("Adding hydro constraints for mode = grid")
    add_hydro_constraints!(Lower(model), T, hydro, maximum_travel_time, "grid")

    # Controlled flows' limits =======> ADD LATER

    # Create objective function
    @info("Creating objective function")
    future_cost_function          = NashEquilibriumElectricityMarkets.create_future_cost_function(Lower(model), T, hydro)
    grid_and_market_cost_function = NashEquilibriumElectricityMarkets.create_grid_market_cost_function(model, T, hydro, thermal, bus, zone, PBidt_EQ, PBidh_EQ, genco)

    @objective(Lower(model), Min, grid_and_market_cost_function - future_cost_function)

end

function audited_costs(system::Dict; T::Int64 = 24)
    model = Model(Gurobi.Optimizer)
    # set_optimizer_attribute(model, "DualReductions", 0)
    clearing!(model, system, "audited_costs"; T = 24)

    optimize!(model)
    println("Termination status: ", termination_status(model))
    println("Number of solutions: ", result_count(model))
    
    if result_count(model) ≥ 1
        println("Audited Costs: optimal or time limit")
        return create_output(system, model, T)
    else
        println("Audited Costs: infeasible or unbounded")
        return nothing
    end
end

function competitive_equilibrium(system::Dict; T::Int64 = 24)
    model = Model(Gurobi.Optimizer)
    # set_optimizer_attribute(model, "DualReductions", 0)
    clearing!(model, system, "competitive"; T = 24)

    optimize!(model)
    println("Termination status: ", termination_status(model))
    println("Number of solutions: ", result_count(model))
    
    if result_count(model) ≥ 1
        println("Competitive Equilibrium: optimal or time limit")
        return create_output(system, model, T)
    else
        println("Competitive Equilibrium: infeasible or unbounded")
        return nothing
    end
end

function optimal_bid(system::Dict, owner::Int64, price::String,
                        PBidt_EQ::Matrix{Float64}, PBidh_EQ::Matrix{Float64},
                        QBidt_EQ::Matrix{Float64}, QBidh_EQ::Matrix{Float64}; 
                        T::Int64 = 24, price_bid_cap::Float64 = 0.4,
                        time_limit::Union{Bool, Int64} = 60 * 60, big_N = 1e7)

    model = BilevelModel(Gurobi.Optimizer, mode = BilevelJuMP.FortunyAmatMcCarlMode(primal_big_M = big_N, dual_big_M = big_N))
    set_optimizer_attribute(model, "NonConvex", 2)

    if typeof(time_limit) == Int64
        set_optimizer_attribute(model, "TimeLimit", time_limit)
    end

    thermal_owner_idx  = findall(i -> i == owner, getfield.(system["thermal"], :owner))
    hydro_owner_idx    = findall(i -> i == owner, getfield.(system["hydro"], :owner))
    thermal_owner      = system["thermal"][thermal_owner_idx]
    hydro_owner        = system["hydro"][hydro_owner_idx]

    thermal        = system["thermal"]
    hydro          = system["hydro"]
    line           = system["line"]
    exchange       = system["exchange"]
    bus            = system["bus"]
    zone           = system["zone"]
    load           = system["load"]

    #PBidt_EQ, PBidh_EQ, QBidt_EQ, QBidh_EQ = NashEquilibriumElectricityMarkets.initialize_bids(system, T)

    @info("---------- Creating upper level variables ----------")
    # Adding decision variables for the upper level model
    @info("Adding price bid variables for genco $(owner)")
    NashEquilibriumElectricityMarkets.add_price_bid_variables!(Upper(model), T, thermal_owner, hydro_owner)
    @info("Adding quantity bid variables for genco $(owner)")
    NashEquilibriumElectricityMarkets.add_quantity_bid_variables!(Upper(model), T, thermal_owner, hydro_owner)

    @info("---------- Creating lower level complete model ----------")
    NashEquilibriumElectricityMarkets.clearing!(model, system, QBidt_EQ, QBidh_EQ, PBidt_EQ, PBidh_EQ, owner; T = T)

    @info("---------- Creating dual variables of balance constraints ----------")
    NashEquilibriumElectricityMarkets.add_shadow_price_variable!(model, T, bus, zone)
                                    
    @info("---------- Creating upper level constraints ----------")
    @info("Adding price bid constraints for genco $(owner)")
    NashEquilibriumElectricityMarkets.add_price_bid_constraints!(Upper(model), price_bid_cap, T, hydro_owner, thermal_owner)
    @info("Adding quantity bid constraints for genco $(owner)")
    NashEquilibriumElectricityMarkets.add_quantity_bid_constraints!(Upper(model), T, hydro_owner, thermal_owner)
    
    @info("---------- Creating upper level obj. function ----------")
    @info("Creating objective function for genco $(owner)")
    revenue_owner = NashEquilibriumElectricityMarkets.create_revenue_function(model, T, hydro_owner, thermal_owner, bus, zone, price)   
                                        
    @objective(Upper(model), Max, revenue_owner)

    optimize!(model)
    println("Termination status: ", termination_status(model))
    println("Number of solutions: ", result_count(model))

    return round.(JuMP.value.(λt).data, digits = 2), round.(JuMP.value.(λh).data, digits = 2), 
                round.(JuMP.value.(μt).data, digits = 2), round.(JuMP.value.(μh).data, digits = 2), result_count(model)
end

function nash(system::Dict; T::Int64 = 24, iteration_max::Int64 = 100, count_revenue_max::Int64 = 5, 
                price::String = "zonal", price_bid_cap::Float64 = 0.4, time_limit = 60 * 60, big_N = 1e6)

    PBidt_EQ, PBidh_EQ, QBidt_EQ, QBidh_EQ = NashEquilibriumElectricityMarkets.initialize_bids(system, T)

    PBidt_carousel = deepcopy(PBidt_EQ)
    PBidh_carousel = deepcopy(PBidh_EQ)
    QBidt_carousel = deepcopy(QBidt_EQ)
    QBidh_carousel = deepcopy(QBidh_EQ)

    price_maker        = getfield.(system["company"], :price_maker)
    n_price_makers     = count(price_maker)
    owner_price_makers = getfield.(system["company"][price_maker], :number)

    flag_keep_bid        = falses(n_price_makers)
    flag_opt             = falses(n_price_makers)
    revenue_price_makers = zeros(n_price_makers)
    iteration            = 0
    count_revenue        = 0
    
    while (flag_keep_bid ≠ trues(n_price_makers) || flag_opt ≠ trues(n_price_makers)) && count_revenue < count_revenue_max && iteration < iteration_max

        flag_keep_bid = falses(n_price_makers)
        flag_opt      = falses(n_price_makers)
        iteration     += 1

        @info("---------- Iteration $iteration ----------")
        for (i, owner) in enumerate(owner_price_makers)

            flag_keep_bid_owner_thermal = false
            flag_keep_bid_owner_hydro   = false

            @info("Price Maker $owner")
            λt, λh, μt, μh, solution_count = optimal_bid(system, owner, price, PBidt_EQ, PBidh_EQ, QBidt_EQ, QBidh_EQ; 
                                            T = T, price_bid_cap = price_bid_cap,
                                            time_limit = time_limit, big_N = big_N)

            thermal_owner_idx  = findall(i -> i == owner, getfield.(system["thermal"], :owner))
            hydro_owner_idx    = findall(i -> i == owner, getfield.(system["hydro"], :owner))

            if solution_count > 0
                @info("Optimal solutions found for owner $(owner)!")
                flag_opt[i] = true
        
                if !isempty(thermal_owner_idx)
                    PBidt_carousel[:, thermal_owner_idx] = μt
                    QBidt_carousel[:, thermal_owner_idx] = λt

                    if (PBidt_carousel[:, thermal_owner_idx] == PBidt_EQ[:, thermal_owner_idx]) && (QBidt_carousel[:, thermal_owner_idx] == QBidt_EQ[:, thermal_owner_idx])
                        flag_keep_bid_owner_thermal = true
                    end
                end

                if !isempty(hydro_owner_idx)
                    PBidh_carousel[:, hydro_owner_idx] = μh
                    QBidh_carousel[:, hydro_owner_idx] = λh

                    if (PBidh_carousel[:, thermal_owner_idx] == PBidh_EQ[:, thermal_owner_idx]) && (QBidh_carousel[:, thermal_owner_idx] == QBidh_EQ[:, thermal_owner_idx])
                        flag_keep_bid_owner_hydro = true
                    end
                end

                if flag_keep_bid_owner_thermal && flag_keep_bid_owner_hydro
                    flag_keep_bid[i] = true
                    @info("Player $(owner) kept bid!")
                else
                    @info("Player $(owner) updated bid!")    
                end
            else
                @info("Optimal bid not found for owner $(owner)!")
            end
        end

        PBidt_EQ[:, :] = PBidt_carousel[:, :]
        PBidh_EQ[:, :] = PBidh_carousel[:, :]
        QBidt_EQ[:, :] = QBidt_carousel[:, :]
        QBidh_EQ[:, :] = QBidh_carousel[:, :]

        model = Model(Gurobi.Optimizer)

        clearing!(model, system, QBidt_EQ, QBidh_EQ, PBidt_EQ, PBidh_EQ, owner; T = T)

        optimize!(model)
        println("Termination status: ", termination_status(model))
        println("Number of solutions: ", result_count(model))
        
        if result_count(model) ≥ 1
            @info("Clearing: optimal or time limit")
            output_clearing = create_output(system, model, T)
        else
            @info("Clearing: infeasible or unbounded")
            output_clearing = nothing
        end
        
        # ARMAZENAR RESULTADOS DO CLEARING
        # CALCULAR count_revenue
        # EXPORTAR RESULTADOS (FAZER TAMBÉM PARA AS FUNÇÕES DE CUSTOS AUDITADOS E EQUILÍBRIO COMPETITIVO)
       
        println("Equilibrium Vector: ", flag_keep_bid)
        println("Optimal Vector: ", flag_opt)
        println("Revenue flag: ", flag_rev)

        if all(flag_opt)

            if all(flag_keep_bid)
                println("Equilibrium point found.")
            else
                if count_revenue ≥ count_revenue_max 
                    println("Equilibrium point found through revenue.")
                else 
                    println("Equilibrium point not found.")
                end
            end
        else

            if all(flag_keep_bid)
                println("Equilibrium candidate point found. New iteration required.")
            else
                if count_revenue ≥ count_revenue_max 
                    println("Equilibrium candidate point found through revenue. New iteration required.")
                    count_revenue_max = count_revenue_max + 1
                else 
                    println("Equilibrium point not found.")
                end
            end
        end

    end

    # Roda carrossel - 2 MATRIZES, UMA DE EQUILÍBRIO OUTRA DO CARROSSEL; NO FINAL DE CADA PASSAGEM DO WHILE, ELAS VÃO SER IGUAIS
    # VAMOS PREENCHENDO A CARROSSEL NAS COLUNAS ASSOCIADAS AO GENCO SENDO OTIMIZADO
    
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

    owner_list = sort(unique(vcat(owner_list_thermal, owner_list_hydro)))

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

            revenue_nodal_k += sum(sum((output_grid.nodal_price[t, hydro[j].bus]) * output_grid.g[t, hydro[j].number] for t in 1:T) 
                                        for j in hydro_k)

            revenue_zonal_k += sum(sum((output_market.zonal_price[t, hydro[j].zone]) * output_grid.g[t, hydro[j].number] for t in 1:T) 
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

    nodal_price = dual.(model[:KCL_grid])
    v_grid      = JuMP.value.(model[:v_grid]).data
    q_grid      = JuMP.value.(model[:q_grid]).data[end - T + 1:end, :]
    s_grid      = JuMP.value.(model[:s_grid]).data[end - T + 1:end, :]
    p_grid      = JuMP.value.(model[:p_grid])
    g_grid      = JuMP.value.(model[:g_grid])
    f_grid      = JuMP.value.(model[:f_grid])
    δ_grid      = JuMP.value.(model[:δ_grid])
    θ           = JuMP.value.(model[:θ])

    zonal_price = dual.(model[:KCL_market])
    v_market    = JuMP.value.(model[:v_market]).data
    q_market    = JuMP.value.(model[:q_market]).data[end - T + 1:end, :]
    s_market    = JuMP.value.(model[:s_market]).data[end - T + 1:end, :]
    p_market    = JuMP.value.(model[:p_market])
    g_market    = JuMP.value.(model[:g_market])
    f_market    = JuMP.value.(model[:f_market])
    δ_market    = JuMP.value.(model[:δ_market])

    output_grid    = OutputGrid(nodal_price, v_grid, q_grid, s_grid, p_grid, g_grid, f_grid, δ_grid, θ)
    output_market  = OutputMarket(zonal_price, v_market, q_market, s_market, p_market, g_market, f_market, δ_market)
    revenue        = calc_revenue(system, output_grid, output_market, T)

    return Output(output_grid, output_market, revenue)

end

function initialize_bids(system::Dict, T::Int64)

    thermal = system["thermal"]
    hydro   = system["hydro"]

    I, J = length(thermal), length(hydro)

    PBidt_EQ = Matrix{Float64}(undef, T, I)
    PBidh_EQ = Matrix{Float64}(undef, T, J)
    QBidt_EQ = Matrix{Float64}(undef, T, I)
    QBidh_EQ = Matrix{Float64}(undef, T, J)

    for i in 1:I
        PBidt_EQ[:, i] = ones(T) .* thermal[i].uvc
        QBidt_EQ[:, i] = thermal[i].g_max
    end

    for j in 1:J
        PBidh_EQ[:, j] = ones(T) .* hydro[j].water_value
        QBidh_EQ[:, j] = ones(T) .* hydro[j].g_max
    end

    return PBidt_EQ, PBidh_EQ, QBidt_EQ, QBidh_EQ
end

