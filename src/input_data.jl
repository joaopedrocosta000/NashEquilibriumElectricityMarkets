"Creates a vector of bus structs."
function create_bus_structs(path::String)

    bus_df = CSV.read(joinpath(path, "Bus.csv"), DataFrame)
    B      = size(bus_df, 1)

    bus    = Vector{Bus}(undef, B)

    for b in 1:B
        bus[b] = Bus(bus_df[b, "Bus"], bus_df[b, "Zone"], bus_df[b, "Bus_Name"], bus_df[b, "Cdef"])
    end

    return bus
end

"Creates a vector of zone structs."
function create_zone_structs(path::String)

    zone_df = CSV.read(joinpath(path, "Bus.csv"), DataFrame)

    unique_number = unique(zone_df[:, "Zone"])
    unique_name   = unique(zone_df[:, "Zone_Name"])
    
    Z    = length(unique_number)
    zone = Vector{Zone}(undef, Z)

    for z in 1:Z
        zone_bus_df = filter(:Zone => isequal(z), zone_df)
        zone[z] = Zone(unique_number[z], unique_name[z], round(mean(zone_bus_df[:, :Cdef]); digits = 2))
    end

    return zone
end

"Creates a vector of lines structs."
function create_line_structs(path::String)

    line_df = CSV.read(joinpath(path, "Lines.csv"), DataFrame)
    L       = size(line_df, 1)

    line    = Vector{Line}(undef, L)

    for l in 1:L
        line[l] = Line(line_df[l, "Line"], line_df[l, "from"], line_df[l, "to"], line_df[l, "Fmin"], line_df[l, "Fmax"], line_df[l, "x"])
    end

    return line
end

"Creates a vector of exchanges structs."
function create_exchange_structs(path::String)

    exchange_df = CSV.read(joinpath(path, "LinesZ.csv"), DataFrame)
    E           = size(exchange_df, 1)

    exchange    = Vector{Exchange}(undef, E)

    for e in 1:E
        exchange[e] = Exchange(exchange_df[e, "Line"], exchange_df[e, "from"], exchange_df[e, "to"], exchange_df[e, "Fmin"], exchange_df[e, "Fmax"])
    end

    return exchange
end

"Creates a vector of load structs."
function create_load_structs(path::String)

    load_df = CSV.read(joinpath(path, "Load.csv"), DataFrame)
    D       = size(load_df, 1)

    load    = Vector{Load}(undef, D)

    for d in 1:D
        load[d] = Load(load_df[d, "Bus"], load_df[d, "Zone"], Vector(load_df[d, 3:end]))
    end

    return load
end

"Creates a vector of thermal structs."
function create_thermal_structs(path::String, T::Int64)

    thermal_df = CSV.read(joinpath(path, "Thermal.csv"), DataFrame)
    I          = size(thermal_df, 1)

    thermal    = Vector{ThermalGenerator}(undef, I)

    for i in 1:I
        thermal[i] = ThermalGenerator(thermal_df[i, "TGenerator"], thermal_df[i, "TGenerator_cod"], thermal_df[i, "Name"],
                                        thermal_df[i, "GENCO"], thermal_df[i, "Gtbus"], thermal_df[i, "Gtzone"], 
                                        Vector(thermal_df[i, 10:10 + T - 1]), Vector(thermal_df[i, 10 + T:10 + 2T - 1]),
                                        thermal_df[i, "Rup"], thermal_df[i, "Rdown"], thermal_df[i, "UVC"], thermal_df[i, "Fuel"])
    end

    return thermal
end

"Creates a vector of hydro structs."
function create_hydro_structs(path::String)

    hydro_df            = CSV.read(joinpath(path, "Hydro.csv"), DataFrame)
    cascade_df          = CSV.read(joinpath(path, "Hydro_Cascade.csv"), DataFrame)[:, 2:end]
    inflow_df           = CSV.read(joinpath(path, "Hydro_Inflow.csv"), DataFrame)[:, 2:end]
    spillage_df         = CSV.read(joinpath(path, "Hydro_SpillageHist.csv"), DataFrame)
    turbined_outflow_df = CSV.read(joinpath(path, "Hydro_TurbinedOutflowHist.csv"), DataFrame)

    J = size(hydro_df, 1)
    maximum_τ = maximum(Matrix(cascade_df)[:, 2:end])

    hydro = Vector{HydroGenerator}(undef, J)
   
    for j in 1:J

        ρ           = hydro_df[j, "p"]
        ρc          = hydro_df[j, "pc"]
        q_max       = hydro_df[j, "Qhmax"]
        q_min       = hydro_df[j, "Qhmin"]
        water_value = hydro_df[j, "WV"]
        c           = 1/(10^6/3600)

        g_max = round(q_max * ρ, digits = 1)
        g_min = round(q_min * ρ, digits = 1)
        γ     = round(water_value * ρc / c, digits = 2)

        dict_cascade = Dict{Int64, Int64}()
        
        cascade_jm = findall(i -> i != 0, cascade_df[:, j])

        for n in cascade_jm
            dict_cascade[n] = cascade_df[n, j]
        end

        hydro[j] = HydroGenerator(hydro_df[j, "HGenerator"], hydro_df[j, "HGenerator_cod"], hydro_df[j, "Name"],
                                    hydro_df[j, "GENCO"], hydro_df[j, "Ghbus"], hydro_df[j, "Ghzone"],
                                    ρ, ρc, g_max, g_min, round(q_max, digits = 2), q_min, hydro_df[j, "Vmax"],
                                    hydro_df[j, "Vmin"], hydro_df[j, "Vi"], Vector(turbined_outflow_df[j, 2:maximum_τ + 1]),
                                    Vector(spillage_df[j, 2:maximum_τ + 1]), hydro_df[j, "Rup"], hydro_df[j, "Rdown"],
                                    water_value, γ, inflow_df[:, j], dict_cascade)
    end

    return hydro
end

"Creates a vector of company structs."
function create_company_structs(path::String, T::Int64, price_makers::Vector{String})

    thermal_df = CSV.read(joinpath(path, "Thermal.csv"), DataFrame)
    hydro_df   = CSV.read(joinpath(path, "Hydro.csv"), DataFrame)

    unique_number = unique(vcat(thermal_df[:, "GENCO"], hydro_df[:, "GENCO"]))
    unique_name   = unique(vcat(thermal_df[:, "GENCO_Name"], hydro_df[:, "GENCO_Name"]))

    company_df = DataFrames.sort(DataFrame(GENCO = unique_number, GENCO_Name = unique_name), "GENCO")
    # company_df = classify_price_makers(thermal_df, hydro_df, company_df, T)
    C       = length(unique_number)
    company = Vector{Company}(undef, C)

    for c in 1:C
        
        if String(company_df[c, "GENCO_Name"]) in price_makers
            is_price_maker = true
        else
            is_price_maker = false
        end

        company[c] = Company(company_df[c, "GENCO"], String(company_df[c, "GENCO_Name"]), is_price_maker)
    end
    
    return company
end

"Verifies if a generator is a small gen."
is_small_gen(name::Sl) where {Sl} = split(name, " ")[1] == "Small"

"Alternative classification for companies as prices makers or price takers."
function classify_price_makers(thermal_df::DataFrame, hydro_df::DataFrame, company_df::DataFrame, T::Int64) 

    ρ           = hydro_df[:, "p"]
    q_max       = hydro_df[:, "Qhmax"]

    g_max = q_max .* ρ

    thermal_genco_df = DataFrame(GENCO = thermal_df[:, "GENCO"], GENCO_Name = thermal_df[:, "GENCO_Name"],
                                    Capacity = maximum(Matrix(thermal_df[:, 10:10 + T - 1]), dims = 2)[:, 1])
    
    hydro_genco_df = DataFrame(GENCO = hydro_df[:, "GENCO"], GENCO_Name = hydro_df[:, "GENCO_Name"],
                                    Capacity = g_max)

    genco_df              = groupby(vcat(thermal_genco_df, hydro_genco_df), "GENCO")
    sum_genco_capacity_df = combine(genco_df, :Capacity => sum)                                

    sum_genco_capacity_df =  DataFrames.sort(rename(hcat(company_df[:, "GENCO_Name"], sum_genco_capacity_df), Dict("x1" => "GENCO_Name")), 
                                                "Capacity_sum", rev = true)

    price_maker_candidates_df = filter(:GENCO_Name => !is_small_gen, sum_genco_capacity_df)
    total_capacity            = sum(price_maker_candidates_df[:, "Capacity_sum"])

    percentual_capacity = 0
    C                   = size(sum_genco_capacity_df, 1)
    price_maker         = Vector{Bool}(undef, C)
    price_maker         .= false
    
    i = 1
    while percentual_capacity < 0.6
        if !is_small_gen(sum_genco_capacity_df[i, "GENCO_Name"])
            percentual_capacity = percentual_capacity + sum_genco_capacity_df[i, "Capacity_sum"] / total_capacity
            price_maker[i] = true
        end
        i = i + 1
    end

    sum_genco_capacity_df = hcat(sum_genco_capacity_df, price_maker)
    DataFrames.sort!(sum_genco_capacity_df, "GENCO")

    return rename(sum_genco_capacity_df, Dict("x1" => "Price_Maker"))
end

"Create a dictionary with the input structures."
function create_input_data(path::String, price_makers::Vector{String}; T::Int64 = 24)

    dict_input = Dict()

    @info("Creating input...")
    dict_input["bus"]      = create_bus_structs(path)
    dict_input["zone"]     = create_zone_structs(path)
    dict_input["line"]     = create_line_structs(path)
    dict_input["exchange"] = create_exchange_structs(path)
    dict_input["load"]     = create_load_structs(path)
    dict_input["thermal"]  = create_thermal_structs(path, T)
    dict_input["hydro"]    = create_hydro_structs(path)
    dict_input["company"]  = create_company_structs(path, T, price_makers)

    return dict_input
end

"Calculates the maximum travel time."
function get_maximum_travel_time(hydro::Vector{HydroGenerator})

    maximum_travel_time = 0
    J = length(hydro)

    for j in 1:J
        if length(values(hydro[j].cascade)) > 0 && (maximum(values(hydro[j].cascade)) > maximum_travel_time)
            maximum_travel_time = maximum(values(hydro[j].cascade))
        end
    end

    return maximum_travel_time
end