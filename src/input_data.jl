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

"Creates a vector of lines structs."
function create_line_structs(path::String)

    line_df = CSV.read(joinpath(path, "Lines.csv"), DataFrame)
    L       = size(line_df, 1)

    line    = Vector{Line}(undef, L)

    for l in 1:L
        line[l] = Line(line_df[l, "Line"], line_df[l, "from"], line_df[l, "to"], line_df[l, "Fmax"], line_df[l, "x"])
    end

    return line
end

"Creates a vector of exchanges structs."
function create_exchange_structs(path::String)

    exchange_df = CSV.read(joinpath(path, "LinesZ.csv"), DataFrame)
    E           = size(line_df, 1)

    exchange    = Vector{Exchange}(undef, E)

    for e in 1:E
        exchange[e] = Exchange(line_df[e, "Line"], line_df[e, "from"], line_df[e, "to"], line_df[e, "Fmax"])
    end

    return exchange
end

"Creates a vector of load structs."
function create_load_structs(path::String)

    load_df = CSV.read(joinpath(path, "Load.csv"), DataFrame)
    D       = size(load_df, 1)

    load    = Vector{Load}(undef, D)

    for d in 1:D
        load[d] = Load(load_df[d, "Bus"], load_df[d, "Zone"], load_df[d, 3:end])
    end

    return load
end

"Creates a vector of thermal structs."
function create_thermal_structs(path::String, T::Int64)

    thermal_df = CSV.read(joinpath(path, "Thermal.csv"), DataFrame)
    I          = size(thermal_df, 1)

    thermal    = Vector{Thermal}(undef, I)

    for i in 1:I
        thermal[i] = ThermalGenerator(thermal_df[i, "TGenerator"], thermal_df[i, "TGenerator_cod"], thermal_df[i, "Name"],
                                        thermal_df[i, "GENCO"], thermal_df[i, "Gtbus"], thermal_df[i, "Gtzone"], 
                                        thermal_df[i, 10:10 + T - 1], thermal_df[i, 10 + T:10 + 2T - 1],
                                        thermal_df[i, "Rup"], thermal_df[i, "Rdown"], thermal_df[i, "UVC"], thermal_df[i, "Fuel"])
    end

    return thermal
end

"Creates a vector of hydro structs."
function create_hydro_structs(path::String)

    hydro_df            = CSV.read(joinpath(path, "Hydro.csv"), DataFrame)
    cascade_df          = CSV.read(joinpath(path, "Hydro_Cascade.csv"), DataFrame)[:, 2:end]
    inflow_df           = CSV.read(joinpath(path, "Hydro_Inflow.csv"), DataFrame)
    spillage_df         = CSV.read(joinpath(path, "Hydro_SpillageHist.csv"), DataFrame)
    turbined_outflow_df = CSV.read(joinpath(path, "Hydro_TurbinedOutflowHist.csv"), DataFrame)

    J = size(hydro_df, 1)
    maximum_τ = maximum(Matrix(cascade_df)[:, 2:end])

    hydro = Vector{HydroGenerator}(undef, J)
   
    for j in 1:J

        ρ           = hydro_df[j, "p"]
        ρc          =  hydro_df[j, "pc"]
        q_max       = hydro_df[j, "Qhmax"]
        q_min       = hydro_df[j, "Qhmin"]
        water_value = hydro_df[j, "WV"]
        c           = 1/(10^6/3600)

        g_max = q_max * ρ
        g_min = q_min * ρ
        γ     = water_value * ρc / c

        dict_cascade = Dict{Int64, Int64}()
        
        cascade_jm = findall(i -> i != 0, cascade_df[:, j])

        for n in cascade_jm
            dict_cascade[n] = cascade_df[n, j]
        end

        hydro[j] = HydroGenerator(hydro_df[j, "HGenerator"], hydro_df[j, "HGenerator_cod"], hydro_df[j, "Name"],
                                    hydro_df[j, "GENCO"], hydro_df[j, "Ghbus"], hydro_df[j, "Ghzone"],
                                    ρ, ρc, g_max, g_min, q_max, q_min, hydro_df[j, "Vmax"],
                                    hydro_df[j, "Vmin"], hydro_df[j, "Vi"], Vector(turbined_outflow_df[j, 2:maximum_τ + 1]),
                                    Vector(spillage_df[j, 2:maximum_τ + 1]), hydro_df[j, "Rup"], hydro_df[j, "Rdown"],
                                    water_value, γ, inflow_df[:, j], dict_cascade)
    end

    return hydro
end
