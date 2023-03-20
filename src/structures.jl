struct ThermalGenerator 
    number::Int64
    id::Int64 # Dessem
    name::String
    owner::Int64
    bus::Int64
    zone::Int64
    g_max::Vector{Float64}
    g_min::Vector{Float64}
    ramp_up::Float64
    ramp_down::Float64
    uvc::Float64
    fuel::String
end

struct HydroGenerator
    number::Int64
    id::Int64 # Dessem
    name::String
    owner::Int64
    bus::Int64
    zone::Int64
    ρ::Float64
    ρ_cascade::Float64
    g_max::Float64
    g_min::Float64
    q_max::Float64
    q_min::Float64
    v_max::Float64
    v_min::Float64
    v_initial::Float64
    q_initial::Vector{Float64}
    s_initial::Vector{Float64}
    ramp_up::Float64
    ramp_down::Float64
    water_value::Float64 # R$/MWh
    γ::Float64 # R$/hm3
    inflow::Vector{Float64}
    cascade::Dict{Int64, Int64} # dicionário com usinas a montante: chave é o número da hidro e valor o tempo de viagem
end

struct Bus
    number::Int64
    zone::Int64
    name::String
    deficit_cost::Float64
end

struct Zone
    number::Int64
    name::String
    deficit_cost::Float64
end

struct Line
    number::Int64
    from::Int64
    to::Int64
    f_min::Float64
    f_max::Float64
    reactance::Float64
end

struct Exchange
    number::Int64
    from::Int64
    to::Int64
    f_min::Float64
    f_max::Float64
end

struct Company
    number::Int64
    name::String
    price_maker::Bool
end

struct Load
    bus::Int64
    zone::Int64
    value::Vector{Float64}
end