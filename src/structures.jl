struct ThermalGenerator 
    id::Int64 # Dessem
    name::String
    owner::Int64
    bus::Int64
    zone::Int64
    g_max::Float64
    g_min::Float64
    ramp_up::Float64
    ramp_down::Float64
    uvc::Float64
    fuel::String
end

struct SolarGenerator
end

struct WindGenerator
end

struct HydroGenerator
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
    s_max::Float64
    v_initial::Float64
    q_initial::Vector{Float64}
    s_initial::Vector{Float64}
    ramp_up::Float64
    ramp_down::Float64
    water_value::Float64 # R$/MWh
    γ::Float64 # R$/hm3
end

struct Bus
end

struct Zone
end

struct Line
end

struct Exchange
end

struct Company
end
