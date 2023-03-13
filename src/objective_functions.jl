"Create future cost function."
function create_future_cost_function(model::Ml, T::Int64, hydro::Vector{HydroGenerator}; nash::Bool = true) where {Ml}

    Jgc = length(hydro)

    if nash

        return sum(hydro[j].γ * (model[:v][T, hydro[j].number] - hydro[j].v_initial) for j = 1:Jgc)
    else

        return sum(hydro[j].γ * (model[:v_grid][T, hydro[j].number] - hydro[j].v_initial) for j = 1:Jgc) + 
                sum(hydro[j].γ * (model[:v_market][T, hydro[j].number] - hydro[j].v_initial) for j = 1:Jgc)
    end
end