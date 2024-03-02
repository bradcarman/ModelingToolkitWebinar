function clean(eq)
    eqs = string(eq)
    eqs = replace(eqs, "Differential(t)"=>"der")
    eqs = replace(eqs, "(t)"=>"")
    eqs = replace(eqs, "ρ"=>"rho")
    eqs = replace(eqs, "β"=>"beta")
    eqs = replace(eqs, "₁"=>"1")
    eqs = replace(eqs, "₂"=>"2")
    eqs = replace(eqs, "₊"=>"_")
    eqs = replace(eqs, "ẋ"=>"dx")
    eqs = replace(eqs, "ṙ"=>"dr")
    eqs = replace(eqs, "ṁ"=>"dm")
    eqs = replace(eqs, "ₒ"=>"_0")
    eqs = replace(eqs, "₀"=>"_0")
    
    
    return eqs
end

function convert_to_modelica(sys::ODESystem, file="modelica.mo")

    defs = ModelingToolkit.defaults(sys)
    # guesses = ModelingToolkit.guesses(sys)
    # defs = merge(guesses, defs)
    vars = states(sys)
    eqs = full_equations(sys)
    pars = parameters(sys)

    code = String[]
    push!(code, "model MTK")
    for par in pars
        push!(code, "\tparameter Real $(clean(par)) = $(ModelingToolkit.fixpoint_sub(par, defs));")
    end

    for st in vars
        push!(code, "\tReal $(clean(st))(start = $(ModelingToolkit.fixpoint_sub(st, defs)));")
    end

    push!(code, "equation")
    for eq in eqs
        push!(code, "\t$(clean(eq.lhs)) = $(clean(eq.rhs));")
    end
    push!(code, "end MTK;")

    open(file, "w") do io
        for line in code
            println(io, line)
        end
    end
    

    return
end