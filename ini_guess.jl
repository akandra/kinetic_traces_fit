"""
Adds columns to the df with initial guesses for a fitting parameter
"""
function set_guess_pars!(df::DataFrame, par_guesses, kinetic_traces, mins, maxs, δs)

    nrows = nrow(df)

    for r in par_guesses

        if !(r[:name] in fit_parnames)
            error("parameter "*r[:name]*" does not exist")
        end

        df[!,r[:name]] = [fitpar() for _ in 1:nrows]

        for i=1:nrows
            if r[:value] isa Number
                df[i,r[:name]].value = r[:value]
            elseif r[:value] isa Tuple
                if r[:value][1] == "maxs"
                    df[i,r[:name]].value = r[:value][2]*maxs[i][1]
                elseif r[:value][1] isa Function
                    df[i,r[:name]].value = 
                        r[:value][1](r[:value][2], kinetic_traces[i], mins[i], maxs[i], δs[i])
                else
                    error("set_guess_par(): the 2nd element of a tuple must be either maxs or function")
                end
            else
                error("set_guess_par(): a value must be a number or a tuple{number, string}")
            end

            haskey(r,:min) && (df[i,r[:name]].min = r[:min])
            haskey(r,:max) && (df[i,r[:name]].max = r[:max])
            haskey(r,:var) && (df[i,r[:name]].var = r[:var])
            haskey(r,:glbl)&& (df[i,r[:name]].glbl= r[:glbl])
        end
        
    end

          
                
end

"""
Adds columns to the df with initial guesses for Arrhenius prefactor and energy
for a global fit
"""
function set_guess_Arrh!(df::DataFrame; name::String, ν::Float64, ϵ::Float64, var::Bool)

    guess_rate!(df, name,  ν, ϵ,  true)

    νname = "ν_"*name
    ϵname = "ϵ_"*name

    df[!,νname] = [fitpar() for _ in 1:nrow(df)]
    df[!,ϵname] = [fitpar() for _ in 1:nrow(df)]

    for i=1:nrow(df)
        df[i,νname].value = ν
        df[i,ϵname].value = ϵ
        df[i,νname].var   = var
        df[i,ϵname].var   = var
        df[i,νname].min   = df[i,νname].value/1024.0
        df[i,ϵname].min   = df[i,ϵname].value/1024.0
        df[i,νname].glbl  = true
        df[i,ϵname].glbl  = true
    end

end

"""
Adds columns to the df with initial guesses for rate for a local fit
"""
function set_guess_rate!(df::DataFrame, rate_guesses)

    nrows = nrow(df)

    for r in rate_guesses

        if !(r[:name] in rate_constants)
            error("parameter "*r[:name]*" does not exist")
        end

        df[!,r[:name]] = [fitpar() for _ in 1:nrows]

        for i=1:nrows
            df[i,r[:name]].value = Arrhenius(df[i,:temperature], r[:ν], r[:ϵ])
            df[i,r[:name]].var   = r[:var]
            df[i,r[:name]].min   = df[i,r[:name]].value/1024.0
            df[i,r[:name]].glbl  = false
        end
    
    end

end


"""
Adds columns to the df with initial guesses for Arrhenius parameters
with detailed specifications for globality and variability

Warning: the function is too flexible, so be careful
"""
function set_guess!(df::DataFrame, 
                    νname::String, ν::Float64, var_ν::Bool, glbl_ν::Bool,
                    ϵname::String, ϵ::Float64, var_ϵ::Bool, glbl_ϵ::Bool)

    if (!glbl_ν & !glbl_ϵ & var_ν & var_ϵ)
        error("guess(): The combination of boolean arguments does not make any sense")
    end

    if νname in names(df)
        error("guess(): duplicate parameter name "*νname)
    end
    if ϵname in names(df)
        error("guess(): duplicate parameter name "*ϵname)
    end

    df[!,νname] = [fitpar() for _ in 1:nrow(df)]
    df[!,ϵname] = [fitpar() for _ in 1:nrow(df)]

    for i=1:nrow(df)
        df[i,νname].value = ν
        df[i,ϵname].value = ϵ
        df[i,νname].var   = var_ν
        df[i,ϵname].var   = var_ϵ
        df[i,νname].min   = df[i,νname].value/1024.0
        df[i,ϵname].min   = df[i,ϵname].value/1024.0
        df[i,νname].glbl  = glbl_ν
        df[i,ϵname].glbl  = glbl_ϵ
    end
    

end
    
