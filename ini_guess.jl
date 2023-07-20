function set_initial_guesses!(df::DataFrame)

    set_guess_rate!(df)    
    set_guess_Arrh!(df)
    set_guess!(df)    
    
end

"""
Adds columns to the df with initial guesses for rate for a local fit
"""
function set_guess_rate!(df::DataFrame)

    nrows = nrow(df)

    for r in rate_guesses

        k_name = rate_constants_base*r[:sfx]

        if !(r[:sfx] in rate_constants_sfx)
            error("parameter with suffix "*r[:sfx]*" does not exist")
        end

        df[!,k_name] = [fitpar() for _ in 1:nrows]

        for i=1:nrows
            df[i,k_name].value = Arrhenius(df[i,:temperature], r[:ν], r[:ϵ])
            df[i,k_name].var   = r[:var]
            df[i,k_name].min   = df[i,k_name].value/1024.0
            df[i,k_name].glbl  = false
        end
    
    end

end

"""
Adds columns to the df with initial guesses for Arrhenius prefactor and energy
for a global fit
"""
function set_guess_Arrh!(df::DataFrame)

    nrows = nrow(df)

    for r in Arrh_guesses

        kname = rate_constants_base*r[:sfx]
        νname = "ν"*r[:sfx]
        ϵname = "ϵ"*r[:sfx]
    
        if !(r[:sfx] in rate_constants_sfx)
            error("parameter with suffix "*r[:sfx]*" does not exist")
        end

        df[!,kname] = [fitpar() for _ in 1:nrows]
        df[!,νname] = [fitpar() for _ in 1:nrows]
        df[!,ϵname] = [fitpar() for _ in 1:nrows]
    
        for i=1:nrows

            df[i,kname].value = Arrhenius(df[i,:temperature], r[:ν], r[:ϵ])
            df[i,kname].glbl  = true

            df[i,νname].value = r[:ν]
            df[i,ϵname].value = r[:ϵ]
            df[i,νname].var   = r[:var]
            df[i,ϵname].var   = r[:var]
            df[i,νname].min   = df[i,νname].value/1024.0
            df[i,ϵname].min   = df[i,ϵname].value/1024.0
            df[i,νname].glbl  = true
            df[i,ϵname].glbl  = true
    
        end
    
    end

end


"""
Adds columns to the df with initial guesses for Arrhenius parameters
with detailed specifications for globality and variability

Warning: the function is too flexible, so be careful
"""
function set_guess!(df::DataFrame)

    nrows = nrow(df)

    for r in guesses

        if (!r[:glbl_ν] & !r[:glbl_ϵ] & r[:var_ν] & r[:var_ϵ])
            error("set_guess(): The combination of boolean arguments for sfx "*r[:sfx]*" does not make any sense")
        end
    
        if !(r[:sfx] in rate_constants_sfx)
            error("parameter with suffix "*r[:sfx]*" does not exist")
        end

        νname = "ν"*r[:sfx]
        ϵname = "ϵ"*r[:sfx]
    
        df[!,νname] = [fitpar() for _ in 1:nrow]
        df[!,ϵname] = [fitpar() for _ in 1:nrow]
    
        for i=1:nrows

            df[i,νname].value = r[:ν]
            df[i,ϵname].value = r[:ϵ]
            df[i,νname].var   = r[:var_ν]
            df[i,ϵname].var   = r[:var_ϵ]
            df[i,νname].min   = df[i,νname].value/1024.0
            df[i,ϵname].min   = df[i,ϵname].value/1024.0
            df[i,νname].glbl  = r[:glbl_ν]
            df[i,ϵname].glbl  = r[:glbl_ϵ]
    
        end
    
    end

end
    

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

        if r[:value] isa Number

            [df[i,r[:name]].value = r[:value] for i=1:nrows]

        elseif r[:value] isa Tuple

            if r[:value][1] == "maxs"
                [df[i,r[:name]].value = r[:value][2]*maxs[i][1] for i=1:nrows]
            elseif r[:value][1] isa Function
                [df[i,r[:name]].value = 
                    r[:value][1](r[:value][2], kinetic_traces[i], mins[i], maxs[i], δs[i])  for i=1:nrows]
            else
                error("set_guess_par(): the 2nd element of a guess_par tuple must be either 'maxs' or a function")
            end

        elseif r[:value] isa String
            
            data = [ get_results_local(output_path, df[i,:ktfname], crit=r[:value]) for i=1:nrows]
            for i=1:nrows
                if isnothing(data[i])
                    println("  ", df[i,:ktfname])
                end
            end
            
            if any(x->isnothing(x), data)
                error("Do local fit for the above traces.")
            end

            [df[i,r[:name]].value = data[i][2][r[:name]] for i=1:nrows]

        else
            error("set_guess_par(): a value must be a number, a string or a tuple{number, string}")
        end

        for i=1:nrows
            haskey(r,:min) && (df[i,r[:name]].min = r[:min])
            haskey(r,:max) && (df[i,r[:name]].max = r[:max])
            haskey(r,:var) && (df[i,r[:name]].var = r[:var])
            haskey(r,:glbl)&& (df[i,r[:name]].glbl= r[:glbl])
        end
        
    end

          
                
end
