"""
Adds columns to the df with initial guesses for a fitting parameter
"""
function set_guess_par!(df::DataFrame; name::String, value,
                    min::Float64   = -Inf,
                    max::Float64   = Inf,
                    var::Bool      = false,
                    glbl::Bool     = false)

    if !(name in [rate_constants; fit_parnames])
        error("guess_par(): parameter "*name*" does not exist")
    end

    if name in names(df)
        error("guess_par(): duplicate parameter name "*name)
    end

    if value isa Vector && size(value,1) != nrow(df)
      error("guess_par():  length of values vector has to be equal to number of data files.")
    end
          
    df[!,name] = [fitpar() for _ in 1:nrow(df)]

    for i=1:nrow(df)
        if value isa Vector#{Float64}
            df[i,name].value = value[i]
        elseif value isa Number
            df[i,name].value = value
        else
            error("guess_par(): a value must be a number or a vector of numbers")
        end
        df[i,name].min   = min
        df[i,name].max   = max
        df[i,name].var   = var
        df[i,name].glbl  = glbl
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
function set_guess_rate!(df::DataFrame; name::String, ν::Float64, ϵ::Float64, var::Bool)

    if !(name in rate_constants)
        error("guess_rate(): parameter "*name*" does not exist")
    end

    if name in names(df)
        error("guess_rate(): duplicate parameter name "*name)
    end

    df[!,name] = [fitpar() for _ in 1:nrow(df)]
    for i=1:nrow(df)
        df[i,name].value = Arrhenius(df[i,:temperature], ν, ϵ)
        df[i,name].var   = var
        df[i,name].min   = df[i,name].value/1024.0
        df[i,name].glbl  = false
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
    
