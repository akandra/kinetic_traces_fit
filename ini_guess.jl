"""
Adds columns to the df with initial guesses for a fitting parameter
"""
function guess_par(df::DataFrame, name::String; val)

end

"""
Adds columns to the df with initial guesses for Arrhenius prefactor and energy
for a global fit
"""
function guess_Arrh!(df::DataFrame, kname::String, ν::Float64, ϵ::Float64, var::Bool)

    guess_rate!(df, kname,  ν, ϵ,  true)

    νname = "ν_"*kname
    ϵname = "ϵ_"*kname

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
function guess_rate!(df::DataFrame, kname::String, ν::Float64, ϵ::Float64, var::Bool)

    if !(kname in keys(rate_constants))
        error("guess_rate(): parameter "*kname*" does not exist")
    end

    if kname in names(df)
        error("guess_rate(): duplicate parameter name "*kname)
    end

    df[!,kname] = [fitpar() for _ in 1:nrow(df)]
    for i=1:nrow(df)
        df[i,kname].value = Arrhenius(df[i,:temperature], ν, ϵ)
        df[i,kname].var   = var
        df[i,kname].min   = df[i,kname].value/1024.0
        df[i,kname].glbl  = false
    end
    
end


"""
Adds columns to the df with initial guesses for Arrhenius parameters
with detailed specifications for globality and variability

Warning: the function is too flexible, so be careful
"""
function guess!(df::DataFrame, 
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
    
