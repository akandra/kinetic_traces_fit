# the Arrhenius initial guess function

function ini_guess!(df::DataFrame, rate_suffix::String, 
                    ν::Float64, var_ν::Bool, glbl_ν::Bool,
                    ϵ::Float64, var_ϵ::Bool, glbl_ϵ::Bool)
# set initial values for the prefactors and activation energy values
    νname = "ν"*rate_suffix 
    ϵname = "ϵ"*rate_suffix 

    if (glbl_ν & glbl_ϵ) # Arrhenius parameters regime

        for i=1:nrow(df)
            df[i,νname].value = ν
            df[i,ϵname].value = ϵ
            df[i,νname].var   = var_ν
            df[i,ϵname].var   = var_ϵ
            df[i,νname].min   = df[i,νname].value/1024.0
            df[i,ϵname].min   = df[i,ϵname].value/1024.0
            df[i,νname].glbl  = true
            df[i,ϵname].glbl  = true
        end
    
    else # rate constant regime

        for i=1:nrow(df)
            df[i,νname].value = glbl_ν ? ν : Arrhenius(df[i,:temperature], ν, ϵ)
            df[i,ϵname].value = 0.0
            df[i,νname].min   = df[i,νname].value/1024.0
            df[i,νname].var   = var_ν
            df[i,ϵname].var   = false
            df[i,νname].glbl  = glbl_ν
            df[i,ϵname].glbl  = false # irrelevant
        end
    
    end

end
    
