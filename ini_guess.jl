# the Arrhenius initial guess function
ksr.ini_guess!(df2fit, (νname,ϵname), fill(ν,n), true,  true, ϵ,  true, true)

ksr.ini_guess!(df2fit, "1", Arrh(ν,ϵ,T), true, false, 0, false, true)

ksr.ini_guess!(df2fit, "1",        ν, false, true, ϵ, false, true)

function ini_guess!(df::DataFrame, rate_suffix::String, 
                    ν::Float64, var_ν::Bool, glbl_ν::Bool,
                    ϵ::Float64, var_ϵ::Bool, glbl_ϵ::Bool)
# set initial values for the prefactors and activation energy values
    νname = "ν"*rate_suffix 
    ϵname = "ϵ"*rate_suffix 

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
    

end
    
