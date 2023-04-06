"""
Returns a vector of product fluxes calculated from a kinetic model

flux_vector(x, p, df2fit, df2fitpar, data, ndata, flg)

x is dummy argument for curve_fit function !!! Consider to get rid of it

p is a vector of fitting parameters

df2fit::DataFrame contains all info, see function create_df for details

df2fitpar::DatFrame a view of df with fitting parameters info

kt::Vector{Array} contains measured kinetic traces [dataset_index] 
                  with time(μs) in the 1st column and the signal in the 2nd one

path::String path to the output folders

iguess::Array initial guess for fitting parameters

wtd::Tuple{String} tells what to do: (action, selectors, ...) 
                    where action can be either "fit" or "analysis, and
                    selectors defines grouping (can be empty)

"""
function H2OProduction(x, p, df2fit, df2fitpar, data, ndata, flg)

    flg[1] += 1
    # update the values of variable parameters in df2fitpar
    j=0
    for d in eachcol(df2fitpar)
        if d[1].glbl
            if d[1].var
                j = j + 1
                [ d[i].value = p[j] for i in 1:ndata ]
            end
        else
            for i in 1:ndata
                if d[i].var
                    j = j + 1
                    d[i].value = p[j]
                end
            end
        end
    end

    #--------------------------------------------------------------------------
    # calculate the model function
    #--------------------------------------------------------------------------
    t0s       = [ df2fitpar.t0[i].value       for i in 1:ndata ]
    as        = [ df2fitpar.a[i].value        for i in 1:ndata ]
    baselines = [ df2fitpar.baseline[i].value for i in 1:ndata ]
    f_trs     = [ df2fitpar.f_tr[i].value     for i in 1:ndata ]
    k_vacs    = [ df2fitpar.k_vac[i].value    for i in 1:ndata ]

    flux=Float64[]
    for (i, d) in enumerate(data)

        y0 = [0.0, df2fit.Oini[i], 0.0, 0.0, 0.0]
        tspan = ( d[begin,1], d[end,1] + 200.0 ) #.- t0s[i]

        νs = [ df2fitpar[:, r][i].value for r in [:ν1,:νm1,:ν2,:ν3,:ν4,:ν5] ]
        ϵs = [ df2fitpar[:, r][i].value for r in [:ϵ1,:ϵm1,:ϵ2,:ϵ3,:ϵ4,:ϵ5] ]
        rates = Arrhenius.( df2fit.temperature[i], νs, ϵs)

        prob = ODEProblem( (ydot,y,r,t) -> eqns!(ydot,y,r,t, df2fit.beampars[i], df2fit.geompars[i]) ,y0,tspan,rates )
        sol = solve(prob,abstol=1e-14)(d[:,1] .- t0s[i])[5,:]
        
        normfactor = maximum(sol)
        if normfactor == 0
            append!(flux, as[i]*sol .+ baselines[i])
        else
    #        append!(flux, as[i]*(1.0 .+ f_trs[i]*exp.(-k_vacs[i]*(d[:,1] .- t0s[i]) ) ) .* sol/normfactor .+ baselines[i])
            sol_acc = zeros(length(sol))
            cumsum!(sol_acc, sol/normfactor)
            #append!(flux, as[i]*(sol + f_trs[i]*exp.(-k_vacs[i]*(d[:,1] .- t0s[i])) .* sol_acc )/normfactor .+ baselines[i] )
            append!(flux, as[i]*(sol/normfactor + f_trs[i]*exp.(-k_vacs[i]*(d[:,1] .- t0s[i])) .* sol_acc/sol_acc[end]) .+ baselines[i] )
        end

    end
     
    return flux
end
