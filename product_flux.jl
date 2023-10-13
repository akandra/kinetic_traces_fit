"""
Returns a vector of product fluxes calculated from a kinetic model

product_flux(x, p, df2fit, df2fitpar, data, ndata, flg)

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
function product_flux(x, p::Vector{Float64}, df2fit, df2fitpar, data, ndata, flg)
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

    # update rate constants for Arrhenius fitting parameters
    for (j,s) in enumerate(rate_constants_sfx)
        if df2fit[1,rate_constants_base*s].glbl
            if T_function_keys[j,1] == 1
                for i in 1:ndata
                    df2fit[i,rate_constants_base*s].value = 
                        Arrhenius(df2fit[i,:temperature], 
                                  df2fitpar[1,"ν"*s].value, df2fitpar[1,"ϵ"*s].value)
                end
            end
            if T_function_keys[j,1] == 2
                for i in 1:ndata
                    df2fit[i,rate_constants_base*s].value = 
                        hTST(df2fit[i,:temperature], 
                             [  df2fitpar[1,"hνR"*s*string(j)].value for j in 1:T_function_keys[j,2] ],
                             [ df2fitpar[1,"hνTS"*s*string(j)].value for j in 1:T_function_keys[j,2]-1 ],
                            df2fitpar[1,"ϵ"*s].value)
                end
            end
        end
    end

    flux=Float64[]
    for (i, d) in enumerate(data)

        # WARNING: think of something better than 200.0
        tspan = ( d[begin,1], d[end,1] + 200.0 ) #.- t0s[i]

        # get values for rate constants
        rates = [ df2fit[i,r].value for r in rate_constants]

        # initialize concentrations
        y0 = zeros(Float64,length(species))
        # get concentration of coverage species
        y0 = cov_species_eq!(Oini::Float64, rates::Vector{Float64}, θs::Float64, Tsurf::Float64)

#        y0[ findfirst(isequal(df2fit.cov_species[i]), species) ] =  df2fit.cov0[i]

        prob = ODEProblem( (ydot,y,r,t) -> 
                            kin_model!(ydot, y, r, t, df2fit.pumppars[i], df2fit.step_density[i], df2fit.temperature[i]), 
                                y0, tspan, rates )
        sol = solve(prob,abstol=1e-14)
 
        fpars  = [ df2fit[i,r].value for r in fit_parnames]
        append!(flux, fit_function(sol, d[:,1], df2fit.step_density[i], rates, fpars))
    end

    return flux
end
