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

#    y0 = [0.0, df2fit.Oini[i], 0.0, 0.0, 0.0]
        y0 = [:H₂ => 0.,:O => df2fit.Oini[i], :H => 0.0, :OH => 0.0, :H₂O => 0.0]

        tspan = ( d[begin,1], d[end,1] + 200.0 ) #.- t0s[i]

        νs = [ df2fitpar[:, r][i].value for r in [:ν1,:νm1,:ν2,:ν3,:ν4,:ν5] ]
        ϵs = [ df2fitpar[:, r][i].value for r in [:ϵ1,:ϵm1,:ϵ2,:ϵ3,:ϵ4,:ϵ5] ]
#        rates = Arrhenius.( df2fit.temperature[i], νs, ϵs)
        rates = Pair.(parameters(df2fit.reaction[i]),Arrhenius.( df2fit.temperature[i], νs, ϵs))

#        prob = ODEProblem( (ydot,y,r,t) -> eqns!(ydot,y,r,t, df2fit.beampars[i], df2fit.geompars[i]) ,y0,tspan,rates )

# # solve ODEs
        prob = ODEProblem(df2fit.reaction[i], y0, tspan, rates)
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
