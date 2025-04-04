function do_local_fit(df2fit::DataFrame, df2fitpar::DataFrame, 
                      kinetic_traces::Vector{Matrix{Float64}}, 
                      iguess::Matrix{Float64})
        
    ndata = nrow(df2fit)

    # initialization
    fit1 = []

    # fit each data set separately
    for i in 1:ndata

        # create single row dataframes, so we can iterate over columns
        df2fit1    = df2fit[i:i,:]
        df2fitpar1 = df2fitpar[i:i,:]

        # construct x and y vectors for curve_fit and plotting
        xdata1 = kinetic_traces[i][1:df2fit1[1,:cutoff],1]
        ydata1 = kinetic_traces[i][1:df2fit1[1,:cutoff],2]

        if wtd[1] == "fit"
            # compose a vector of initial values, lower and upper bounds for the fitting parameters
            pini1= Float64[]
            plb1 = Float64[]
            pub1 = Float64[]
            for d in eachcol(df2fitpar1)
                if d[1].var
                    push!(pini1, d[1].value)
                    push!(plb1 , d[1].min)
                    push!(pub1 , d[1].max)
                end
            end

            best_fit_pars = Float64[]
            χ2 = 0.0
            if length(pini1) == 0
                println("No fitting parameters found for data set ",i,". Using initial guess.")
                
            else
                call_counts = [0]
                @time push!(fit1, 
                curve_fit( (x,p)->product_flux(x, p, df2fit1, df2fitpar1, [ kinetic_traces[i][1:df2fit1[1,:cutoff]] ], 1, call_counts), 
                        xdata1, ydata1, pini1, lower=plb1, upper=pub1; 
                        lambda = 10,min_step_quality=1e-3, maxIter=1000)
                        )
                best_fit_pars = fit1[i].param
                println("Number of function calls: ",call_counts)
                println("Fit pars: ",fit1[i].param)
                println("Done fitting date #",i," of ",ndata)

                # write out parameter-related data
                mkpath(output_path)
                open( joinpath(output_path, df2fit[i,:].ktfname), "a") do io
                    write(io, "# parameter names\n")
                    writedlm(io, hcat(names(df2fitpar)...) )
                    write(io, "# initial guess\n")
                    writedlm(io, transpose(iguess[i,:]))
                    write(io,"# best fit\n")
                    writedlm(io, transpose([ d.value for d in df2fitpar[i,:] ]))
                    write(io,"# standard error\n")
                    sef = try 
                            stderror(fit1[i])
                        catch e
                            fill(NaN,length(fit1[i].param))
                        end
                    se = zeros(ncol(df2fitpar))
                    j = 0
                    for d in 1:ncol(df2fitpar)
                        if df2fitpar[i,d].var
                            j+=1
                            se[d] = sef[j]
                        end
                    end
                    writedlm(io, transpose(se))
                    write(io, "# reduced χ²\n")
                    χ2 = sum(fit1[i].resid .^ 2)/(length(fit1[i].resid) - j)
                    writedlm(io, χ2)

                end
            end        

        elseif wtd[1] == "analysis"

            data = get_results_local(output_path,df2fit[i,:ktfname], crit="best")
            isnothing(data) && println("Error: Fit results file ", joinpath(output_path,df2fit[i,:ktfname]),
                                       " does not exist.\n     Please, run the fit for these data.")

            best_fit_pars = Float64[]
            for d in names(df2fitpar)
                if data[3][d] != 0.0
                    push!(best_fit_pars, data[2][d])
                end
            end
            χ2 = data[4][1]
            # update fit parameter structure in df2fitpar
            for p in data[2]
                df2fitpar[i,p[1]].value = p[2]
                df2fitpar[i,p[1]].var = (data[3][p[1]] != 0.0)
            end

        elseif wtd == "clean"

            println("Tell me what and how to clean!")

        else

            println("Error: Don't know wtd")

        end

        # plotting individual kinetic trace i together with the local fit to it
        call_counts = [0]
        yfit1 = product_flux(0, best_fit_pars, df2fit1, df2fitpar1, 
                    [kinetic_traces[i]],1,call_counts)

        fbase= split(df2fit.ktfname[i],".")[1]
        tag = split(fbase,"-")[1]
        facet = "("*split(fbase,"-")[2]*")"
        tsurf = "T = "*split(fbase,"-")[3]*" K"
        rrr = "RRR = "*split(fbase,"-")[4]*"/"*split(fbase,"-")[5]
        display(plot( kinetic_traces[i][:,1], [kinetic_traces[i][:,2], yfit1, [df2fit1.baseline[1].value],[xdata1[df2fit1[1,:cutoff]]]],
            seriestype=[:scatter :line :hline :vline], framestyle=:box, label=["data" "fit" "baseline" "cutoff"],
            xlabel="time (μs)", ylabel="flux (a.u.)", 
            title=tag*": "*facet*", "*tsurf*", "*rrr,
            size=(600,600),
            ylimits=(minimum(kinetic_traces[i][:,2]),maximum(kinetic_traces[i][:,2])),
            ann=append!(
                [((0.7,0.85-0.05*j), ann_par(rate_constants_base*"_{"*k*"}",  df2fit[:,rate_constants[j]][i])) for (j,k) in enumerate(rate_constants_sfx)],
                [((0.7,0.84-0.05*size(rate_constants,1)-0.05*j), ann_par(f,  df2fit[:,f][i])) for (j,f) in enumerate(fit_parnames[1:3])],
                [((0.7,0.33), text("\$\\chi^2 = "*string(round(χ2/maximum(ydata1)^2,sigdigits=3))*"\$",:left)) ]) ) )
                                        
    end

    # extract data for Arrhenius fit
    if wtd[1] == "analysis"

        df_grouped = groupby(df2fit, [wtd[2:end]...])

        color_list = [:black, :red, :blue, :green, :orange, :gray, :magenta, :cyan]

        # create a list of wtd[2] unique values
        crit_list = unique(df2fit[:,wtd[2]])
        crit_dict = Dict(crit_list .=> 1:length(crit_list))

        if T_cutoff > 0
            plot_template = plot([1/T_cutoff],framestyle=:box, 
                xlabel="1000/T (1/K)", ylabel="log(rate/μs)", seriestype=:vline, 
                label="\$ T_\\textrm{max}="*string(T_cutoff)*"\$ K", color=:gray,line=(:dash))
        else
            plot_template = plot(framestyle=:box, xlabel="1000/T (1/K)", ylabel="log(rate/μs)")
        end
        plot_dict = Dict()

        for (i,dfg) in enumerate(df_grouped)

            # get names of rate constants that are allowed to vary
            k_sfxs = [] 
            for (j,sfx) in enumerate(rate_constants_sfx)
                if dfg[1,rate_constants[j]].var
                    push!(k_sfxs, sfx)
                end
            end

            # fit local rates to Arrhenius
            Tinv = Vector{Float64}[]
            rlog = Vector{Float64}[]
            fitArrh = []
            for sfx in k_sfxs
                push!(Tinv, 1.0 ./ dfg.temperature)
                push!(rlog, [ log(d.value) for d in dfg[:,rate_constants_base*sfx] ])
                # pick up the data subset and do fit
                ipos = findall(x->x>1/T_cutoff,Tinv[end])
                push!(fitArrh, curve_fit((x,p) -> p[1] .- p[2]*x, Tinv[end][ipos], rlog[end][ipos], [1.0, 1.0]))
            end

            # do global fit for the data in ith group

            # get unbound copy of parameter part of dfg
            dfgpar = deepcopy(dfg[:,names(dfg,fitpar)])

            # adjust the values of glbl and var
            for (i,sfx) in enumerate(k_sfxs)
                # remove a rate constant column from dfgpar
                select!(dfgpar,Not(Symbol(rate_constants_base*sfx)))
                # get column names for the Arrhenius parameters                
                νname = "ν"*sfx
                ϵname = "ϵ"*sfx
                # add them to the dfgpar
                dfgpar[!,νname] = [fitpar() for _ in eachrow(dfgpar)]
                dfgpar[!,ϵname] = [fitpar() for _ in eachrow(dfgpar)]
                # fill them with Arrhenius fit values
                # and set the glbl and var fields to true
                for r in eachrow(dfgpar)
                    r[νname].glbl = true 
                    r[ϵname].glbl = true 
                    r[νname].var  = true 
                    r[ϵname].var  = true 
                    r[νname].value = exp(fitArrh[i].param[1]) 
                    r[ϵname].value = fitArrh[i].param[2]/11604.5
                end
            end

            # fix the local values for f_tr, k_vac
            for r in eachrow(dfgpar)
                r[:f_t].var  = false
                r[:k_vac].var = false
            end

            # get indices for rows with relevant kinetic traces
            kt_idx = Int[]
            for f in dfg.ktfname
                push!(kt_idx, findfirst( ==(f), df2fit.ktfname) )
            end

            # set initial guess for prefactor and energy to those from local fit
            iguessg = zeros(nrow(dfgpar), ncol(dfgpar))
            for k in 1:nrow(dfgpar)
                for j in 1:ncol(dfgpar)
                    iguessg[k,j] = dfgpar[k,j].value
                end
            end

            # Trick code to believe it's in the global fit mode

            # set wtd to "fit"
            wtd[1] = "fit"
            # set glbl to true and Arrhenius T_function_key for all variable rate constants in dfg
            for (i,sfx) in enumerate(rate_constants_sfx)
                if sfx in k_sfxs
                    dfg[1,rate_constants_base*sfx].glbl = true
                    T_function_keys[i,1] = 1
                end
            end    
            # get global fit for the group dfg
            do_global_fit(dfg, dfgpar, kinetic_traces[kt_idx], iguessg)
            # back to normal mode
            wtd[1] = "analysis"
            for (i,sfx) in enumerate(rate_constants_sfx)
                if sfx in k_sfxs
                    dfg[1,rate_constants_base*sfx].glbl = false
                    T_function_keys[i,1] = 0
                end
            end    
 
            # plotting arrheniusly

            subtitle = length(wtd) == 1 ? "" : " for "*wtd[2]*"="*string(dfg[1,wtd[2]])
            ic = mod(i-1,length(color_list)) + 1

            for(j,sfx) in enumerate(k_sfxs)

                plot_title = "\$"*rate_constants_base*"_{"*sfx*"}(T)\$"*subtitle
                iplot = crit_dict[dfg[1,wtd[2]]]
                if get(plot_dict, sfx, 0) == 0 
                    plot_dict[sfx] = deepcopy(plot_template)
                end

                legend_label = [""]
                for f in wtd[3:end]
                    legend_label[1] = legend_label[1]*string(dfg[1,f])*" "
                end

                scatter!(plot_dict[sfx], Tinv[j],rlog[j], color=color_list[ic], 
                        title = plot_title, label=legend_label[1])

                ν = "ν" .* sfx
                ϵ = "ϵ" .* sfx
                ygfit = log(dfgpar[1,ν].value) .- dfgpar[1,ϵ].value*11604.5*Tinv[j]
                ylfit = fitArrh[j].param[1] .- fitArrh[j].param[2]*Tinv[j]

                plot!(plot_dict[sfx] , Tinv[j], [ylfit], #[ygfit], 
                    title = plot_title, label=["local fit" "global fit"], color=color_list[ic],
                    line=[:solid :dash],
                    ann=[(
                        (0.05,0.6 - (plot_dict[sfx].n - 1)*0.05), 
                        text("\$E_a^{(\\ell)}="
                        *string(round(fitArrh[j].param[2]*8.61733326e-5,sigdigits=3))
                        *"\\, \\textrm{eV},\$"*" \$\\nu^{(\\ell)}="
                        *string(round(1e6*exp(fitArrh[j].param[1]),sigdigits=2))*"\\, / \\textrm{s}\$", 
                        color=color_list[ic], :8,:left)),

                        ( (0.05,0.6 - (plot_dict[sfx].n - 1)*0.05 - 0.07), 
                        text("\$E_a^{(g)}="
                        *string(round(dfgpar[1,ϵ].value,sigdigits=3))
                        *"\\, \\textrm{eV},\$"*" \$\\nu^{(g)}="
                        *string(round(1e6*dfgpar[1,ν].value,sigdigits=2))*"\\, / \\textrm{s}\$", 
                        color=color_list[ic], :8,:left) )
                        ] 
                )

            end

        end
println(sort(collect(keys(plot_dict))))
        for kv in sort(collect(keys(plot_dict)))
            display(plot_dict[kv])
        end
        #savefig(plots[1],"test.png")
    end

end;
