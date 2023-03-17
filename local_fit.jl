function local_fit(df2fit, df2fitpar, kinetic_traces, results_path, iguess; 
                   wtd = ("analysis"), max_T_fit = 10000.0 )
                       
    # initialization
    fit1 = []
    ndata = nrow(df2fit)
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
            H2OProdflg = [0]
            @time push!(fit1, 
                curve_fit( (x,p)->H2OProduction(x, p, df2fit1, df2fitpar1, [ kinetic_traces[i][1:df2fit1[1,:cutoff]] ], 1,H2OProdflg), 
                        xdata1, ydata1, pini1, lower=plb1, upper=pub1; 
                        lambda = 10,min_step_quality=1e-3, maxIter=1000))
            best_fit_pars = fit1[i].param
            println("Number of function calls: ",H2OProdflg)
            println("Fit pars: ",fit1[i].param)
            println("Done fitting date #",i," of ",ndata)


            # write out parameter-related data
            open(results_path*df2fit[i,:].ktfname, "a") do io
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

        elseif wtd[1] == "analysis"
            #%% 
            data = get_results_local(results_path, df2fit[i,:ktfname],crit="best")
            isnothing(data) && println("Error: Fit results file ", results_path*df2fit[i,:ktfname],
                                       " does not exist.\n     Please, run the fit for these data.")
            best_fit_pars = Float64[]
            for (k, se) in enumerate(data[3])
                if se > 0.0
                    push!(best_fit_pars, data[2][k])
                end
            end
            χ2 = data[4][1]
            # update fit parameter structure in df2fitpar
            for j in 1:ncol(df2fitpar)
                df2fitpar[i,j].value = data[2][j]
                df2fitpar[i,j].var = data[3][j] > 0.0
            end
            #%%        

        elseif wtd == "clean"
            println("Tell me what and how to clean!")
        else
            println("Error: Don't know wtd")
        end
        
        # plotting individual kinetic trace i together with the local fit to it
        H2OProdflg = [0]
        yfit1 = H2OProduction(0, best_fit_pars, df2fit1, df2fitpar1, 
                             [kinetic_traces[i]],1,H2OProdflg)
        
        fbase= split(df2fit.ktfname[i],".")[1]
        facet = "("*split(fbase,"-")[2]*")"
        tsurf = "T = "*split(fbase,"-")[3]*" K"
        rrr = "RRR = "*split(fbase,"-")[4]*"/"*split(fbase,"-")[5]
        display(plot( kinetic_traces[i][:,1], [kinetic_traces[i][:,2], yfit1, [df2fit1.baseline[1].value],[xdata1[df2fit1[1,:cutoff]]]],
                    seriestype=[:scatter :line :hline :vline], framestyle=:box, label=["data" "fit" "baseline" "cutoff"],
                    xlabel="time (μs)", ylabel="H₂O flux (a.u.)", 
                    title="#"*string(findfirst(df.ktfname .== df2fit[i,:].ktfname))*": "*facet*", "*tsurf*", "*rrr,
                    size=(600,600),
                    ann=[((0.7,0.80),ann_par("\\kappa_1",  df2fit[:,:ν1  ][i])),
                        ((0.7,0.75),ann_par("\\kappa_{-1}",df2fit[:,:νm1 ][i])),
                        ((0.7,0.70),ann_par("\\kappa_2",   df2fit[:,:ν2  ][i])),
                        ((0.7,0.65),ann_par("\\kappa_3",   df2fit[:,:ν3  ][i])),
                        ((0.7,0.60),ann_par("\\kappa_4",   df2fit[:,:ν4  ][i])),
                        ((0.7,0.55),ann_par("a",           df2fit[:,:a   ][i])),
                        ((0.7,0.50),ann_par("t_0",         df2fit[:,:t0  ][i])),
                        ((0.7,0.45),ann_par("f_{tr}",      df2fit[:,:f_tr][i])),
                        ((0.7,0.4), text("\$\\chi^2 = "*string(round(χ2/maximum(ydata1)^2,sigdigits=3))*"\$",:left)) ] ) )
                                                
    end

    # extract data for Arrhenius fit
    if wtd[1] == "analysis"

        df_grouped = groupby(df2fit, [wtd[2:end]...])

        color_list = [:black, :red, :blue, :green, :orange, :gray, :magenta, :cyan]


        crit_list = unique(df2fit[:,wtd[2]])
        crit_dict = Dict(crit_list .=> 1:length(crit_list))

        if max_T_fit > 0
            plot_template = plot([1/max_T_fit],framestyle=:box, 
                    xlabel="1000/T (1/K)", ylabel="log(rate/μs)", seriestype=:vline, 
                    label="\$ T_\\textrm{max}="*string(max_T_fit)*"\$ K", color=:gray,line=(:dash))
        else
            plot_template = plot(framestyle=:box, xlabel="1000/T (1/K)", ylabel="log(rate/μs)")
        end
        plots = [ deepcopy(plot_template) for i in 1:length(crit_list)]

        for (i,dfg) in enumerate(df_grouped)

            # get names of prefactors
            ν_names = [] 
            for n in names(dfg[1,r"ν"])
                if dfg[1,n].var
                    push!(ν_names, n)
                end
            end
            
            # fit local rates to Arrhenius
            Tinv = Vector{Float64}[]
            rlog = Vector{Float64}[]
            fitArrh = []
            for ν in ν_names
                push!(Tinv, 1.0 ./ dfg.temperature)
                push!(rlog, [ log(d.value) for d in dfg[:,ν] ])
                # pick up the data subset and do fit
                ipos = findall(x->x>1/max_T_fit,Tinv[end])
                push!(fitArrh, curve_fit((x,p) -> p[1] .- p[2]*x, Tinv[end][ipos], rlog[end][ipos], [1.0, 1.0]))
            end

            # do global fit for the data in i-th group

            # get unbound copy of parameter part of dfg
            dfgpar = deepcopy(dfg[:,names(dfg,fitpar)])

            # adjust the values of glbl and var
            for (i,ν) in enumerate(ν_names) 
                ϵ = replace(ν,"ν"=>"ϵ") 
                for r in eachrow(dfgpar[:, [ν, ϵ] ])
                    r[1].glbl = true
                    r[2].glbl = true
                    r[2].var =  true 
                    r[1].value = exp(fitArrh[i].param[1]) 
                    r[2].value = fitArrh[i].param[2]/11604.5
                end
            end
            # fix the local values for t0, f_tr, k_vac
            for r in eachrow(dfgpar)
                # r[:t_0].var   = false
                r[:f_tr].var  = false
                r[:k_vac].var = false
            end

            # select relevant kinetic traces
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
            
            # get global fit for the group dfg
            global_fit(dfg, dfgpar, kinetic_traces[kt_idx], results_path, 
                        iguessg, wtd = ("fit", ))
            
            # plotting arrheniusly

            subtitle = length(wtd) == 1 ? "" : " for "*wtd[2]*"="*string(dfg[1,wtd[2]])
            ic = mod(i-1,length(color_list)) + 1

            for(j,ν) in enumerate(ν_names)

                plot_title = replace(ν,"ν"=>"κ")*"(T)"*subtitle
                iplot = crit_dict[dfg[1,wtd[2]]]
        
                legend_label = [""]
                for f in wtd[3:end]
                    legend_label[1] = legend_label[1]*string(dfg[1,f])*" "
                end

                scatter!(plots[iplot], Tinv[j],rlog[j], color=color_list[ic], 
                         title = plot_title, label=legend_label[1])

                ϵ = replace(ν,"ν"=>"ϵ")
                ygfit = log(dfgpar[1,ν].value) .- dfgpar[1,ϵ].value*11604.5*Tinv[j]
                ylfit = fitArrh[j].param[1] .- fitArrh[j].param[2]*Tinv[j]

                plot!(plots[iplot], Tinv[j], [ylfit, ygfit], 
                    title = plot_title, label=["local fit" "global fit"], color=color_list[ic],
                    line=[:solid :dash],
                    ann=[(
                        (0.05,0.6 - (plots[iplot].n - 1)*0.05), 
                        text("\$E_a^{(\\ell)}="
                        *string(round(fitArrh[j].param[2]*8.61733326e-5,sigdigits=3))
                        *"\\, \\textrm{eV},\$"*" \$\\nu^{(\\ell)}="
                        *string(round(1e6*exp(fitArrh[j].param[1]),sigdigits=2))*"\\, / \\textrm{s}\$", 
                        color=color_list[ic], :8,:left)),

                        ( (0.05,0.6 - (plots[iplot].n - 1)*0.05 - 0.07), 
                        text("\$E_a^{(g)}="
                        *string(round(dfgpar[1,ϵ].value,sigdigits=3))
                        *"\\, \\textrm{eV},\$"*" \$\\nu^{(g)}="
                        *string(round(1e6*dfgpar[1,ν].value,sigdigits=2))*"\\, / \\textrm{s}\$", 
                        color=color_list[ic], :8,:left) )
                        ] 
                )
            
            end

        end
        
        for p in plots
            display(p)
        end
#        savefig(plots[1],"test.png")
    end

end;
