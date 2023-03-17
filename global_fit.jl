"""
Do global fit to a set of kinetic traces

global_fit(df, dfpar, kt, path, iguess;  wtd)

df::DataFrame contains all info, see function create_df for details

dfpar::DatFrame a view of df with fitting parameters info

kt::Vector{Array} contains measured kinetic traces [dataset_index] 
                  with time(μs) in the 1st column and signal in the 2nd open

path::String path to the output folders

iguess::Array initial guess for fitting parameters

wtd::Tuple{String} tells what to do: (action, selectors, ...) 
                    where action can be either "fit" or "analysis, and
                    selectors defines grouping (can be empty)

"""
function global_fit(df2fit, df2fitpar, kinetic_traces, results_path, iguess; wtd = ("analysis",) )

    if wtd[1] == "fit"

        # number of data sets
        ndata = nrow(df2fit)

        # println("ν2=",df2fitpar[1,:ν2].value, " is "*(df2fitpar[1,:ν2].var ? "variable" : "fixed"))
        # println("ϵ2=",df2fitpar[1,:ϵ2].value, " is "*(df2fitpar[1,:ϵ2].var ? "variable" : "fixed"))
        # println("ν2=",df2fitpar[3,:ν2].value, " is "*(df2fitpar[3,:ν2].var ? "variable" : "fixed"))
        # println("ϵ2=",df2fitpar[1,:ϵ2].value, " is "*(df2fitpar[1,:ϵ2].var ? "variable" : "fixed"))
        # println("ν3=",df2fitpar[1,:ν3].value, " is "*(df2fitpar[1,:ν3].var ? "variable" : "fixed"))
        # println("ϵ3=",df2fitpar[1,:ϵ3].value, " is "*(df2fitpar[1,:ϵ3].var ? "variable" : "fixed"))
        # println("ν3=",df2fitpar[3,:ν3].value, " is "*(df2fitpar[3,:ν3].var ? "variable" : "fixed"))
        # println("ϵ3=",df2fitpar[1,:ϵ3].value, " is "*(df2fitpar[1,:ϵ3].var ? "variable" : "fixed"))
        # println("a=",df2fitpar[1,:a].value, " is "*(df2fitpar[1,:a].var ? "variable" : "fixed"))
        # println("t0=",df2fitpar[1,:t0].value, " is "*(df2fitpar[1,:t0].var ? "variable" : "fixed"))
        # println("f_tr=",df2fitpar[1,:f_tr].value, " is "*(df2fitpar[1,:f_tr].var ? "variable" : "fixed"))
        # println("k_vac=",df2fitpar[1,:k_vac].value, " is "*(df2fitpar[1,:k_vac].var ? "variable" : "fixed"))
        # println("baseline=",df2fitpar[1,:baseline].value, " is "*(df2fitpar[1,:baseline].var ? "variable" : "fixed"))

        print("Do global fit on ", ndata, " data sets...")

        # compose x and y data for curve_fit
        xdata = collect(Iterators.flatten([ kt[:,1] for kt in kinetic_traces ]))
        ydata = collect(Iterators.flatten([ kt[:,2] for kt in kinetic_traces ]))

        # compose a vector of initial values, lower and upper bounds for the fitting parameters
        pini= Float64[]
        plb = Float64[]
        pub = Float64[]
        for d in eachcol(df2fitpar)
            n = d[1].glbl ? 1 : ndata
            if d[1].var
                for i=1:n
                    push!(pini, d[1].value)
                    push!(plb , d[1].min)
                    push!(pub , d[1].max)
                end
            end
        end
        H2OProdflg = [0]
        fit = curve_fit( (x,p)->H2OProduction(x, p, df2fit, df2fitpar, kinetic_traces, ndata, H2OProdflg), 
                        xdata, ydata, pini, lower=plb, upper=pub; 
                        maxIter=1000)
        best_fit_pars = fit.param

        println(" It took ", H2OProdflg[1], " function calls.")
#        println("Fit pars: ",fit.param)

        # save the parameter-related data to a file
        open(results_path*"global_fit.dat", "a") do io
            write(io, "# file names\n")
            writedlm(io, df2fit.ktfname)
            write(io, "# initial guess\n")
            for i in 1:ndata
                writedlm(io, transpose(iguess[i,:]))
            end
            write(io,"# best fit\n")
            for i in 1:ndata
                writedlm(io, transpose([ d.value for d in df2fitpar[i,:] ]))
            end
            write(io,"# standard error\n")
            sef = try 
                    stderror(fit)
                catch e
                    fill(NaN,length(fit.param))
                end
            se = zeros( ndata, ncol(df2fitpar) )
            j = 1
            for (id,d) in enumerate(eachcol(df2fitpar))
                if d[1].var
                    if d[1].glbl
                        se[:,id] .= sef[j]
                        j+=1
                    else
                        se[:,id] = sef[j:j+ndata-1]
                        j+=ndata
                    end
                end
            end
            writedlm(io, se)
            write(io, "# reduced χ²\n")
            χ2 = sum(fit.resid .^ 2)/(length(fit.resid) - j + 1)
            writedlm(io, χ2)
        end  

    elseif wtd[1] == "analysis"
        #%% 
        data = get_results_global(results_path, df2fit.ktfname)
        best_fit_pars = Float64[]
        for k in 1:ncol(df2fitpar)
            if data[3][1,k] > 0.0 # is variable?
                if all(data[2][1,k] .== data[2][:,k]) # is global?
                    push!(best_fit_pars, data[2][1,k])
                else
                    push!(best_fit_pars, data[2][:,k]...)
                end
            end
        end
        χ2 = data[4][1]
        # update fit parameter values in df2fitpar
        for i in 1:nrow(df2fitpar)
            for j in 1:ncol(df2fitpar)
                df2fitpar[i,j].value = data[2][i,j]
            end
        end
        #%%        
    elseif wtd[1] == "clean"
        println("Tell me what and how to clean!")
    else
        println("Error: Don't know wtd")

    end

    # plotting data
    H2OProdflg = [0]
    yfit = H2OProduction(0, best_fit_pars, df2fit, df2fitpar, kinetic_traces, ndata, H2OProdflg)
    
    # find global dataset characteristics
    ch_list = ["rrr", "facet", "tag"]
    current_y = [0.95] 
    ann = [((0.5,current_y[1]), text("Global fit:"))]
    current_y = [0.9] 
    shift = 0.065
    for t in ch_list
        if length(unique(df2fit[:,t]))==1
            current_y[1] -= shift
            push!(ann, ( (0.1, current_y[1] ), 
                          text(t*"="*string(df2fit[1,t]),:left,10) ) )
        end
    end
    # find global fitting parameters
    current_y[1] -= shift
    for n in names(df2fitpar)
        if df2fitpar[1,n].glbl
            current_y[1] -= shift
            push!(ann, ( (0.1,current_y[1]), 
                  text(n*"="*string( round(df2fitpar[1,n].value, sigdigits=3) ),:left,10) ) )
        end
    end

    plots = []

    # 1st panel with global fit info
    push!(plots, plot(framestyle=:box,ticks=false,ann=ann) )

    icounter = [1]
    for (i,d) in enumerate(kinetic_traces)
        fbase= split(df2fit.ktfname[i],".")[1]
        facet = "("*split(fbase,"-")[2]*")"
        tsurf = "T = "*split(fbase,"-")[3]*" K"
        rrr = "RRR = "*split(fbase,"-")[4]*"/"*split(fbase,"-")[5]
        push!(plots, 
                plot(d[:,1], [d[:,2], yfit[icounter[1]:icounter[1]+size(d,1)-1]],
                seriestype=[:scatter :line], framestyle=:box, label=["data" "global fit"],
                xlabel="time (μs)", ylabel="H₂O flux (a.u.)", 
                title="#"*string(findfirst(df.ktfname .== df2fit[i,:].ktfname))
                 )
            )
        icounter[1] = icounter[1] + size(d,1)
    end

    page_layout = (3,4)
    nslots = prod(page_layout)
    for p in 1:nslots-mod(length(plots),nslots)
        push!(plots, plot(showaxis=false,ticks=false))
    end
    for p in 1:nslots:length(plots)
        display( plot!(plots[p:p+nslots-1]..., layout = page_layout, size=(1200,900)) )
    end

end;
