using DelimitedFiles
using DifferentialEquations
using LsqFit
using Parameters
using BenchmarkTools
using DataFrames
using Plots
# using Printf
# using LaTeXStrings
# using FiniteDiff

# structure to deal with fitting parameters
@with_kw mutable struct fitpar

    value::Float64 = 0.0
    min::Float64   = -Inf
    max::Float64   = Inf
    var::Bool      = false
    glbl::Bool     = false

end

function ann_par(name, p::fitpar) 
    # s = latexstring(name,L"=",string(round(p.value,sigdigits=3)))
    s = "\$"*name*" = "*string(round(p.value,sigdigits=3))*"\$"
    #s = @sprintf "\$ %.2f \\cdot 10^1\$" significand(p.value)#*2.0^exponent(p.value)
    text(s,:left, p.var ? :black : :gray)
end

function get_data(data_path, filenames)::Vector{Matrix{Float64}}

    data = Vector{Matrix{Float64}}
    data= []
    for f in filenames
        push!(data, readdlm(data_path * f, Float64, comments=true))
    end
    return data
end

"""
Gets the parameters related to the best local fit

get_results_local(path_to_the_data::AbstractString, local_fit_filename::AbstractString)

Returns:

a tuple of initial guesses, best fit parameters, standard error and χ2 for the fit with the smallest value of χ2 
if local_fit_filename exists
otherwise nothing
"""
function get_results_local(path, filename; crit = "best")

    if isfile(path*filename)

        ini_guess = Float64[]
        best_fit  = Float64[]
        se        = Float64[]
        chi2      = Float64

        data = readdlm(path*filename,comments=true) 
        if crit == "best" 
            chi2ind = argmin(data[4:4:end,1])  # get the fit with minimal chi2
        elseif crit == "last" 
            chi2ind = lastindex(data[4:4:end,1])  # get the last fit
        else 
            println("get_results_local: crit is unknown, takeing the default value.")
        end
        ini_guess = data[1+4*(chi2ind-1),:]
        best_fit  = data[2+4*(chi2ind-1),:]
        se        = data[3+4*(chi2ind-1),:]
        chi2      = data[4+4*(chi2ind-1),1]

        return ini_guess, best_fit, se, chi2
    else
        return nothing
    end
end

function get_results_global(path, filenames)

    filename = "global_fit.dat"

    filenames = sort(filenames)
    fnames    = String[]
    ini_guess = Float64[]
    best_fit  = Float64[]
    se        = Float64[]
    chi2      = Inf64
    ind       = 0 

    
    data      = readdlm(path*filename)
    ifn = first.( Tuple.(findall(s-> s=="file",data)) ) .+ 1     # indices for file names
    iig = first.( Tuple.(findall(s-> s=="initial",data)) ) .+ 1  # indices for initial guesses
    ibf = first.( Tuple.(findall(s-> s=="best",data)) ) .+ 1     # indices for best fit parameters
    ise = first.( Tuple.(findall(s-> s=="standard",data)) ) .+ 1 # indices for standard errors
    iχ2 = first.( Tuple.(findall(s-> s=="reduced",data)) ) .+ 1  # indices for reduced χ2
    nds = iig .- ifn .- 1
    # find a dataset that matches filenames and has the smallest χ2
    for i in 1:length(ifn)
        if filenames == sort(data[ifn[i]:ifn[i]+nds[i]-1,1])
            if data[iχ2[i]] < chi2 
                chi2 = data[iχ2[i]]
                ind = i
            end
        end
    end

    ini_guess = data[iig[ind]:iig[ind]+nds[ind]-1,:]
    best_fit  = data[ibf[ind]:ibf[ind]+nds[ind]-1,:]
    se        = data[ise[ind]:ise[ind]+nds[ind]-1,:]
    chi2      = data[iχ2[ind]]

    return ini_guess, best_fit, se, chi2
end

#function Arrhenius(temperature::Float64, prefactor::Float64, energy::Float64)::Float64
#        return prefactor*exp(-11604.5*energy/temperature)
#end
function Arrhenius(temperature, prefactor, energy)
    return prefactor*exp(-11604.5*energy/temperature)
end

# this instance of the function is for the local rate guess
function ini_guess!(df::DataFrame, rate_suffix::String,  
                    ν::Float64, var_ν::Bool, glbl_ν::Bool, ϵ::Float64)
# set initial values for the prefactors and activation energy
# considering ν as thermal (glbl_ν=false) or temperature-independent (glbl_ν=true) rate constant 
# (implying ϵ = 0)
    νname = "ν"*rate_suffix 
    ϵname = "ϵ"*rate_suffix 
    for i=1:nrow(df)
        df[i,νname].value = glbl_ν ? ν : Arrhenius(df[i,:temperature], ν, ϵ)
        df[i,ϵname].value = 0.0
        df[i,νname].min   = df[i,νname].value/1024.0
        df[i,νname].var   = var_ν
        df[i,ϵname].var   = false
        df[i,νname].glbl  = glbl_ν
    end
end

# this instance of the function is for the Arrhenius guess
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
        df[i,νname].glbl  = glbl_ν
        df[i,ϵname].glbl  = glbl_ϵ
    end
end
    
function H2Pulse(t::Float64, a::Vector{Float64}, fwhm::Vector{Float64}, tcenter::Vector{Float64})

    sigma = fwhm / ( 2.0*sqrt(2.0*log(2.0)) )
    
    gausses = 0
    for i=1:3
        gausses = gausses + a[i]/(sqrt(2*π)*sigma[i])*exp( -(t - tcenter[i])^2/(2*sigma[i]^2) )
    end

    return gausses
end

function get_beampars(fnames)
# ALARM!!!!! Fix me
return  [
    [[0.000120897, 0.000170694, 0.000105154],   [30.3872, 15.1266, 75.293],  [51.9994, 52.4238, 85.3074]], 
    [[0.0000148181, 0.000212279, 0.000128483],  [14.8624, 18.5067, 68.9036], [35.5675, 52.2248, 72.4033]], 
    [[0.000256086, 0.0000951222, 0.000241312],  [22.4236, 12.2085, 65.8626], [47.6367, 49.1664, 71.6434]], 
    # the weird 4th beam is replaced by the regular 6th beam
    # [[0.0000402161, 0.000242335, 0.000130398],  [11.7742, 20.1696, 177.909], [82.8323, 97.182, 96.105]], 
    [[0.0000754665, 0.0000707978, 0.000257094], [196.966, 16.4958, 25.2458], [156.164, 100.393, 111.913]],
    [[0.0000645202, 0.000256638, 0.0000599968], [11.6545, 22.2035, 52.3407], [100.279, 112.345, 123.511]], 
    [[0.0000754665, 0.0000707978, 0.000257094], [196.966, 16.4958, 25.2458], [156.164, 100.393, 111.913]]
        ]
end

function eqns!(ydot,y,p,t, beampars, geompars)

    # All-step model from Theo (everything occurs on steps)
    # (1)  H₂ + O <- km1, k1 -> HO + H
    # (2)  H + OH -k2> H₂O 
    # (3) OH + OH -k3> H₂O + O
    # (4)  H + H  -k4> H₂ 
    # (2)  H₂O -k5> H₂O(gas)  

    θ_s, σO, σOH, σH2O, σH  = geompars
    κ1, κm1, κ2, κ3, κ4, κ5 = p
    a, fwhm, tcenter        = beampars
    
    ydot[1] = -κ1*y[1]*y[2]*θ_s/σO + θ_s*κ4*y[3]*y[3]/(σH*σH) + θ_s*κm1*y[3]*y[4]/(σH*σOH) + θ_s*H2Pulse(t,a,fwhm,tcenter)
    ydot[2] = -κ1*y[1]*y[2] + κm1*y[3]*y[4]*σO/(σH*σOH) +  κ3*y[4]*y[4]*σO/(σOH*σOH)
    ydot[3] =  κ1*y[1]*y[2]/σO - 2.0*κ4*y[3]*y[3]/σH - κ2*y[3]*y[4]/σOH - κm1*y[3]*y[4]/σOH
    ydot[4] =  κ1*y[1]*y[2]/σO - κ2*y[3]*y[4]/σH - κm1*y[3]*y[4]/σH - 2.0*κ3*y[4]*y[4]/σOH
    ydot[5] =  κ2*y[3]*y[4]*σH2O/(σH*σOH) + κ3*y[4]*y[4]*σH2O/(σOH*σOH) - κ5*y[5]   

end

function fpopsol(model,u0,tspan,p)
    prob = ODEProblem(model,u0,tspan,p)
    solve(prob,Rosenbrock23())
end;

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
end;

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

function create_df(data_path)

    datafilenames = readdir(data_path)
    beamfnames      = filter(x-> occursin("beam",x) ,datafilenames)
    Oinifnames      = filter(x-> occursin("Oini",x) ,datafilenames)
    kintracesfnames = filter(x->!(occursin("beam",x) || occursin("Oini",x)) ,datafilenames)

    # create kinetic traces data frame
    tags =[ map( f->split( splitext(f)[1], "-")[i], kintracesfnames ) for i in 1:5 ]
    df = DataFrame(tags, [:tag,:facet,:temperature,:rrH2,:rrO2])
    df[!,:temperature] = parse.(Float64,df[!,:temperature])
    df[!,:rrH2] = parse.(Float64,df[!,:rrH2])
    df[!,:rrO2] = parse.(Float64,df[!,:rrO2])
    df[!,:rrr]  = df[!,:rrO2] ./ df[!,:rrH2]
    df[!,:ktfname] = kintracesfnames
    df[!,:beamfname] = tags[1] .* "-beam.dat"

    # create H2-pulse data frame
    dfbeam = DataFrame( beamfname = beamfnames, beampars = get_beampars(beamfnames) )
    df = innerjoin(df,dfbeam, on=:beamfname)

    # create geometry parameters data frame
    dfgeom = DataFrame( facet = String[], geompars = Vector{Float64}[] )
    push!(dfgeom, ( "332", [1.0/6.0, 2, 1, 1, 1]) )
    push!(dfgeom, ( "111", [  0.005, 2, 1, 1, 1]) )
    df = innerjoin(df,dfgeom, on=:facet)

    # create [O]_ini data frame
    Oinidata = get_data(data_path, Oinifnames)
    tagsOini = map( f->split( splitext(f)[1], "-")[1], Oinifnames )
    tags1 = vcat(fill.(tagsOini,size.(Oinidata,1))...)
    dfOini = DataFrame(vcat(Oinidata...),[:temperature,:rrH2,:rrO2,:Oini])
    dfOini[!,:tag] = tags1
    # join above data frames
    df = innerjoin(df, dfOini, on = [:tag, :temperature, :rrH2, :rrO2])

    # set the names for fit parameter columns
    fitparsnames = ["ν1","ϵ1","νm1","ϵm1","ν2","ϵ2","ν3","ϵ3","ν4","ϵ4","ν5","ϵ5",
                     "a", "t0", "baseline", "f_tr", "k_vac" ]
    # we avoid using fill(), see rebind_vs_mutate.jl in juliaFun to find out why
    [ df[!,n] = [fitpar() for _ in 1:nrow(df)] for n in fitparsnames]

    return df

end;

function load_kinetic_traces(df2fit)

    println("Loading Kinetic Traces...")
    [ println(" ",filename) for filename in df2fit.ktfname ]
    kinetic_traces = get_data(data_path, df2fit.ktfname)
    ndata = length(kinetic_traces)
    maxs = [ findmax(kt[:,2]) for kt in kinetic_traces ]
    mins = [ findmin(kinetic_traces[i][1:maxs[i][2],2]) for i in 1:ndata ]
    δs   = [ maxs[i][1] - mins[i][1] for i=1:ndata ]

    return kinetic_traces, maxs, mins, δs

end;

# ==========================================================================
# Main part
# ==========================================================================

# set paths to the data and results folders
path         = "../../Dropbox/Kinetics of Surface Reactions/"
data_path    = path * "data/"
results_path = path * "results/" * "nu3_1/"

# create dataframe
df = create_df(data_path)

# select kinetic trace data to fit
ktfnames = [
            "20170828-332-373-100-200.dat", 
             "20170828-332-423-100-200.dat", 
            "20170828-332-473-100-200.dat",
            "20201117-111-373-25-50.dat", 
            "20201117-111-423-25-50.dat", 
            "20201117-111-473-25-50.dat"
           ]

condition = df.ktfname .∈ [ktfnames[4:6]]
#df2fit = filter( :ktfname => in(ktfnames),df)
#condition = (df.rrr .== 2) # .| (df.rrr .== 8.0)
#condition = (df.tag .!= "20201103") .| (df.tag .!= "20201118") 
#df2fit = filter([:facet,:rrH2,:rrO2] => (f,h,o) ->  f=="332" && h==100 && o==100 ,df)
# dataframe with data to fit
#condition = 55:55
df2fit = df[condition,:]
ndata = nrow(df2fit)

# load kinetic traces
kinetic_traces, maxs, mins, δs = load_kinetic_traces(df2fit)

# initial guesses for the Arrhenius parameters cooked by Theo
ν_theoguess = [ 1.0*10^5, 1.0*10^5, 1.0*10^5, 3.0*10^4, 1.0*10^6, 1.0*10^10 ] # μs
ϵ_theoguess = [ 0.2, 2.0, 0.4, 0.36, 0.75, 0.5 ] # eV

# Set what is to be done
# fit_is_local = true or false
# what_to_do   = "fit" or "analysis"

fit_is_local = true
what_to_do = ("fit",      "rrr", "facet", "tag")
what_to_do = ("analysis", "rrr", "facet", "tag")
T_cutoff = 480.0 # max temperature for Arrhenius fit

# set initial values for the fitting parameters and other defaults
# units: μs⁻¹ for prefactors and rates; eV for energy
#

if fit_is_local
# -------------------------------------------------------------------------------
# LOCAL FIT
# -------------------------------------------------------------------------------

    # rates
    ini_guess!(df2fit, "1", ν_theoguess[1], false, false, ϵ_theoguess[1])
    ini_guess!(df2fit,"m1", ν_theoguess[2], false, false, ϵ_theoguess[2])
    ini_guess!(df2fit, "2", ν_theoguess[3], false, false, ϵ_theoguess[3])
    ini_guess!(df2fit, "3", ν_theoguess[4],  true, false, ϵ_theoguess[4])
    ini_guess!(df2fit, "4", ν_theoguess[5], false, false, ϵ_theoguess[5])
    ini_guess!(df2fit, "5", ν_theoguess[6], false, false, ϵ_theoguess[6])

    # amplitudes, t_0 etc.
    for i=1:ndata

        df2fit.a[i].value = maxs[i][1]*0.25
        df2fit.a[i].min   = 0.001
        df2fit.a[i].var   = true
    
        df2fit.t0[i].value = -100.0
        df2fit.t0[i].min   = -200.0
        df2fit.t0[i].max   =  200.0
        df2fit.t0[i].var   = true
    
        df2fit.f_tr[i].value = 0.001
        df2fit.f_tr[i].min   = 0.0
        df2fit.f_tr[i].max   = Inf
        df2fit.f_tr[i].glbl  = false
        df2fit.f_tr[i].var   = true
    
        df2fit.k_vac[i].value = 1e-5
        df2fit.k_vac[i].glbl  = false
        df2fit.k_vac[i].var   = false
    end
    
else
# -------------------------------------------------------------------------------
# GLOBAL FIT
# -------------------------------------------------------------------------------

    # prefactors and activation energies
    ini_guess!(df2fit, "1", ν_theoguess[1], false, false, ϵ_theoguess[1], false, false)
    ini_guess!(df2fit,"m1", ν_theoguess[2], false, false, ϵ_theoguess[2], false, false)
    ini_guess!(df2fit, "2", ν_theoguess[3], false, false, ϵ_theoguess[3], false, false)
    ini_guess!(df2fit, "3", ν_theoguess[4],  true,  true, ϵ_theoguess[4],  true,  true)
    ini_guess!(df2fit, "3", 4.4204050640348025e7,  true,  true, 0.5327034287892817,  true,  true)
    ini_guess!(df2fit, "4", ν_theoguess[5], false, false, ϵ_theoguess[5], false, false)
    ini_guess!(df2fit, "5", ν_theoguess[6], false, false, ϵ_theoguess[6], false, false)

    # amplitudes, t_0s etc.
    for i=1:ndata

        # get initial guesses from local fit if it exists
        data = get_results_local(results_path, df2fit[i,:ktfname],crit="best")
        if !isnothing(data)
            # put the values into the dataframe
            for (j,p) in enumerate(names(df2fit,fitpar))
                if p == "a"
                    df2fit[i,p].value = data[2][j]
                    df2fit[i,p].min   = 0.001
                    df2fit[i,p].var   = true
                elseif p == "t0"
                    df2fit[i,p].value = data[2][j]
                    df2fit[i,p].var   = true
                elseif (p == "f_tr") || (p == "k_vac")
                    df2fit[i,p].value = data[2][j]
                    df2fit[i,p].var   = false
                end
            end
        # or set initial guesses manually
        else
            df2fit.a[i].value = maxs[i][1]*0.25
            df2fit.a[i].min   = 0.001
            df2fit.a[i].var   = true
        
            df2fit.t0[i].value = -100.0
            df2fit.t0[i].min   = -200.0
            df2fit.t0[i].max   =  200.0
            df2fit.t0[i].var   = true
        
            df2fit.f_tr[i].value = 0.001
            df2fit.f_tr[i].min   = 0.0
            df2fit.f_tr[i].max   = Inf
            df2fit.f_tr[i].var   = true
        
            df2fit.k_vac[i].value = 1e-5
            df2fit.k_vac[i].var   = false
        end

    end
    
end

# set initial values for baselines
for i=1:ndata
    bl_range = findfirst( x -> x>mins[i][1]+δs[i]/2 ,kinetic_traces[i][1:maxs[i][2],2]) - 7
    df2fit.baseline[i].value = mean(kinetic_traces[i][1:bl_range,2])
    #df2fit.baseline[i].value = 0.0
end

# set data cutoffs
cutoff = zeros(Int32, ndata)
for (k,kt) in enumerate(kinetic_traces)
    vmax, imax = maxs[k]
    iend = findfirst(kt[imax:end,2] .< -10.0*vmax )
    cutoff[k] = typeof(iend) == Int ? iend + imax : size(kt,1) 
end
df2fit[!,:cutoff] = cutoff

# select df columns of type fitpar
df2fitpar = df2fit[!,names(df2fit,fitpar)]
# save initial guesses
iguess = zeros(ndata, ncol(df2fitpar))
for i in 1:ndata
    for j in 1:ncol(df2fitpar)
        iguess[i,j] = df2fitpar[i,j].value
    end
end


# if any parameter is global then do a fit to all the data simultaneously
# dataset 1 is chosen, since all .glbl's are the same

if  any( [x[1].glbl for x in eachcol(df2fitpar)] )
    global_fit(df2fit,df2fitpar, kinetic_traces, results_path, iguess, 
              wtd = what_to_do )
else # fit each dataset separately
    local_fit(df2fit,df2fitpar, kinetic_traces, results_path, iguess, 
              wtd = what_to_do, max_T_fit = T_cutoff)
end