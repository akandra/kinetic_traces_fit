function do_it()

    # create dataframe
    df = create_df(data_path, fn_delim, fn_pump_sfx, fn_cov_sfx)

    df2fit = df
    # apply data selection conditions
    for c in conds
        if c[1] isa Pair
            df2fit = filter(c[1], df2fit)
        elseif all( isa.(c,String) )
            df2fit = df2fit[df2fit.ktfname .∈ [c],:]
        elseif all( isa.(c,Int) )
            df2fit = df2fit[[c...],:]
        elseif c[1] isa UnitRange{Int}
            df2fit = df2fit[c[1],:]
        else 
            error("condition ",c," unknown")
        end
    end

    # load kinetic traces
    kinetic_traces, maxs, mins, δs = load_kinetic_traces(df2fit, data_path)
    # get a vector of baselines for data

    # set data cutoffs
    cutoff_idx = zeros(Int32, size(kinetic_traces,1))
    for (k,kt) in enumerate(kinetic_traces)
        vmax, imax = maxs[k]
        if data_cutoff_fraction == "off"
            cutoff_idx[k] = size(kt,1)
        else
            iend = findfirst(kt[imax:end,2] .< data_cutoff_fraction*vmax )
            cutoff_idx[k] = typeof(iend) == Int ? iend + imax : size(kt,1) 
        end
    end
    # add cutoff value to the data frame
    df2fit[!,:cutoff] = cutoff_idx
    
    # add fit parameters to data frame
    set_guess_rate!(df2fit, rate_guesses)
    set_guess_pars!(df2fit,  par_guesses, kinetic_traces, mins, maxs, δs)
global debug = df2fit
    process(df2fit, kinetic_traces)

    # if any parameter is global then do a fit to all the data simultaneously
    # dataset 1 is chosen, since all .glbl's are the same
    # if  any( [x[1].glbl for x in eachcol(df2fitpar)] )
    #     ksr.global_fit(df2fit,df2fitpar, kinetic_traces, results_path, iguess, 
    #                    wtd = what_to_do )
    # else # fit each dataset separately
    #     ksr.local_fit(df2fit,df2fitpar, kinetic_traces, results_path, iguess, 
    #                   wtd = what_to_do, max_T_fit = T_cutoff)
    # end

end

#NEXTTIME: adjust the kt files and the df to our new understanding the things    
# Consider for the future development
#
# 1. Put all references to the df into the ksr-module
# 2. Change the way how the kinetic traces files are named and supplied with ancillary info.
#    Consider 3 ways of doing that:
#        - create a single ancillary file
#        - create an ancillary file for each kt file
#        - put ancillary information as a header into a kt file