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

    #ndata = ksr.nrow(df2fit)
    # load kinetic traces
    #kinetic_traces, maxs, mins, δs = ksr.load_kinetic_traces(df2fit,data_path, cutoff=-Inf)

        # # select df columns of type fitpar
        # df2fitpar = df2fit[!,names(df2fit,ksr.fitpar)]
        # # save initial guesses
        # iguess = zeros(ndata, ksr.ncol(df2fitpar))
        # for i in 1:ndata
        #     for j in 1:ksr.ncol(df2fitpar)
        #         iguess[i,j] = df2fitpar[i,j].value
        #     end
        # end

    # create geometry parameters data frame
    # dfgeom = DataFrame( facet = String[], geompars = Vector{Float64}[] )
    # push!(dfgeom, ( "332", [1.0/6.0, 2, 1, 1, 1]) )
    # push!(dfgeom, ( "111", [  0.005, 2, 1, 1, 1]) )
    #df = innerjoin(df,dfgeom, on=:facet)

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