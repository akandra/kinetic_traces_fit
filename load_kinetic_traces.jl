function load_kinetic_traces(df2fit, data_path; cutoff::Number)

    println("Loading Kinetic Traces...")
    [ println(" ",filename) for filename in df2fit.ktfname ]

    kinetic_traces = get_data(data_path, df2fit.ktfname)
    ndata = length(kinetic_traces)

    # finding maximum intensity of a signal
    maxs = [ findmax(kt[:,2]) for kt in kinetic_traces ]
    
    # finding minimum intensity at the l.h.s. from a maximum
    mins = [ findmin(kinetic_traces[i][1:maxs[i][2],2]) for i in 1:ndata ]
    
    # get signal intensities range
    δs   = [ maxs[i][1] - mins[i][1] for i=1:ndata ]

    # set data cutoffs
    cutoff_idx = zeros(Int32, ndata)
    for (k,kt) in enumerate(kinetic_traces)
        vmax, imax = maxs[k]
        iend = findfirst(kt[imax:end,2] .< cutoff*vmax )
        cutoff_idx[k] = typeof(iend) == Int ? iend + imax : size(kt,1) 
    end
    df2fit[!,:cutoff] = cutoff_idx

    return kinetic_traces, maxs, mins, δs

end