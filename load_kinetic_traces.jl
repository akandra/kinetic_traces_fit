function load_kinetic_traces(df2fit, data_path)

    println("Loading Kinetic Traces...")
    [ println(" ",filename) for filename in df2fit.ktfname ]
    println("")

    kinetic_traces = get_data(data_path, df2fit.ktfname)
    ndata = length(kinetic_traces)

    # finding maximum intensity of a signal
    maxs = [ findmax(kt[:,2]) for kt in kinetic_traces ]
    
    # finding minimum intensity at the l.h.s. from a maximum
    mins = [ findmin(kinetic_traces[i][1:maxs[i][2],2]) for i in 1:ndata ]
    
    # get signal intensities range
    δs   = [ maxs[i][1] - mins[i][1] for i=1:ndata ]

    return kinetic_traces, maxs, mins, δs

end