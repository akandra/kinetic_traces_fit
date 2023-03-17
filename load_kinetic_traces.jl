function load_kinetic_traces(df2fit, data_path)

    println("Loading Kinetic Traces...")
    [ println(" ",filename) for filename in df2fit.ktfname ]
    kinetic_traces = get_data(data_path, df2fit.ktfname)
    ndata = length(kinetic_traces)
    maxs = [ findmax(kt[:,2]) for kt in kinetic_traces ]
    mins = [ findmin(kinetic_traces[i][1:maxs[i][2],2]) for i in 1:ndata ]
    δs   = [ maxs[i][1] - mins[i][1] for i=1:ndata ]

    return kinetic_traces, maxs, mins, δs

end