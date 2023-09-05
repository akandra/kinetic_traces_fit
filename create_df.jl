function create_df(data_path::String, delim::String, pump_sfx::String, cov_sfx::String)

    datafilenames = readdir(data_path)
    pumpfnames      = filter(x-> occursin(pump_sfx,x) ,datafilenames)
    covfnames       = filter(x-> occursin(cov_sfx,x) ,datafilenames)
    kintracesfnames = filter(x->!(occursin(pump_sfx,x) || occursin(cov_sfx,x)) ,datafilenames)

    # create kinetic traces data frame
    # Warning:  after changing the storing of the tag information
    #           do not forget to make the corresponding adjustments in the following code
    tags =[ map( f->split( splitext(f)[1], delim)[i], kintracesfnames ) for i in 1:5 ]
    df = DataFrame(tags, [:tag,:facet,:temperature,:rr_pump,:rr_cov])
    df[!,:temperature] = parse.(Float64,df[!,:temperature])
    df[!,:rr_pump] = parse.(Float64,df[!,:rr_pump])
    df[!,:rr_cov] = parse.(Float64,df[!,:rr_cov])
    df[!,:rrr] = df.rr_cov ./ df.rr_pump

    df[!,:ktfname] = kintracesfnames
    df[!,:pump_file] = tags[1] .* (delim * pump_sfx * ".dat")

    # create pump-pulse data frame
    dfpump = get_pump_beam(pumpfnames)
    df = innerjoin(df,dfpump, on=:pump_file)

    # create initial coverage data frame
    cov_name, covdata = get_data_cov(data_path, covfnames)
    tagscov = map( f->split( splitext(f)[1], delim)[1], covfnames )
    tags1 = vcat(fill.(tagscov,size.(covdata,1))...)
    dfcov = DataFrame(vcat(covdata...),[:temperature,:rr_pump,:rr_cov,:cov0])
    dfcov[!,:cov_species] .= cov_name
    dfcov[!,:tag] = tags1
    # join above data frames
    df = innerjoin(df, dfcov, on = [:tag, :temperature, :rr_pump, :rr_cov])

    # create step density data frame
    dfθs = DataFrame( facet =        [ k[:facet] for k in θs ],
                      step_density = [ k[:value] for k in θs ])
    # join to the total data frame
    df = innerjoin(df,dfθs, on=:facet)

    return df

end