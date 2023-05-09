function create_df(data_path)

    # Warning: don't forget to remove selection
    datafilenames = readdir(data_path)
    pumpfnames      = filter(x-> occursin("beam",x) ,datafilenames)
    Oinifnames      = filter(x-> occursin("Oini",x) ,datafilenames)
    kintracesfnames = filter(x->!(occursin("beam",x) || occursin("Oini",x)) ,datafilenames)

    # create kinetic traces data frame
    # Warning:  after changing the storing of the tag information
    #           do not forget to make the corresponding adjustments in the following code
    tags =[ map( f->split( splitext(f)[1], "-")[i], kintracesfnames ) for i in 1:5 ]
    df = DataFrame(tags, [:tag,:facet,:temperature,:rr_pump,:rr_cov])
    df[!,:temperature] = parse.(Float64,df[!,:temperature])
    df[!,:rr_pump] = parse.(Float64,df[!,:rr_pump])
    df[!,:rr_cov] = parse.(Float64,df[!,:rr_cov])
    df[!,:rrr] = df.rr_cov ./ df.rr_pump

    df[!,:ktfname] = kintracesfnames
    df[!,:pump_file] = tags[1] .* "-beam.dat"

    # create pump-pulse data frame
    dfpump = DataFrame( pump_file = pumpfnames, pumppars = get_pump_beam(pumpfnames) )
    df = innerjoin(df,dfpump, on=:pump_file)

    # create [O]_ini data frame
    Oinidata = get_data(data_path, Oinifnames)
    tagsOini = map( f->split( splitext(f)[1], "-")[1], Oinifnames )
    tags1 = vcat(fill.(tagsOini,size.(Oinidata,1))...)
    dfOini = DataFrame(vcat(Oinidata...),[:temperature,:rrH2,:rrO2,:Oini])
    dfOini[!,:tag] = tags1
    # join above data frames
    df = innerjoin(df, dfOini, on = [:tag, :temperature, :rrH2, :rrO2])

    display(df)
    stop

    # create geometry parameters data frame
    dfgeom = DataFrame( facet = String[], geompars = Vector{Float64}[] )
    push!(dfgeom, ( "332", [1.0/6.0, 2, 1, 1, 1]) )
    push!(dfgeom, ( "111", [  0.005, 2, 1, 1, 1]) )
    df = innerjoin(df,dfgeom, on=:facet)

    # set the names for fit parameter columns
    fitparsnames = ["ν1","ϵ1","νm1","ϵm1","ν2","ϵ2","ν3","ϵ3","ν4","ϵ4","ν5","ϵ5",
                     "a", "t0", "baseline", "f_tr", "k_vac" ]
    # we avoid using fill(), see rebind_vs_mutate.jl in juliaFun to find out why
    [ df[!,n] = [fitpar() for _ in 1:nrow(df)] for n in fitparsnames]

    return df

end