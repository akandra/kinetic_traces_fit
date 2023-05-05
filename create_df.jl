function create_df(data_path)

    datafilenames = readdir(data_path)
    beamfnames      = filter(x-> occursin("beam",x) ,datafilenames)
    Oinifnames      = filter(x-> occursin("Oini",x) ,datafilenames)
    kintracesfnames = filter(x->!(occursin("beam",x) || occursin("Oini",x)) ,datafilenames)[8:12] # Warning: don't forget to remove selection

    # create kinetic traces data frame
    tags =[ map( f->split( splitext(f)[1], "-")[i], kintracesfnames ) for i in 1:5 ]
    df = DataFrame(tags[1:3], [:tag,:facet,:temperature])
    df[!,:temperature] = parse.(Float64,df[!,:temperature])
    df[!,:reprates] = [ parse.(Float64,[tags[4][i],tags[5][i]]) for i in 1:length(tags[4]) ]

    println(df)
    stop

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

end