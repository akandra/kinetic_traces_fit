
using DelimitedFiles
using DataFrames

function get_data(filenames)

      #dropbox_path = "c:/users/akandra/Dropbox/"
      dropbox_path = "../../"
      data_path = dropbox_path * "Kinetics of Surface Reactions/data/"
  
      data = Vector{Matrix{Float64}}
      data= []
      for f in filenames
          push!(data, readdlm(data_path * f, Float64,comments=true))
      end
      data
  end
  
dropbox_path = "../../"
data_path = dropbox_path * "Kinetics of Surface Reactions/data/"

datafilenames = readdir(data_path)

beamfnames      = filter(x-> occursin("beam",x) ,datafilenames)
Oinifnames      = filter(x-> occursin("Oini",x) ,datafilenames)
kintracesfnames = filter(x->!(occursin("beam",x) || occursin("Oini",x)) ,datafilenames)

# create kinetic traces data frame
tags =[ map( f->split( splitext(f)[1], "-")[i], kintracesfnames ) for i in 1:5 ]
dfkt = DataFrame(tags, [:tag,:facet,:temperature,:rrH2,:rrO2])
dfkt[!,:temperature] = parse.(Float64,dfkt[!,:temperature])
dfkt[!,:rrH2] = parse.(Float64,dfkt[!,:rrH2])
dfkt[!,:rrO2] = parse.(Float64,dfkt[!,:rrO2])
dfkt[!,:ktfname] = kintracesfnames
dfkt[!,:beamfname] = tags[1] .* "-beam.dat"

# create [O]_ini data frame
Oinidata = get_data(Oinifnames)
tagsOini = map( f->split( splitext(f)[1], "-")[1], Oinifnames )
tags1 = vcat(fill.(tagsOini,size.(Oinidata,1))...)
dfOini = DataFrame(vcat(Oinidata...),[:temperature,:rrH2,:rrO2,:Oini])
dfOini[!,:tag] = tags1

# join above data frames
df = innerjoin(dfkt, dfOini, on = [:tag, :temperature, :rrH2, :rrO2])
