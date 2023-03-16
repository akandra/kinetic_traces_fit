
using DataFrames

dropbox_path = "../../"
data_path = dropbox_path * "Kinetics of Surface Reactions/data/"

datafilenames = readdir(data_path)

beamfnames      = filter(x-> occursin("beam",x) ,datafilenames)
Oinifnames      = filter(x-> occursin("Oini",x) ,datafilenames)
kintracesfnames = filter(x->!(occursin("beam",x) || occursin("Oini",x)) ,datafilenames)

tags =[ map( f->split( splitext(f)[1], "-")[i], kintracesfnames ) for i in 1:5 ]

df = DataFrame(tag=tags[1], facet=tags[2], temperature=tags[3], rr1=tags[4], rr2=tags[5], 
                  ktfname=kintracesfnames, 
                  beamfname=[ filter(x-> occursin(i,x), beamfnames)[1] for i in tags[1] ],
                   )
 
data4fit = [ 
            ["332", "373", "100", "200" ],
            ["332", "423", "100", "200" ],
            ["332", "473", "100", "200" ],
            ["111", "373",  "25",  "50" ],
            ["111", "423",  "25",  "50" ],
            ["111", "473",  "25",  "50" ]
            ]

[ filter([:facet,:temperature,:rr1,:rr2] =>
        (f,t,r1,r2) -> f == i[1] && t==i[2] && r1==i[3] && r2==i[4] , df).ktfname
  for i in data4fit ]
