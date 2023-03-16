include("structures.jl")
include("load_data.jl")

using .Structs
using .Load


# set paths to the data and results folders
path         = "../../Dropbox/Kinetics of Surface Reactions/"
data_path    = path * "data/"
results_path = path * "results/" * "nu3_1/"

# create dataframe
df = create_df(data_path)






pars = Structs.fitpar()

load_kinetic_traces(df2fit)