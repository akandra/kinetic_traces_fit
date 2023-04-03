module Kinetics_of_Surface_Reactions

using DataFrames
using DelimitedFiles
using Parameters

# in/out functions
include("get_beampars.jl")
include("get_data.jl")
include("load_kinetic_traces.jl")

# initialization functions
include("Arrhenius.jl")
include("ini_guess.jl")

# structure to deal with fitting parameters
@with_kw mutable struct fitpar

    value::Float64 = 0.0
    min::Float64   = -Inf
    max::Float64   = Inf
    var::Bool      = false
    glbl::Bool     = false

end

include("create_df.jl")

end