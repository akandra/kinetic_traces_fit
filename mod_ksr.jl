module Kinetics_of_Surface_Reactions

using DataFrames
using DelimitedFiles
using Parameters
using DifferentialEquations
using LsqFit
using Plots

# structure to deal with fitting parameters
@with_kw mutable struct fitpar

    value::Float64 = 0.0
    min::Float64   = -Inf
    max::Float64   = Inf
    var::Bool      = false
    glbl::Bool     = false

end

# in/out functions
include("get_pump_beam.jl")
include("get_data.jl")
include("load_kinetic_traces.jl")

# utilities
include("Arrhenius.jl")
include("ann_par.jl") # annotation for plots

# initialization functions
include("ini_guess.jl")

# functions getting results of previous fits
include("get_results_local.jl")
include("get_results_global.jl")

# functions doing fits
include("global_fit.jl")
include("local_fit.jl")

# function calculating the fitting function
include("H2OProduction.jl")

# function defining the kinetic model
include("kinetic_model.jl")
include("H2Pulse.jl")

# function creating the dataframe
include("create_df.jl")

end