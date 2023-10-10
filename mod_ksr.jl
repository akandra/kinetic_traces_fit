module Kinetics_of_Surface_Reactions

using DataFrames
using DelimitedFiles
using Parameters
using DifferentialEquations
using LsqFit
using Plots
using Statistics: mean

# debugging output tag
debug = 0

# important constants
kB = 8.61733*10^-5 # conversion from K to eV
hP = 4.13567*10^-9 # Planck constant in eV⋅μs

# structure to deal with fitting parameters
@with_kw mutable struct fitpar

    value::Float64 = 0.0
    min::Float64   = -Inf
    max::Float64   = Inf
    var::Bool      = false
    glbl::Bool     = false

end

# functions logging things from input
include("log_stuff.jl")

# in/out functions
include("get_pump_beam.jl")
include("get_data.jl")
include("load_kinetic_traces.jl")

# utilities
include("temperature_functions.jl")
include("ann_par.jl") # annotation for plots

# initialization functions
include("ini_guess.jl")

# functions getting results of previous fits
include("get_results_local.jl")
include("get_results_global.jl")

# functions doing fits
include("do_local_fit.jl")
include("do_global_fit.jl")

# function calculating the fitting function
include("product_flux.jl")

# function creating the dataframe
include("create_df.jl")

# do-it function
include("do_it.jl")

end