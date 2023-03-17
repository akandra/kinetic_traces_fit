using DelimitedFiles
using DifferentialEquations
using LsqFit
using Parameters
using BenchmarkTools
using DataFrames
using Plots
# using Printf
# using LaTeXStrings
# using FiniteDiff
# ==========================================================================
# Main part
# ==========================================================================

# if any parameter is global then do a fit to all the data simultaneously
# dataset 1 is chosen, since all .glbl's are the same

if  any( [x[1].glbl for x in eachcol(df2fitpar)] )
    global_fit(df2fit,df2fitpar, kinetic_traces, results_path, iguess, 
              wtd = what_to_do )
else # fit each dataset separately
    local_fit(df2fit,df2fitpar, kinetic_traces, results_path, iguess, 
              wtd = what_to_do, max_T_fit = T_cutoff)
end

mean