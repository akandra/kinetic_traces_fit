# Built-In Fitting Model
"""
Returns a product flux calculated with a fitting function based on solution to the kinetic model

fit_function(sol,dfpar,times)

sol::Vector{Float64} is the solution of kinetic model

dfpar::DataFrame is the df with values of fitting parameters

times::Vector{Float64} is times from the kintetic trace data

A*{x(t-t0)/X + f_tr*exp[-k_vac*(t-t0)]*c(t-t0)/c(tmax)} + b,

where 

x(t) is the product flux from kinetic model

X = max{x} is the normalization factor

f_tr is the fraction of the trace contributing to the density build-up 

k_vac is the evacuation rate

c(t) = int_0^t x(t) is the amount of product produced to the time instant t

b is a baseline value

"""
# Set the names for fit parameter dataframe columns
fit_parnames = ["a", "t_0", "f_t", "k_vac", "baseline", "cutoff" ]

function fit_function(ODEsol, dfpar::DataFrame, times::Vector{Float64})

    # get values of fitting parameters from the data frame
    # NB: the names have to be consistent with fit_parnames

    t0       = dfpar[1,"t_0"].value
    a        = dfpar[1,"a"].value
    baseline = dfpar[1,"baseline"].value
    f_tr     = dfpar[1,"f_t"].value
    k_vac    = dfpar[1,"k_vac"].value

    # produce a solution vector for a product
    sol = ODEsol(times .- t0)[species[product_species],:]

    # define the model function

    normfactor = maximum(sol)
    if normfactor == 0
        return a*sol .+ baseline
    else
        sol_acc = zeros(length(sol))
        cumsum!(sol_acc, sol/normfactor)
        return a*(sol/normfactor + f_tr*exp.(-k_vac*(times .- t0)) .* sol_acc/sol_acc[end]) .+ baseline
    end

end

# set initial values for baselines
function set_baseline(shift::Int, kinetic_trace::Matrix{Float64}, 
    mins::Tuple{Number,Number}, 
    maxs::Tuple{Number,Number}, 
    δs::Float64)

bl_range = findfirst( x -> x>mins[1]+δs/2 , kinetic_trace[1:maxs[2],2]) - shift
return mean(kinetic_trace[1:bl_range,2])

end 

