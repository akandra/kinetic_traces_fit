# Built-In Fitting Model

# Set the names for fit parameter dataframe columns
fit_parnames = ["a", "t0", "f_tr", "k_vac", "baseline", "cutoff" ]

# set initial values for baselines
function set_baseline(shift::Int, kinetic_trace::Matrix{Float64}, 
                                  mins::Tuple{Number,Number}, 
                                  maxs::Tuple{Number,Number}, 
                                  δs::Float64)

    bl_range = findfirst( x -> x>mins[1]+δs/2 , kinetic_trace[1:maxs[2],2]) - shift
    return mean(kinetic_trace[1:bl_range,2])

end 


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
function fit_function(sol::Vector{Float64} ,dfpar::DataFrame, times::Vector{Float64})

    # get values of fitting parameters from the data frame

    t0       = dfpar.t0[1].value
    a        = dfpar.a[1].value
    baseline = dfpar.baseline[1].value
    f_tr     = dfpar.f_tr[1].value
    k_vac    = dfpar.k_vac[1].value

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
