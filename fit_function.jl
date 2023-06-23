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

function fit_function()

end
