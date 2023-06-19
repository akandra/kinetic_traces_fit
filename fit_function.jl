# Built-In Fitting Model

# Set the names for fit parameter dataframe columns
fit_parnames = ["a", "t0", "f_tr", "k_vac", "baseline", "cutoff" ]

# set initial values for baselines
function set_baseline(ndata::Int)

    baseline = zeros(ndata)
    for i=1:ndata
        bl_range = findfirst( x -> x>mins[i][1]+δs[i]/2 ,kinetic_traces[i][1:maxs[i][2],2]) - 7
        baseline[i] = mean(kinetic_traces[i][1:bl_range,2])
    end
    return baseline
end 

function fit_function()

end
