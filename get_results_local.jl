"""
Gets the parameters related to the best local fit

get_results_local(fit_filename::AbstractString)

Returns:

a tuple of initial guesses, best fit parameters, standard error and χ2 for the fit with the smallest value of χ2 
if local_fit_filename exists
otherwise nothing
"""
function get_results_local(path, filename; crit = "best")

    fname = joinpath(path,filename)
    if isfile(fname)

        # ini_guess = Float64[]
        # best_fit  = Float64[]
        # se        = Float64[]
        # chi2      = Float64


        data = readdlm(fname,comments=true) 
        if crit == "best" 
            chi2ind = argmin(data[5:5:end,1])  # get the fit with minimal chi2
        elseif crit == "last" 
            chi2ind = lastindex(data[5:5:end,1])  # get the last fit
        else 
            error("get_results_local: crit key is unknown.")
        end
        
        par_names = data[1+5*(chi2ind-1),:]
        ini_guess = Dict(par_names .=> data[2+5*(chi2ind-1),:])
        best_fit  = Dict(par_names .=> data[3+5*(chi2ind-1),:])
        se        = Dict(par_names .=> data[4+5*(chi2ind-1),:])
        chi2      = data[5+5*(chi2ind-1),1]

        return ini_guess, best_fit, se, chi2
    else
        return nothing
    end
end
