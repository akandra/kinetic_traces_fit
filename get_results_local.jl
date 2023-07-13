"""
Gets the parameters related to the best local fit

get_results_local(fit_filename::AbstractString)

Returns:

a tuple of initial guesses, best fit parameters, standard error and χ2 for the fit with the smallest value of χ2 
if local_fit_filename exists
otherwise nothing
"""
function get_results_local(filename; crit = "best")

    if isfile(filename)

        ini_guess = Float64[]
        best_fit  = Float64[]
        se        = Float64[]
        chi2      = Float64

        data = readdlm(filename,comments=true) 
        if crit == "best" 
            chi2ind = argmin(data[4:4:end,1])  # get the fit with minimal chi2
        elseif crit == "last" 
            chi2ind = lastindex(data[4:4:end,1])  # get the last fit
        else 
            println("get_results_local: crit is unknown, takeing the default value.")
        end
        ini_guess = data[1+4*(chi2ind-1),:]
        best_fit  = data[2+4*(chi2ind-1),:]
        se        = data[3+4*(chi2ind-1),:]
        chi2      = data[4+4*(chi2ind-1),1]

        return ini_guess, best_fit, se, chi2
    else
        return nothing
    end
end
