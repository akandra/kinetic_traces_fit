function get_results_global(path, filenames)

    filename = "global_fit.dat"

    filenames = sort(filenames)
    fnames    = String[]
    ini_guess = Float64[]
    best_fit  = Float64[]
    se        = Float64[]
    chi2      = Inf64
    ind       = 0 

    
    data      = readdlm(path*filename)
    ifn = first.( Tuple.(findall(s-> s=="file",data)) ) .+ 1     # indices for file names
    iig = first.( Tuple.(findall(s-> s=="initial",data)) ) .+ 1  # indices for initial guesses
    ibf = first.( Tuple.(findall(s-> s=="best",data)) ) .+ 1     # indices for best fit parameters
    ise = first.( Tuple.(findall(s-> s=="standard",data)) ) .+ 1 # indices for standard errors
    iχ2 = first.( Tuple.(findall(s-> s=="reduced",data)) ) .+ 1  # indices for reduced χ2
    nds = iig .- ifn .- 1
    # find a dataset that matches filenames and has the smallest χ2
    for i in 1:length(ifn)
        if filenames == sort(data[ifn[i]:ifn[i]+nds[i]-1,1])
            if data[iχ2[i]] < chi2 
                chi2 = data[iχ2[i]]
                ind = i
            end
        end
    end

    ini_guess = data[iig[ind]:iig[ind]+nds[ind]-1,:]
    best_fit  = data[ibf[ind]:ibf[ind]+nds[ind]-1,:]
    se        = data[ise[ind]:ise[ind]+nds[ind]-1,:]
    chi2      = data[iχ2[ind]]

    return ini_guess, best_fit, se, chi2
end
