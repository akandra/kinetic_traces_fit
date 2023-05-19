# Built-In Fitting Model

# Set the names for fit parameter dataframe columns
fitparsnames_fit = ["a", "t0", "baseline", "f_tr", "k_vac" ]
# we avoid using fill(), see rebind_vs_mutate.jl in juliaFun to find out why
#[ df[!,n] = [fitpar() for _ in 1:nrow(df)] for n in fitparsnames]


function fit_function()

end
