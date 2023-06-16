using Statistics: mean

include("mod_ksr.jl")
import .Kinetics_of_Surface_Reactions as ksr

# Set paths
path         = "../../Dropbox/Kinetics of Surface Reactions/H-Oxidation-Pt/"
data_path    = path * "data/"
model_path   = "./"#path * "models/"

# Julia code containing a kinetic model
model_fn = "kinetic_model"#"Theos_All_Step_Model"
include(model_path * model_fn * ".jl")

# Julia code containing a fitting function
fit_function_fn = "./fit_function.jl"
include(fit_function_fn)

results_path = path * "results/" * model_fn

# tag delimiter in file names
delim    = "-"
# set suffices for beam file names
pump_sfx = "beam"
cov_sfx  = "Oini"

# create dataframe
df = ksr.create_df(data_path, delim, pump_sfx, cov_sfx)

# select kinetic trace data to fit
ktfnames = [
            "20170828-332-373-100-200.dat", 
             "20170828-332-423-100-200.dat", 
            "20170828-332-473-100-200.dat",
            "20201117-111-373-25-50.dat", 
            "20201117-111-423-25-50.dat", 
            "20201117-111-473-25-50.dat"
           ]

condition = df.ktfname .∈ [ktfnames[1:6]]
#df2fit = filter( :ktfname => in(ktfnames),df)
#condition = (df.rrr .== 2) # .| (df.rrr .== 8.0)
#condition = (df.tag .!= "20201103") .| (df.tag .!= "20201118") 
#df2fit = filter([:facet,:rrH2,:rrO2] => (f,h,o) ->  f=="332" && h==100 && o==100 ,df)
# dataframe with data to fit
#condition = 55:55

df2fit = df[condition,:]
ndata = ksr.nrow(df2fit)

# load kinetic traces
kinetic_traces, maxs, mins, δs = ksr.load_kinetic_traces(df2fit,data_path)

# Set what is to be done
# what_to_do   = "fit" or "analysis"

what_to_do = ("fit",      "rrr", "facet", "tag")
what_to_do = ("analysis", "rrr", "facet", "tag")
T_cutoff = 480.0 # max temperature for Arrhenius fit

# Construct a list of fitting parameter names
#fitparsnames = [ fitparsnames_model; fitparsnames_fit]

# set initial values for the fitting parameters and other defaults
# units: μs⁻¹ for prefactors and rates; eV for energy

# ksr.guess_Arrh!(df2fit, "k1",  1.0*10^5, 0.2, true)
# ksr.guess!(df2fit, "ν1",  ν_theoguess[1], true, false, "ϵ1", ϵ_theoguess[1], true, false)

#                       name         ν      ϵ     var
ksr.guess_rate!(df2fit, "k1",  1.0*10^5,  0.2,  false)
ksr.guess_rate!(df2fit,"km1",  1.0*10^5,  2.0,  false)
ksr.guess_rate!(df2fit, "k2",  1.0*10^5,  0.4,  false)
ksr.guess_rate!(df2fit, "k3",  3.0*10^4,  0.36, true)
ksr.guess_rate!(df2fit, "k4",  1.0*10^6,  0.75, false)
ksr.guess_rate!(df2fit, "k5",  1.0*10^10, 0.5,  false)

NEXTIME: make guess with optional keyword arguments

ksr.guess_par(df2fit, "a", val= maxs[:][1]*0.25, min=0.001, var=true)
ksr.guess_par(df2fit, "t0", val= -100.0, min=-200.0, max=200.0, var=true)

ksr.guess_a(df2fit, 0.25, 0.001, true)

# create df columns for fit function pars
[ df2fit[!,n] = [ksr.fitpar() for _ in 1:ndata] for n in fit_parnames]

for i=1:ndata

    df2fit.a[i].value = maxs[i][1]*0.25
    df2fit.a[i].min   = 0.001
    df2fit.a[i].var   = true

    df2fit.t0[i].value = -100.0
    df2fit.t0[i].min   = -200.0
    df2fit.t0[i].max   =  200.0
    df2fit.t0[i].var   = true

    df2fit.f_tr[i].value = 0.001
    df2fit.f_tr[i].min   = 0.0
    df2fit.f_tr[i].max   = Inf
    df2fit.f_tr[i].glbl  = false
    df2fit.f_tr[i].var   = true

    df2fit.k_vac[i].value = 1e-5
    df2fit.k_vac[i].glbl  = false
    df2fit.k_vac[i].var   = false
end

    
# set initial values for baselines
for i=1:ndata
    bl_range = findfirst( x -> x>mins[i][1]+δs[i]/2 ,kinetic_traces[i][1:maxs[i][2],2]) - 7
    df2fit.baseline[i].value = mean(kinetic_traces[i][1:bl_range,2])
    #df2fit.baseline[i].value = 0.0
end
    
ksr.set_cutoff(-10.0)
# set data cutoffs
cutoff = zeros(Int32, ndata)
for (k,kt) in enumerate(kinetic_traces)
    vmax, imax = maxs[k]
    iend = findfirst(kt[imax:end,2] .< -10.0*vmax )
    cutoff[k] = typeof(iend) == Int ? iend + imax : size(kt,1) 
end
df2fit[!,:cutoff] = cutoff
    
    # select df columns of type fitpar
    df2fitpar = df2fit[!,names(df2fit,ksr.fitpar)]
    # save initial guesses
    iguess = zeros(ndata, ksr.ncol(df2fitpar))
    for i in 1:ndata
        for j in 1:ksr.ncol(df2fitpar)
            iguess[i,j] = df2fitpar[i,j].value
        end
    end

# create geometry parameters data frame
dfgeom = DataFrame( facet = String[], geompars = Vector{Float64}[] )
push!(dfgeom, ( "332", [1.0/6.0, 2, 1, 1, 1]) )
push!(dfgeom, ( "111", [  0.005, 2, 1, 1, 1]) )
#df = innerjoin(df,dfgeom, on=:facet)

# if any parameter is global then do a fit to all the data simultaneously
# dataset 1 is chosen, since all .glbl's are the same

if  any( [x[1].glbl for x in eachcol(df2fitpar)] )
    ksr.global_fit(df2fit,df2fitpar, kinetic_traces, results_path, iguess, 
                   wtd = what_to_do )
else # fit each dataset separately
    ksr.local_fit(df2fit,df2fitpar, kinetic_traces, results_path, iguess, 
                  wtd = what_to_do, max_T_fit = T_cutoff)
end

#NEXTTIME: adjust the kt files and the df to our new understanding the things    
# Consider for the future development
#
# 1. Put all references to the df into the ksr-module
# 2. Change the way how the kinetic traces files are named and supplied with ancillary info.
#    Consider 3 ways of doing that:
#        - create a single ancillary file
#        - create an ancillary file for each kt file
#        - put ancillary information as a header into a kt file