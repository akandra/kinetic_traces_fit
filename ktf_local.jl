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

# Set what is to be done
# what_to_do   = "fit" or "analysis"
what_to_do = ("fit",      "rrr", "facet", "tag")
what_to_do = ("analysis", "rrr", "facet", "tag")
T_cutoff = 480.0 # max temperature for Arrhenius fit

# load kinetic traces
kinetic_traces, maxs, mins, δs = ksr.load_kinetic_traces(df2fit,data_path, cutoff=-Inf)

# set initial values for the fitting parameters and other defaults
# units: μs⁻¹ for prefactors and rates; eV for energy

ksr.guess_rate!(df2fit, name= "k1", ν=1.0*10^5,  ϵ=0.2,  var=false)
ksr.guess_rate!(df2fit, name="km1", ν=1.0*10^5,  ϵ=2.0,  var=false)
ksr.guess_rate!(df2fit, name= "k2", ν=1.0*10^5,  ϵ=0.4,  var=false)
ksr.guess_rate!(df2fit, name= "k3", ν=3.0*10^4,  ϵ=0.36, var=true)
ksr.guess_rate!(df2fit, name= "k4", ν=1.0*10^6,  ϵ=0.75, var=false)
ksr.guess_rate!(df2fit, name= "k5", ν=1.0*10^10, ϵ=0.5,  var=false)

ksr.guess_par!(df2fit, name= "a",   value= 0.25*first.(maxs), min=0.001, var=true, glbl=false)
ksr.guess_par!(df2fit, name="t0",   value=-100.0,  min=-200.0, max=200.0, var=true, glbl=false)
ksr.guess_par!(df2fit, name="f_tr", value=1e-3, min=0.0, var=true, glbl=false)
ksr.guess_par!(df2fit, name="k_vac",value=1e-5, var=false, glbl=false)

ksr.guess_par!(df2fit, name="baseline",value=set_baseline(ndata), var=false, glbl=false)
    

NEXTTIME: gloriously go on!
   
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