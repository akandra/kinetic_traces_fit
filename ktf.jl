using Statistics: mean

include("mod_ksr.jl")
import .Kinetics_of_Surface_Reactions as ksr

# Set paths
path         = "../../Dropbox/Kinetics of Surface Reactions/H-Oxidation-Pt/"
data_path    = path * "data/"
model_path   = "./"#path * "models/"

# Julia code containing a kinetic model
model_fn = "eqns"#"Theos_All_Step_Model"
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
# fit_is_local = true or false
# what_to_do   = "fit" or "analysis"

fit_is_local = true
what_to_do = ("fit",      "rrr", "facet", "tag")
what_to_do = ("analysis", "rrr", "facet", "tag")
T_cutoff = 480.0 # max temperature for Arrhenius fit

# Construct a list of fitting parameter names
fitparsnames = [ fitparsnames_model; fitparsnames_fit]


# initial guesses for the Arrhenius parameters cooked by Theo
ν_theoguess = [ 1.0*10^5, 1.0*10^5, 1.0*10^5, 3.0*10^4, 1.0*10^6, 1.0*10^10 ] # μs
ϵ_theoguess = [ 0.2, 2.0, 0.4, 0.36, 0.75, 0.5 ] # eV

# set initial values for the fitting parameters and other defaults
# units: μs⁻¹ for prefactors and rates; eV for energy
#

ksr.ini_guess!(df2fit, "1", ν_theoguess[1], false, true, ϵ_theoguess[1], false, true)

if fit_is_local
    # -------------------------------------------------------------------------------
    # LOCAL FIT
    # -------------------------------------------------------------------------------
    
        # rates
        ksr.ini_guess!(df2fit, "1", ν_theoguess[1], false, false, ϵ_theoguess[1])
        ksr.ini_guess!(df2fit,"m1", ν_theoguess[2], false, false, ϵ_theoguess[2])
        ksr.ini_guess!(df2fit, "2", ν_theoguess[3], false, false, ϵ_theoguess[3])
        ksr.ini_guess!(df2fit, "3", ν_theoguess[4],  true, false, ϵ_theoguess[4])
        ksr.ini_guess!(df2fit, "4", ν_theoguess[5], false, false, ϵ_theoguess[5])
        ksr.ini_guess!(df2fit, "5", ν_theoguess[6], false, false, ϵ_theoguess[6])
    
        # amplitudes, t_0 etc.
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
        
    else
    # -------------------------------------------------------------------------------
    # GLOBAL FIT
    # -------------------------------------------------------------------------------
    
        # prefactors and activation energies
        ksr.ini_guess!(df2fit, "1", ν_theoguess[1], false, false, ϵ_theoguess[1], false, false)
        ksr.ini_guess!(df2fit,"m1", ν_theoguess[2], false, false, ϵ_theoguess[2], false, false)
        ksr.ini_guess!(df2fit, "2", ν_theoguess[3], false, false, ϵ_theoguess[3], false, false)
        ksr.ini_guess!(df2fit, "3", ν_theoguess[4],  true,  true, ϵ_theoguess[4],  true,  true)
        ksr.ini_guess!(df2fit, "3", 4.4204050640348025e7,  true,  true, 0.5327034287892817,  true,  true)
        ksr.ini_guess!(df2fit, "4", ν_theoguess[5], false, false, ϵ_theoguess[5], false, false)
        ksr.ini_guess!(df2fit, "5", ν_theoguess[6], false, false, ϵ_theoguess[6], false, false)
    
        # amplitudes, t_0s etc.
        for i=1:ndata
    
            # get initial guesses from local fit if it exists
            data = get_results_local(results_path, df2fit[i,:ktfname],crit="best")
            if !isnothing(data)
                # put the values into the dataframe
                for (j,p) in enumerate(names(df2fit,fitpar))
                    if p == "a"
                        df2fit[i,p].value = data[2][j]
                        df2fit[i,p].min   = 0.001
                        df2fit[i,p].var   = true
                    elseif p == "t0"
                        df2fit[i,p].value = data[2][j]
                        df2fit[i,p].var   = true
                    elseif (p == "f_tr") || (p == "k_vac")
                        df2fit[i,p].value = data[2][j]
                        df2fit[i,p].var   = false
                    end
                end
    
            # or set initial guesses manually
            else
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
                df2fit.f_tr[i].var   = true
            
                df2fit.k_vac[i].value = 1e-5
                df2fit.k_vac[i].var   = false
            end
    
        end
        
    end
    
    # set initial values for baselines
    for i=1:ndata
        bl_range = findfirst( x -> x>mins[i][1]+δs[i]/2 ,kinetic_traces[i][1:maxs[i][2],2]) - 7
        df2fit.baseline[i].value = mean(kinetic_traces[i][1:bl_range,2])
        #df2fit.baseline[i].value = 0.0
    end
    
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