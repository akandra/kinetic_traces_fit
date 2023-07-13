include("mod_ksr.jl")
import .Kinetics_of_Surface_Reactions as ksr

# Set paths
ksr.working_dir("../../Dropbox/Kinetics of Surface Reactions/H-Oxidation-Pt")
ksr.path_to_data("data")
ksr.path_to_model("",".")
ksr.path_to_fit_function("",".")

# Julia code containing a kinetic model and a fitting function
ksr.model("kinetic_model")
ksr.fit_function("fit_function")

ksr.path_to_output("out", @__FILE__)

# tag delimiter in file names
ksr.delim("-")
# set suffices for beam file names
ksr.pump_sfx("beam")
ksr.cov_sfx("Oini")

# set conditions selecting data
# ksr.conditions( :rrr => ==(2) )
ksr.conditions( "20170828-332-373-100-200.dat", 
                "20170828-332-423-100-200.dat", 
                "20170828-332-473-100-200.dat",
                "20201117-111-373-25-50.dat", 
                "20201117-111-423-25-50.dat", 
                "20201117-111-473-25-50.dat")
# ksr.conditions(2,3)
#ksr.conditions(10:20)

# Set what is to be done
# what_to_do   = "fit" or "analysis"
# !!!Warning: consider replacing what_to_do with tag_hierarchy
ksr.what_to_do("fit",      "rrr", "facet", "tag")
ksr.what_to_do("analysis", "rrr", "facet", "tag")

# Set maximum temperature for Arrhenius fits
ksr.Arrh_Tmax(480.0)

# Set kinetic trace data cutoff fraction
ksr.cutoff_fraction("off")

# set initial values for the fitting parameters and other defaults
# units: μs⁻¹ for prefactors and rates; eV for energy

ksr.guess_rate(name= "k1", ν=1.0*10^5,  ϵ=0.2,  var=false)
ksr.guess_rate(name="km1", ν=1.0*10^5,  ϵ=2.0,  var=false)
ksr.guess_rate(name= "k2", ν=1.0*10^5,  ϵ=0.4,  var=false)
ksr.guess_rate(name= "k3", ν=3.0*10^4,  ϵ=0.36, var=true)
ksr.guess_rate(name= "k4", ν=1.0*10^6,  ϵ=0.75, var=false)
ksr.guess_rate(name= "k5", ν=1.0*10^10, ϵ=0.5,  var=false)

ksr.guess_par(name= "a",   value= ("maxs", 0.25), min=0.001, var=true, glbl=false)

ksr.guess_par(name="t0",   value=-100,  min=-200.0, max=200.0, var=true, glbl=false)
ksr.guess_par(name="f_tr", value=1e-3, min=0.0, var=true, glbl=false)
ksr.guess_par(name="k_vac",value=1e-5, var=false, glbl=false)

ksr.guess_par(name="baseline",value=(ksr.set_baseline,7), var=false, glbl=false)

ksr.step_density(facet="332", value=1/6)
ksr.step_density(facet="111", value=0.005)

ksr.occupancies(species=  "O", value=2)
ksr.occupancies(species=  "H", value=1)
ksr.occupancies(species= "OH", value=1)
ksr.occupancies(species="H₂O", value=1)

ksr.products("H₂O")

ksr.do_it()