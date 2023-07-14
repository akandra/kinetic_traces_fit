# Theo's All Step Model

# Set the names for fit parameter dataframe columns
# Note: use "ν" to start a name of prefactor parameter,
#       and "ϵ" to start a name of energy parameter
#fitparsnames_model = ["ν1","ϵ1","νm1","ϵm1","ν2","ϵ2","ν3","ϵ3","ν4","ϵ4","ν5","ϵ5"]

species = Dict("H₂"=>1, "O"=>2, "OH"=>3, "H"=>4, "H₂O"=>5)

#rate_constants = Dict("k1"=>1, "km1"=>2, "k2"=>3, "k3"=>4, "k4"=>5, "k5"=>6)
rate_constants_base = "k"
rate_constants_sfx = ["1", "m1", "2", "3", "4", "5"]

# For the future: make it possible to use any name for the function
function kin_model!(ydot,y,p,t, beampars, θs)

    # All-step model from Theo (everything occurs on steps)
    # (1)  H₂ + O <- km1, k1 -> OH + H
    # (2)  H + OH -k2> H₂O 
    # (3) OH + OH -k3> H₂O + O
    # (4)  H + H  -k4> H₂ 
    # (2)  H₂O -k5> H₂O(gas)  

    iH   = species["H"]
    iH2  = species["H₂"]
    iH2O = species["H₂O"]
    iO   = species["O"]
    iOH  = species["OH"]

    σO   = occ_factors["O"]
    σOH  = occ_factors["OH"]
    σH2O = occ_factors["H₂O"]
    σH   = occ_factors["H"]

    κ1, κm1, κ2, κ3, κ4, κ5 = p

    a, fwhm, tcenter        = beampars
    

    ydot[iH2] = -κ1*y[iH2]*y[iO]*θs/σO + θs*κ4*y[iOH]*y[iOH]/(σH*σH) + θs*κm1*y[iOH]*y[iH]/(σH*σOH) + θs*pump_pulse(t,a,fwhm,tcenter)
    ydot[iO]  = -κ1*y[iH2]*y[iO] + κm1*y[iOH]*y[iH]*σO/(σH*σOH) +  κ3*y[iH]*y[iH]*σO/(σOH*σOH)
    ydot[iOH] =  κ1*y[iH2]*y[iO]/σO - 2.0*κ4*y[iOH]*y[iOH]/σH - κ2*y[iOH]*y[iH]/σOH - κm1*y[iOH]*y[iH]/σOH
    ydot[iH]  =  κ1*y[iH2]*y[iO]/σO - κ2*y[iOH]*y[iH]/σH - κm1*y[iOH]*y[iH]/σH - 2.0*κ3*y[iH]*y[iH]/σOH
    ydot[iH2O]=  κ2*y[iOH]*y[iH]*σH2O/(σH*σOH) + κ3*y[iH]*y[iH]*σH2O/(σOH*σOH) - κ5*y[iH2O]   

end
