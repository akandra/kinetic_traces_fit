# Theo's All Step Model

# Set the names for fit parameter dataframe columns
# Note: use "ν" to start a name of prefactor parameter,
#       and "ϵ" to start a name of energy parameter
#fitparsnames_model = ["ν1","ϵ1","νm1","ϵm1","ν2","ϵ2","ν3","ϵ3","ν4","ϵ4","ν5","ϵ5"]

species_model = ["H₂", "O", "OH", "H", "H₂O"]

rate_constants_model = ["k1","km1","k2","k3","k4","k5"]


function eqns!(ydot,y,p,t, beampars, geompars)

    # All-step model from Theo (everything occurs on steps)
    # (1)  H₂ + O <- km1, k1 -> OH + H
    # (2)  H + OH -k2> H₂O 
    # (3) OH + OH -k3> H₂O + O
    # (4)  H + H  -k4> H₂ 
    # (2)  H₂O -k5> H₂O(gas)  

    θ_s, σO, σOH, σH2O, σH  = geompars
    κ1, κm1, κ2, κ3, κ4, κ5 = p
    a, fwhm, tcenter        = beampars
    
    ydot[1] = -κ1*y[1]*y[2]*θ_s/σO + θ_s*κ4*y[3]*y[3]/(σH*σH) + θ_s*κm1*y[3]*y[4]/(σH*σOH) + θ_s*H2Pulse(t,a,fwhm,tcenter)
    ydot[2] = -κ1*y[1]*y[2] + κm1*y[3]*y[4]*σO/(σH*σOH) +  κ3*y[4]*y[4]*σO/(σOH*σOH)
    ydot[3] =  κ1*y[1]*y[2]/σO - 2.0*κ4*y[3]*y[3]/σH - κ2*y[3]*y[4]/σOH - κm1*y[3]*y[4]/σOH
    ydot[4] =  κ1*y[1]*y[2]/σO - κ2*y[3]*y[4]/σH - κm1*y[3]*y[4]/σH - 2.0*κ3*y[4]*y[4]/σOH
    ydot[5] =  κ2*y[3]*y[4]*σH2O/(σH*σOH) + κ3*y[4]*y[4]*σH2O/(σOH*σOH) - κ5*y[5]   

end
