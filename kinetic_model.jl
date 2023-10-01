# Theo's All Step Model

# NOTE FOR TOMORROW: how about DataFrame?
species = [
            adsorbate( "H₂",  "", 0),
            adsorbate(  "O", "s", 2),
            adsorbate( "OH", "s", 1),
            adsorbate(  "H", "s", 1),
            adsorbate("H₂O", "s", 1)
          ]

# consider to get the indices automatically from the above structure
iH2  :: Int = 1
iOs  :: Int = 2
iOHs :: Int = 3
iHs  :: Int = 4
iH2O :: Int = 5


rate_constants_base = "k"
rate_constants_sfx = ["1", "m1", "2", "3", "4", "5"]

# For the future: make it possible to use any name for the function
function kin_model!(ydot::Vector{Float64},y::Vector{Float64},p::Vector{Float64},t::Float64, beampars::Vector{Vector{Float64}}, θs)
    # All-step model from Theo (everything occurs on steps)
    # (1)  H₂ + O <- km1, k1 -> OH + H
    # (2)  H + OH -k2> H₂O 
    # (3) OH + OH -k3> H₂O + O
    # (4)  H + H  -k4> H₂ 
    # (2)  H₂O -k5> H₂O(gas)  

    σO  :: Int = occ_factors["O"]
    σOH :: Int = occ_factors["OH"]
    σH2O:: Int = occ_factors["H₂O"]
    σH  :: Int = occ_factors["H"]

    κ1, κm1, κ2, κ3, κ4, κ5 = p

    a, fwhm, tcenter        = beampars
    
    beam = pump_pulse(t,a,fwhm,tcenter)

    ydot[iH2] = -κ1*y[iH2]*y[iO]*θs/σO + θs*κ4*y[iOH]*y[iOH]/(σH*σH) + θs*κm1*y[iOH]*y[iH]/(σH*σOH) + θs*beam
    ydot[iO]  = -κ1*y[iH2]*y[iO] + κm1*y[iOH]*y[iH]*σO/(σH*σOH) +  κ3*y[iH]*y[iH]*σO/(σOH*σOH)
    ydot[iOH] =  κ1*y[iH2]*y[iO]/σO - 2.0*κ4*y[iOH]*y[iOH]/σH - κ2*y[iOH]*y[iH]/σOH - κm1*y[iOH]*y[iH]/σOH
    ydot[iH]  =  κ1*y[iH2]*y[iO]/σO - κ2*y[iOH]*y[iH]/σH - κm1*y[iOH]*y[iH]/σH - 2.0*κ3*y[iH]*y[iH]/σOH
    ydot[iH2O]=  κ2*y[iOH]*y[iH]*σH2O/(σH*σOH) + κ3*y[iH]*y[iH]*σH2O/(σOH*σOH) - κ5*y[iH2O]   
    
end 
