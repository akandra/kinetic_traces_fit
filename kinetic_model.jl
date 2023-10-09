# Theo's All Step Model

# species list
species = ["H₂", "O", "OH", "H", "H₂O"]
# species indices
iH2 ::Int8 = 1
iO  ::Int8 = 2
iOH ::Int8 = 3
iH  ::Int8 = 4
iH2O::Int8 = 5
# occupancies list
σO  :: Int8 = 2
σOH :: Int8 = 1
σH2O:: Int8 = 1
σH  :: Int8 = 1

rate_constants_base = "k"
rate_constants_sfx = ["1", "m1", "2", "3", "4", "5"]

# For the future: make it possible to use any name for the function
function kin_model!(ydot::Vector{Float64},y::Vector{Float64},p::Vector{Float64},t::Float64, 
                    beampars::Vector{Vector{Float64}}, θs::Float64, Tsurf::Float64)
    # All-step model from Theo (everything occurs on steps)
    # (1)  H₂ + O <- km1, k1 -> OH + H
    # (2)  H + OH -k2> H₂O 
    # (3) OH + OH -k3> H₂O + O
    # (4)  H + H  -k4> H₂ 
    # (2)  H₂O -k5> H₂O(gas)  

    κ1, κm1, κ2, κ3, κ4, κ5 = p

    a, fwhm, tcenter = beampars
    
    beam = pump_pulse(t,a,fwhm,tcenter)

    ydot[iH2] = -κ1*y[iH2]*y[iO]*θs/σO + θs*κ4*y[iOH]*y[iOH]/(σH*σH) + θs*κm1*y[iOH]*y[iH]/(σH*σOH) + θs*beam
    ydot[iO]  = -κ1*y[iH2]*y[iO] + κm1*y[iOH]*y[iH]*σO/(σH*σOH) +  κ3*y[iH]*y[iH]*σO/(σOH*σOH)
    ydot[iOH] =  κ1*y[iH2]*y[iO]/σO - 2.0*κ4*y[iOH]*y[iOH]/σH - κ2*y[iOH]*y[iH]/σOH - κm1*y[iOH]*y[iH]/σOH
    ydot[iH]  =  κ1*y[iH2]*y[iO]/σO - κ2*y[iOH]*y[iH]/σH - κm1*y[iOH]*y[iH]/σH - 2.0*κ3*y[iH]*y[iH]/σOH
    ydot[iH2O]=  κ2*y[iOH]*y[iH]*σH2O/(σH*σOH) + κ3*y[iH]*y[iH]*σH2O/(σOH*σOH) - κ5*y[iH2O]   
end 

function prod_flux!(sol, θs::Float64, rates::Vector{Float64})

    # For future reference γ is available
    # γ = θs/(2 - θs)

    κ5 = rates[6]
    return κ5*sol[iH2O,:]
    
end 
