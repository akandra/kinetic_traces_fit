# Theo's All Step Model extended with O terrace diffusion

# species list
species = [
    adsorbate("H₂", "gas", "beam"),
    adsorbate(name="O",  site="step",    σ = 2, type = "coverage"),
    adsorbate(name="O",  site="terrace", σ = 4, type = "coverage"),
    adsorbate(name="OH", site="step",    σ = 1, type = "intermediate"),
    adsorbate(name="H",  site="step",    σ = 1, type = "intermediate"),
    adsorbate(name="H₂O",site="step",    σ = 1, type = "product")
]


# species indices
iH2 ::Int8 = 1
iOs ::Int8 = 2
iOt ::Int8 = 3
iOH ::Int8 = 4
iH  ::Int8 = 5
iH2O::Int8 = 6

# occupancies list
σO  :: Int8 = 2
σOH :: Int8 = 1
σH2O:: Int8 = 1
σH  :: Int8 = 1

rate_constants_base = "k"
rate_constants_sfx = ["1", "m1", "2", "3", "4", "5", "m5", "6"]

# For the future: make it possible to use any name for the function
function kin_model!(ydot::Vector{Float64},y::Vector{Float64},p::Vector{Float64},t::Float64, 
                    beampars::Vector{Vector{Float64}}, θs::Float64, Tsurf::Float64)
    # All-step O-diffusion modified model from Theo (everything occurs on steps)
    # (1)  H₂ + Os <- km1, k1 -> OH + H
    # (2)  H + OH -k2> H₂O 
    # (3) OH + OH -k3> H₂O + O
    # (4)  H + H  -k4> H₂ 
    # (5)  Os <- km5, k5 -> Ot  
    # (6)  H₂O -k6> H₂O(gas)  

    κ1, κm1, κ2, κ3, κ4, κ5, κm5, κ6 = p

    a, fwhm, tcenter = beampars
    
    beam = pump_pulse(t,a,fwhm,tcenter)

    ydot[iH2] = -κ1*y[iH2]*y[iO]*θs/σO + θs*κ4*y[iOH]*y[iOH]/(σH*σH) + θs*κm1*y[iOH]*y[iH]/(σH*σOH) + θs*beam
    ydot[iOs] = -κ1*y[iH2]*y[iO] + κm1*y[iOH]*y[iH]*σO/(σH*σOH) +  κ3*y[iH]*y[iH]*σO/(σOH*σOH)
                - (κ5*y[iOs] - κm5*y[iOt])
    ydot[iOt] =  κ5*y[iOs] - κm5*y[iOt]
    ydot[iOH] =  κ1*y[iH2]*y[iO]/σO - 2.0*κ4*y[iOH]*y[iOH]/σH - κ2*y[iOH]*y[iH]/σOH - κm1*y[iOH]*y[iH]/σOH
    ydot[iH]  =  κ1*y[iH2]*y[iO]/σO - κ2*y[iOH]*y[iH]/σH - κm1*y[iOH]*y[iH]/σH - 2.0*κ3*y[iH]*y[iH]/σOH
    ydot[iH2O]=  κ2*y[iOH]*y[iH]*σH2O/(σH*σOH) + κ3*y[iH]*y[iH]*σH2O/(σOH*σOH) - κ6*y[iH2O]   
end 

function cov_species_eq!(Oini::Float64, rates::Vector{Float64}, θs::Float64, Tsurf::Float64)

    # Terrace-step O-diffusion
    # Os <- km5, k5 -> Ot  
    
    κ5  = rates[6]
    κm5 = rates[7]

    cOs = Oini/(κ5/κm5 + 1)
    cOt = Oini*(1 - 1/(κ5/κm5 + 1))

    return cOs, cOt
end 

function prod_flux!(sol, θs::Float64, rates::Vector{Float64})

    # For future reference γ is available
    # γ = θs/(2 - θs)

    κ6 = rates[8]
    return κ6*sol[iH2O,:]
    
end 
