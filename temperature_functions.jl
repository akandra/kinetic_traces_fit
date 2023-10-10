function Arrhenius(temperature, prefactor, energy)
    return prefactor*exp(-11604.5*energy/temperature)
end

function hTST(temperature::Float64, hνR::Vector{Float64}, hνTS::Vector{Float64}, energy::Float64)

    kBTinv = 1.0/(kB*temperature) 
    QRinv  = prod(1.0 .- exp.(-hνR *kBTinv))
    QTSinv = prod(1.0 .- exp.(-hνTS*kBTinv))

    return (kB*temperature/hP)*QRinv/QTSinv*exp(-energy*kBTinv)
end