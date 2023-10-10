function hTST(temperature::Float64, hνR::Vector{Float64}, hνTS::Vector{Float64}, energy::Float64)

    kB = 8.61733*10^-5
    hP = 4.13567*10^-9
    kBTinv = 1.0/(kB*temperature) 
    QRinv  = prod(1.0 .- exp.(-hνR *kBTinv))
    QTSinv = prod(1.0 .- exp.(-hνTS*kBTinv))

    return (kB*temperature/hP)*QRinv/QTSinv*exp(-energy*kBTinv)
end

hTST(300.0,[0.01,0.1, 0.12 ], [0.11,0.13],0.5)