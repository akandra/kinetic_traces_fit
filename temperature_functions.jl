function Arrhenius(temperature, prefactor, energy)
    return prefactor*exp(-11604.5*energy/temperature)
end

function hTST(temperature::Float64, hνR::Vector{Float64}, hνTS::Vector{Float64}, energy::Float64)

    kBTinv = 1.0/(kB*temperature) 
    QRinv  = prod(1.0 .- exp.(-hνR *kBTinv))
    QTSinv = prod(1.0 .- exp.(-hνTS*kBTinv))

    return (kB*temperature/hP)*QRinv/QTSinv*exp(-energy*kBTinv)
end

# set the temperature function keys
function set_T_function_keys()

    global T_function_keys = zeros(Int,length(rate_constants),2)

    for (i,r) in enumerate(rate_constants_sfx)
        
        if r in [ v[:sfx]  for v in guess_Arrh_global ]
            T_function_keys[i,1] = 1
        end

        for v in guess_hTST_global 
            if r == v[:sfx]
                T_function_keys[i,1] = 2
                T_function_keys[i,2] = length(v[:hνR])
            end
        end
    
    end    

end
