#function Arrhenius(temperature::Float64, prefactor::Float64, energy::Float64)::Float64
#        return prefactor*exp(-11604.5*energy/temperature)
#end
function Arrhenius(temperature, prefactor, energy)
    return prefactor*exp(-11604.5*energy/temperature)
end