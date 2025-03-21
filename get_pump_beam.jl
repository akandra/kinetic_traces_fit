# Consider the following:
#  1. Make this file user-definable
#  2. Make beam function an argument of kinetic model

function get_pump_beam(fnames)

    return DataFrame(pump_file=fnames, pumppars =
      [
        [[0.000120897, 0.000170694, 0.000105154],   [30.3872, 15.1266, 75.293],  [51.9994, 52.4238, 85.3074]], 
        #[[0.0000148181, 0.000212279, 0.000128483],  [14.8624, 18.5067, 68.9036], [35.5675, 52.2248, 72.4033]], 
        #[[0.000256086, 0.0000951222, 0.000241312],  [22.4236, 12.2085, 65.8626], [47.6367, 49.1664, 71.6434]], 
        # the weird 4th beam is replaced by the regular 6th beam
        # [[0.0000402161, 0.000242335, 0.000130398],  [11.7742, 20.1696, 177.909], [82.8323, 97.182, 96.105]], 
        #[[0.0000754665, 0.0000707978, 0.000257094], [196.966, 16.4958, 25.2458], [156.164, 100.393, 111.913]],
        #[[0.0000645202, 0.000256638, 0.0000599968], [11.6545, 22.2035, 52.3407], [100.279, 112.345, 123.511]], 
        #[[0.0000754665, 0.0000707978, 0.000257094], [196.966, 16.4958, 25.2458], [156.164, 100.393, 111.913]]
            ] )
end

function pump_pulse(t::Float64, a::Vector{Float64}, fwhm::Vector{Float64}, tcenter::Vector{Float64})
    
    factor::Float64 = 1.0 / ( 2.0*sqrt(2.0*log(2.0)) )
    
    gausses::Float64 = 0.0
    for i::Int=1:3
        gausses = gausses + a[i]/(sqrt(2*π)*fwhm[i]*factor)*exp( -(t - tcenter[i])^2/(2*(fwhm[i]*factor)^2) )
    end

    return gausses
end
