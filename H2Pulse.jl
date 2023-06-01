function H2Pulse(t, a::Vector{Float64}, fwhm::Vector{Float64}, tcenter::Vector{Float64})

    sigma = fwhm / ( 2.0*sqrt(2.0*log(2.0)) )
    
    gausses = 0
    for i=1:3
        gausses = gausses + a[i]/(sqrt(2*π)*sigma[i])*exp( -(t - tcenter[i])^2/(2*sigma[i]^2) )
    end

    return gausses
end
