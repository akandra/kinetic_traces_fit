using Plots
using LsqFit
using Parameters
using BenchmarkTools

function expDecays(x, p, pars, data)

    # update the values of variable parameters in pars
    j=0
    for (k,d) in fitpars
        if d.glbl 
            if d.var[1] 
                j = j + 1
                d.value[1:length(data)] .= p[j]
            end
        else
            for i=1:length(data)
                if d.var[i] 
                    j = j + 1
                    d.value[i] = p[j]
                end
            end
        end
    end
    
    # calculate a model function
    collect(Iterators.flatten(
        [   pars[:a].value[i].*exp.(-pars[:k].value[i].*datum[:,1]) .+ pars[:b].value[i]
        for (i, datum) in enumerate(data) ]
    ))
end

ndata = 2

@with_kw mutable struct FitPar
    value:: Vector{Float64} = zeros(ndata)
    min::   Vector{Float64} = fill(-Inf,ndata)
    max::   Vector{Float64} = fill( Inf,ndata)
    var::   Vector{Bool}    = fill(false,ndata)
    glbl::  Bool            = false 
end

fitpars = Dict{Symbol,FitPar}()
fitpars[:a] = FitPar(value=[0.1,0.2],var=[true, true])
fitpars[:k] = FitPar(value=[1.0,2.0],var=[true, true],glbl=true)
fitpars[:b] = FitPar(value=[0.1,0.2],var=[false, false])

println(length(filter(d->last(d).glbl && any(y->y!=last(d).value[1],last(d).value) ,fitpars)))
length(filter(d->last(d).glbl && any(y->y!=last(d).value[1],last(d).value) ,fitpars))>0 ? println("global parameter values fixed") : println("ok")

# construct the data mimicking the kinetic traces
data1 =  [ [x, 1.0*exp(-3*x) + 0.1 + 0.05*randn() ] for x in 0:0.01:1.0 ]
data2 =  [ [x, 1.5*exp(-4*x) + 0.2 + 0.05*randn() ] for x in 0:0.02:1.5 ]
data = Vector{Matrix{Float64}}
data = [ transpose(reduce(hcat, data1)), transpose(reduce(hcat, data2)) ]

xdata = collect(Iterators.flatten([ x[:,1] for x in data ]))
ydata = collect(Iterators.flatten([ x[:,2] for x in data ]))

pini=Float64[]
for (k,d) in fitpars
    n = d.glbl ? 1 : length(data)
    for i=1:n
        if d.var[i] 
            push!(pini, d.value[i])
        end
    end
end

println("Fitting...")
fit = curve_fit( (x,p)->expDecays(xdata, p, fitpars, data), xdata, ydata, pini)

yfit = expDecays(xdata, fit.param, fitpars, data)

println("Plotting...")
plots = []
icounter = 1
for (i,d) in enumerate(data)
    push!(plots, 
          plot(d[:,1], [d[:,2], yfit[icounter:icounter+size(d,1)-1]],
            seriestype=[:scatter :line], framestyle=:box, label=["data" "fit"],
            xlabel="ξ", ylabel="Ψ", title="data set "*string(i))
    )
    global icounter = icounter + size(d,1)
end
display( plot(plots...) )

