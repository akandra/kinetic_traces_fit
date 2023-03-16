using Plots
using LsqFit
using Optim
using BenchmarkTools

 @. expDecay(x, p, a) = p[1]*exp(-p[2]*x) + a

x = [0:0.05:1;]
pfix = 0.5
data = expDecay(x,[1,3],pfix) + 0.05*randn(length(x))

a0 = 0.1
k0 = 1.0

println("Fitting...")
@time fit = curve_fit( (x, p) -> expDecay(x,p,pfix), x, data, [a0,k0])

println("Plotting...")
display(plot(x,[ data, expDecay(x,coef(fit),pfix) ],seriestype=[:scatter :line :line]))
#display(plot(kinetic_traces[:,1],
#    [kinetic_traces[:,2], H2OProduction(kinetic_traces[:,1],coef(fit)) ],
#    seriestype = [:scatter :line], framestyle=:box,label=["data" "fit"]))

