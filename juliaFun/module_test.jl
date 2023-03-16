include("m1.jl")
#push!(JULIA_LOAD_PATH,joinpath(pwd(),"juliaFun"))
using .m

m.hello("regards to George.")