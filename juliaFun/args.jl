function f(args...)
push!(conds,args)
return args
end 

conds= []

f(1,2,"a","s")
f("2nd")
println([ size(c,1) for c in conds])