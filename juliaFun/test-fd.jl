using FiniteDifferences
using FiniteDiff

f(x) = exp(x[1]x[2]^2)

function Avv!(dir_deriv,p,george)
    for i=1:length(george)
        dir_deriv[i] = transpose(george) * 
        FiniteDiff.finite_difference_hessian( f, p[i]) * george
    end
end

println("grad from FiniteDifferences is ", 
    grad(central_fdm(5,1), f, [1.0,1.0])[1]) 

println("jacobian from FiniteDifferences is ",
    jacobian(central_fdm(5,2), f, [1.0,1.0])[1])

println("jacobian from FiniteDiff is ",
    FiniteDiff.finite_difference_jacobian( f,  [1.0,1.0]))

p = [ [1.0,1.0],[0.0,1.0] ]   
println("hessian in ", p[1]," is ",
    FiniteDiff.finite_difference_hessian( f, p[1]) )
println("hessian in ", p[2]," is ",
    FiniteDiff.finite_difference_hessian( f, p[2]) )

output = zeros(2)
Avv!(output, p, [0,1])
println("directional 2nd derivatives are ",output)

