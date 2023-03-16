using DifferentialEquations
using Plots

module m1
    using DelimitedFiles

    function get_data(filename)

        #dropbox_path = "c:/users/akandra/Dropbox/"
        dropbox_path = "../../"
        data_path = dropbox_path * "Kinetics of Surface Reactions/data/"
        data_fn = data_path * filename
        readdlm(data_fn, Float64)
    end
end

function parameterized_lorenz!(du,u,p,t)
    du[1] = p[1]*(u[2]-u[1])
    du[2] = u[1]*(p[2]-u[3]) - u[2]
    du[3] = u[1]*u[2] - p[3]*u[3]
end

println("invoking m1.get_data")
data = m1.get_data("20170825-beam.dat");
println("first three rows")
data[1:3,:]
println("typeof(data)", typeof(data))
data

u0 = [1.0,0.0,0.0]
tspan = (0.0,100.0)
p = [10.0,28.0,8/3]
prob = ODEProblem(parameterized_lorenz!,u0,tspan,p)
println("Solving...")
sol = solve(prob)

println("Plotting...")
plot(sol,vars=(1,2,3))

