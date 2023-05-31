using Catalyst, DifferentialEquations, Plots
using Latexify

s(t,t₀,σ) = exp(-(t - t₀)^2/(2*σ^2))/(σ*sqrt(2*π))

reqs = @reaction_network begin
  @parameters t₀ σ
  @species A(t)
    s(t,t₀,σ), 0 --> A
end

latexify(reqs, env=:chemical) |> render
latexify(reqs,form=:ode) |> render

p = (:t₀=>1.0, :σ=>.1)
tspan = (0.,4.)
u0 = [:A => 0.0]
  
# solve ODEs
oprob = ODEProblem(reqs, u0, tspan, p)
osol  = solve(oprob)

plot(osol; title = "Adsorption")
