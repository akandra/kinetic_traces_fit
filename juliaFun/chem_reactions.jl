using Catalyst, DifferentialEquations, Plots
using Latexify

reqs = @reaction_network begin
    (k₁,kᵦ), H₂ +  O <--> OH + H
#    kᵦ, OH +  H --> H₂ + O
    k₂,  H + OH --> H₂O
    k₃, OH + OH --> H₂O
    k₄,  H +  H --> H₂
    k₅,  H₂O    --> 0
  end

latexify(reqs, env=:chemical) |> render
latexify(reqs,form=:ode) |> render

  p = (:k₁ => 1.0, :kᵦ => 1e-3,:k₂ => 1.0,:k₃ => 10.0,:k₄ => 1.0,:k₅ => 1.0)
  tspan = (0.,10.)
  u0 = [:H₂ => 1.,:O => 1., :H => 0.0, :OH => 0.0, :H₂O => 0.0]
  
# solve ODEs
oprob = ODEProblem(reqs, u0, tspan, p)
osol  = solve(oprob, Tsit5())

plot(osol; title = "Theo's All Step Mechanism")
