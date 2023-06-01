using Catalyst, DifferentialEquations, Plots
using Latexify

f(t,t₀,σ::Vector{Float64}) = exp(-(t - t₀)^2/(2*σ[1]^2))/(σ[1]*sqrt(2*π))+σ[2]

reqs = @reaction_network begin
  @parameters t₀ σ[1:2] k₁ kᵦ k₂ k₃ k₄ k₅
#  @species A(t)
 
    f(t,t₀,σ), ∅ → H₂
    (k₁,kᵦ),   H₂ + O ↔ OH + H
    k₂,       H + OH → H₂O
    k₃,      OH + OH → H₂O
    k₄,       H +  H → H₂
    k₅,     H₂O      → 0
  end

latexify(reqs, env=:chemical; mathjax = true) |> render
latexify(reqs,form=:ode) |> render

  p = (:t₀=>1.0, :(σ[1])=>.1, :(σ[2])=>0.0, :k₁ => 1.0, :kᵦ => 1e-3,:k₂ => 1.0,:k₃ => 10.0,:k₄ => 1.0,:k₅ => 100.0)
#  p = (1.0,.1, 0.0, 1.0, 1e-3,1.0,10.0,1.0,100.0)
  tspan = (0.,10.)
  u0 = [:H₂ => 0.,:O => 1., :H => 0.0, :OH => 0.0, :H₂O => 0.0]
  
# solve ODEs
oprob = ODEProblem(reqs, u0, tspan, p)
osol  = solve(oprob)

plot(osol; title = "Theo's All Step Mechanism")
