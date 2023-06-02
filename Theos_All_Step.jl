function react()

    return @reaction_network begin
        @parameters a1 a2 a3 σ1 σ2 σ3 tc1 tc2 tc3 k₁ k₋₁ k₂ k₃ k₄ k₅

        H2Pulse(t,[a1,a2,a3],[σ1,σ2,σ3],[tc1,tc2,tc3]), ∅ → H₂
        (k₁, k₋₁), H₂ + O  ↔ OH + H
        k₂,         H + OH → H₂O
        k₃,        OH + OH → H₂O
        k₄,         H +  H → H₂
        k₅,       H₂O      → ∅
    end

end

# function react(pump_pulse::Function)

#     return @reaction_network begin

#         pump_pulse(t), ∅ → H₂
#         (k₁, k₋₁), H₂ + O  ↔ OH + H
#         k₂,         H + OH → H₂O
#         k₃,        OH + OH → H₂O
#         k₄,         H +  H → H₂
#         k₅,       H₂O      → ∅
#     end

# end
# latexify(reqs, env=:chemical; mathjax = true) |> render
# latexify(reqs,form=:ode) |> render

#   p = (:t₀=>1.0, :σ=>.1, :k₁ => 1.0, :kᵦ => 1e-3,:k₂ => 1.0,:k₃ => 10.0,:k₄ => 1.0,:k₅ => 100.0)
#   tspan = (0.,10.)
#   u0 = [:H₂ => 0.,:O => 1., :H => 0.0, :OH => 0.0, :H₂O => 0.0]
  
# # solve ODEs
# oprob = ODEProblem(reqs, u0, tspan, p)
# osol  = solve(oprob)

# plot(osol; title = "Theo's All Step Mechanism")
