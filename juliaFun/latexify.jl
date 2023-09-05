using Catalyst, Latexify
rn = @reaction_network begin
    (a,b), 0 <--> X
end
latexify(rn; mathjax = true) |> render