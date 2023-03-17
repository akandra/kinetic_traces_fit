function ann_par(name, p::fitpar) 
    # s = latexstring(name,L"=",string(round(p.value,sigdigits=3)))
    s = "\$"*name*" = "*string(round(p.value,sigdigits=3))*"\$"
    #s = @sprintf "\$ %.2f \\cdot 10^1\$" significand(p.value)#*2.0^exponent(p.value)
    text(s,:left, p.var ? :black : :gray)
end