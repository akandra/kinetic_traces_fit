function f(x,y)
a = 1/3
y[1] = x[1] * a
y[2] = x[2] * a
return x
end

function do1()
    out = zeros(2)
    [@time f([1,2],out) for _ in 1:3]
end

do1()