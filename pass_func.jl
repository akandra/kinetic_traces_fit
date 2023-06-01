myprod2(x,a)=a*x
mysum3(x,a,b)=x + a + b

f2(f::Function, x) = f(x)

println(f2( (x) -> myprod2(x, 5) , 2 ))

println(f2( (x) -> mysum3(x, 5, 1) , 2 ))