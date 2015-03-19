d = 1000
A = randn(10*d,d)
A = A'*A
b = randn(d,1)
x = zeros(d,1)

tol = 1e-6

nb = vecnorm(b)
r = b - A*x
rm = r
p = r
iter = 0
loop = true

while loop
    iter += 1
    Qp = A*p

    alpha = vecnorm(r)^2/dot(vec(Qp),vec(p))
    x = x + alpha*p
    r = r - alpha*Qp
    betaa = (vecnorm(r)/vecnorm(rm))^2
    p = r + betaa*p

    crit = vecnorm(r)
    println(crit)
    if iter > 10000 || crit < tol
             loop = false
    end
    rm = r
end
return x
