include("../prox.jl")

n = 64
p = 21
mu = 3.

image = rand(n,n).^2
mypsf = randn(p,p).^2
mypsfadj = flipdim(flipdim(mypsf,1),2)

b = imfilter_fft(imfilter_fft(image,mypsf,"circular"),mypsfadj,"circular")
b = b + 3*image

x0 = zeros(n,n)
@time xest = conjgrad(x0,b,mypsf,mypsfadj,mu;tol = 1e-5,itermax = 1e3)
Qres = imfilter_fft(imfilter_fft(xest,mypsf,"circular"),mypsfadj,"circular")
Qres = Qres + 3*xest
println(norm(Qres-b))
