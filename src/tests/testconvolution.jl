n = 20
p = 5

image = rand(n,n).^2
mypsf = randn(p,p).^2
mypsfadj = flipdim(flipdim(mypsf,1),2)

b0 = imfilter_fft(image,mypsf,"circular")
b0t = imfilter_fft(image,mypsfadj,"circular")
b = imfilter_fft(imfilter_fft(image,mypsf,"circular"),mypsfadj,"circular")

n2 = n^2
H = zeros(Float64,n2,n2)
Ht = zeros(Float64,n2,n2)
Q = zeros(Float64,n2,n2)

iter = 1
for i = 1:n, j = 1:n
    x = zeros(n,n)
    x[j,i] = 1


    H[:,iter] = vec(imfilter_fft(x,mypsf,"circular"))
    Ht[:,iter] = vec(imfilter_fft(x,mypsfadj,"circular"))
    Q[:,iter] =  vec(imfilter_fft(imfilter_fft(x,mypsf,"circular"),mypsfadj,"circular"))

    iter += 1
end

println(norm(vec(b0) - H*vec(image)))
println(norm(vec(b0t) - Ht*vec(image)))
println(norm(vec(b) - Q*vec(image)))
println(norm(H-Ht'))
println(norm(Ht*H - Q))
