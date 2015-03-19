
using Images

n=100
n2=n*n
x = zeros(Float64,n,n)
xv = vec(x)
#psf = [0 2 1 ;1 2 1; 2 2 1]
psf = mypsf[:,:,1]
psfadj = flipdim(flipdim(psf,1),2)

H = zeros(Float64,n2,n2)
Ht = zeros(Float64,n2,n2)
Q = zeros(Float64,n2,n2)

iter = 1
for i = 1:n, j = 1:n
    x[:,:] = 0
    x[j,i] = 1


    H[:,iter] = vec(imfilter_fft(x,psf,"circular"))
    Ht[:,iter] = vec(imfilter_fft(x,psfadj,"circular"))
    Q[:,iter] =  vec(imfilter_fft(imfilter_fft(x,psf,"circular"),psfadj,"circular"))

    iter += 1
end
