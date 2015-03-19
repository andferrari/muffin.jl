n = 20
p = 5
mu = 3.

image = rand(n,n).^2
mypsf = randn(p,p).^2
mypsfadj = flipdim(flipdim(mypsf,1),2)

b = imfilter_fft(imfilter_fft(image,mypsf,"circular"),mypsfadj,"circular")
b = b + 3*image

conjgrad(zeros(n,n),b,mypsf,mypsfadj,mu;tol = 1e-6,itermax = 1e3);
