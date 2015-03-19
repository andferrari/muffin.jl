include("init.jl")
include("prox.jl")

nfreq = 10
mydata = datacube[32:96,32:96,1:nfreq]
mypsf = float64(psfcube[:,:,1:nfreq])
mypsfadj = float64(flipdim(flipdim(mypsf,1),2))

# precompute

fty = cubefilter(mydata,mypsfadj)

# main admm loop

niter = 0
nbitermax = 20

tol1 = Float64[]
tol2 = Float64[]
loop = true

nfty = size(fty)[1]

rhop = 1.0
muesp = 1.0
mu = muesp + rhop
taup = zeros(Float64,nfty,nfty,nfreq)


p = zeros(Float64,nfty,nfty,nfreq)
x = zeros(Float64,nfty,nfty,nfreq)
xmm = zeros(Float64,nfty,nfty,nfreq)


while loop
    niter +=1
    println("ADMM iteration: $niter")

    # update x
    for z = 1:nfreq
        b = fty[:,:,z] + taup[:,:,z] + rhop*p[:,:,z]
        x[:,:,z] = conjgrad(x[:,:,z],b,mypsf[:,:,z],mypsfadj[:,:,z],mu,tol=1e-4,itermax = 1e3)
    end

    # prox positivity
    tmp = x-taup/rhop
    p = max(0,tmp)

    # update of Lagrange multipliers
    taup = taup + rhop*(p-x)

    # computer residues
    push!(tol1,vecnorm(x - xmm, 2))
    push!(tol2,vecnorm(x - p, 2))

    # stopping rule
    if (niter >= nbitermax) || ((tol1[niter] < 1E-3) && (tol2[niter] < 1E-3))
        loop = false
    end

    xmm = x

    @printf("| - error ||x - xm||: %02.04e \n", tol1[niter])
    @printf("| - error ||x - xp||: %02.04e \n", tol2[niter])

end

figure(2)
for z = 1:nfreq
    subplot(5,2,z)
    imshow(x[:,:,z])
end
