include("init.jl")

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

rhop = 1.0
nfty = size(fty)[1]
taup = zeros(Float64,nfty,nfty,nfreq)

p = zeros(Float64,nfty,nfty,nfreq)
x = zeros(Float64,nfty,nfty,10)

while loop
    tic()
    niter +=1

    # update x


    # while vecnorm(x[:,:,j+1] - x[:,:,j], 2) > 1E-8
    #     j += 1
    #     alpha[:,:,j] = (r[:,:,j]*r[:,:,j])/((cubefilter(cubefilter(p[:,:,1], mypsf), mypsfadj))*p[:,:,j])
    #     x[:,:,j+1] = x[:,:,j] + alpha[:,:,j]*p[:,:,j]
    #     r[:,:,j+1] = r[:,:,j] - alpha[:,:,j]*(cubefilter(cubefilter(p[:,:,1], mypsf), mypsfadj))
    #     beta[:,:,j] = (r[:,:,j+1]*r[:,:,j+1])/(r[:,:,j]*r[:,:,j])
    #     p[:,:,j+1] = r[:,:,j+1] + beta[:,:,j]*p[:,:,j]
    #     println(vecnorm(x[:,:,j+1] - x[:,:,j], 2))
    # end




    println("up_x")


    # prox positivity
    tmp = x-taup/rhop
    p = max(0,tmp)


    # update of Lagrange multipliers
    println("up_lagr")

    taup = taup + rhop*(p-x)

    push!(tol1,vecnorm(x - xm, 2))
    push!(tol2,vecnorm(x - p, 2))


    if (niter >= nbitermax) || ((tol1[niter] < 1E-8) && (tol2[niter] < 1E-8))
        loop = false
    end
    println(niter)
    toc()
    println(niter,"  ",tol1[niter],"  ",tol2[niter])

end
