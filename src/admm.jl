include("init.jl")

nfreq = 10
mydata = datacube[32:96,32:96,1:nfreq]
mypsf = psfcube[:,:,1:nfreq]
mypsfadj = flipdim(flipdim(mypsf,1),2)

# precompute

fty = cubefilter(mydata,mypsfadj)

# main admm loop

niter = 0
nbitermax = 20

tol1 = Float64[]
tol2 = Float64[]
loop = true

rhop = 1.0
taup


while loop
    tic()
    niter +=1

    # update x
    b = fty + taup + rhop*p
    
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
