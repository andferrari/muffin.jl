include("init.jl")
include("prox.jl")

nfreq = 10
#nfreq = size(psfcube)[3]
mydata = datacube[32:96,32:96,1:nfreq]
mypsf = float64(psfcube[:,:,1:nfreq])
mypsfadj = float64(flipdim(flipdim(mypsf,1),2))

# precompute

fty = cubefilter(mydata,mypsfadj)
nfty = size(fty)[1]

# main admm loop

niter = 0
nbitermax = 1000

tol1 = Float64[]
tol2 = Float64[]

loop = true

rhop = 0.10
muesp = 1.0
mu = muesp + rhop

taup = zeros(Float64,nfty,nfty,nfreq)
p = zeros(Float64,nfty,nfty,nfreq)
x = SharedArray(Float64,nfty,nfty,nfreq)
xmm = zeros(Float64,nfty,nfty,nfreq)

errorrec = zeros(Float64,nfty,nfty,nfreq)
errorest = zeros(Float64,nfreq)
errorraw = zeros(Float64,nfreq)

tic()
    while loop
        tic()
        niter +=1
        println("")
        println("ADMM iteration: $niter")

        # #update x

        ####################################
        ####################################
        # listmap_b =  { ( [fty[:,:,z] + taup[:,:,z] + rhop*p[:,:,z]] )  for z=1:nfreq}
        # listmap = { ( ([x[:,:,z]],[listmap_b[z]],[mypsf[:,:,z]],[mypsfadj[:,:,z]],mu) ) for z=1:nfreq}
        #
        # fmap(M) = conjgrad(M[1],M[2],M[3],M[4],M[5])
        # bla = pmap(fmap, listmap)
        #
        # @sync @parallel  for z = 1:nfreq
        #                     x[:,:,z] = bla[z]
        #                  end
        ####################################
        ####################################
        # x = sum(reshape(pmap(fmap, listmap),1,nfreq),2)

        ####################################
        ####################################
        @sync @parallel  for z = 1:nfreq
                           b = fty[:,:,z] + taup[:,:,z] + rhop*p[:,:,z]
                           x[:,:,z] = conjgrad(x[:,:,z],b,mypsf[:,:,z],mypsfadj[:,:,z],mu,tol=1e-4,itermax = 1e3)
                         end
        ####################################
        ####################################


        # prox positivity
        tmp = x-taup/rhop
        p = max(0,tmp)

        # update of Lagrange multipliers
        taup = taup + rhop*(p-x)

        # computer residues
        push!(tol1,vecnorm(x - xmm, 2))
        push!(tol2,vecnorm(x - p, 2))

        # stopping rule
        if (niter >= nbitermax) || ((tol1[niter] < 1E-3) && (tol2[niter] < 1E-2))
            loop = false
        end

        xmm[:] = x


        @printf("| - error ||x - xm||: %02.04e \n", tol1[niter])
        @printf("| - error ||x - xp||: %02.04e \n", tol2[niter])
        @printf("time for iteration : %f seconds \n", toq())

    end
println("")
@printf("time for ADMM : %f seconds \n", toq())

figure(2)
for z = 1:nfreq
    subplot(5,2,z)
    imshow(x[:,:,z])
end

figure(3)
for z = 1:nfreq
    errorrec[:,:,z] = cluster[32:96,32:96] - x[:,:,z]
    errorest[z] =  vecnorm(cluster[32:96,32:96] - x[:,:,z])/vecnorm(cluster[32:96,32:96])
    errorraw[z] =  vecnorm(mydata[:,:,z] - x[:,:,z])/vecnorm(mydata[:,:,z])
    subplot(5,2,z)
    imshow(errorrec[:,:,z])
end

figure(4)
plot([1:nfreq],errorest,color="blue")
plot([1:nfreq],errorraw,color="red")
