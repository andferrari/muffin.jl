include("init.jl")
include("prox.jl")


#nfreq = 1
nfreq = size(psfcube)[3]
mydata = datacube[:,:,1:nfreq]
mypsf = psfcube[:,:,1:nfreq]
mypsfadj = float64(flipdim(flipdim(mypsf,1),2))

spatialwlt  = [WT.db1,WT.db2,WT.db3,WT.db4,WT.db5,WT.db6,WT.db7,WT.db8,WT.haar]
nspat = length(spatialwlt)

nspec = 1


# precompute

fty = cubefilter(mydata,mypsfadj)
nfty = size(fty)[1]

spectralwlt = zeros(Float64,nfty,nfty,nfreq)



# main admm loop

niter = 0
lastiter = 0
nbitermax = 1000

tol1 = Float64[]
tol2 = Float64[]
tol3 = Float64[]
tol4 = Float64[]
tol5 = Float64[]



err = Array(Float64,nbitermax,nfreq)


loop = true

rhop = 2
rhot = 1
rhov = 5
rhos = 5
μt = 1.0
μv = 1.0
muesp = 1.0
tt = rhot*nspat
mu = muesp + rhop + tt +rhos

s = zeros(Float64,nfty,nfty,nfreq)
taus = zeros(Float64,nfty,nfty,nfreq)
sh = zeros(Float64,nfty,nfty,nfreq)

taup = zeros(Float64,nfty,nfty,nfreq)
p = zeros(Float64,nfty,nfty,nfreq)

tauv = zeros(Float64,nfty,nfty,nfreq)
v = zeros(Float64,nfty,nfty,nfreq)


t = zeros(Float64,nfty,nfty,nfreq,nspat)
taut = zeros(Float64,nfty,nfty,nfreq,nspat)
wlt = SharedArray(Float64,nfty,nfty,nfreq)

x = SharedArray(Float64,nfty,nfty,nfreq)
Hx = SharedArray(Float64,nfty,nfty,nfreq,nspat)
xmm = zeros(Float64,nfty,nfty,nfreq)

errorrec = zeros(Float64,nfty,nfty,nfreq)
errorest = zeros(Float64,nfreq)
errorraw = zeros(Float64,nfreq)

#figure(1)
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
        for z = 1:nfreq, b = 1:nspat
            wlt[:,:,z] = sum(idwt(taut[:,:,z,b] + rhot*t[:,:,z,b],wavelet(spatialwlt[b])),4)
        end

        @sync @parallel  for z = 1:nfreq
        println(z)
                           b = fty[:,:,z] + taup[:,:,z] + rhop*p[:,:,z] + wlt[z] + taus[:,:,z] + rhos*s[:,:,z]
                           x[:,:,z] = conjgrad(x[:,:,z],b,mypsf[:,:,z],mypsfadj[:,:,z],mu,tol=1e-4,itermax = 1e3)
                         end
        ####################################
        ####################################

        # prox spat
        for z = 1:nfreq, b = 1:nspat
            Hx[:,:,z,b] = dwt(x[:,:,z],wavelet(spatialwlt[b]))
        end

        tmp = Hx - taut/rhot
        t = prox_u(tmp,μt)

        # prox positivity
        tmp = x-taup/rhop
        p = max(0,tmp)

        # prox spec
        tmp = permutedims(tauv + rhov*v,[3,1,2])
        s = estime_s(s,tmp)
        sh = estime_sh(s)

        tmp = sh - tauv/rhov
        v = prox_u(tmp,μv)

        # update of Lagrange multipliers
        taup = taup + rhop*(p-x)
        taut = taut + rhot*(t-Hx)
        tauv = tauv + rhov*(v-sh)
        taus = taus + rhos*(s-x)




        # computer residues
        push!(tol1,vecnorm(x - xmm, 2)^2)
        push!(tol2,vecnorm(x - p, 2)^2)
        push!(tol3,vecnorm(Hx - t, 2)^2)
        push!(tol4,vecnorm(x - s, 2)^2)
        push!(tol5,vecnorm(sh - v, 2)^2)


        # plot
            for z = 1:nfreq
                err[niter,z] = vecnorm(x[:,:,z] - sky[:,:,z], 2)^2
            end
            #
            #
            #
            # clf()
            # for z = 1:nfreq
            #     subplot(5,2,z)
            #     plot(err[1:niter,z])
            # end


        # stopping rule
        if (niter >= nbitermax) || ((tol1[niter] < 1E-3) && (tol2[niter] < 1E-2))
            loop = false
            lastiter = niter
        end

        xmm[:] = x


        @printf("| - error ||x - xm||: %02.04e \n", tol1[niter])
        @printf("| - error ||x - xp||: %02.04e \n", tol2[niter])
        @printf("| - error ||Hx - t||: %02.04e \n", tol3[niter])
        @printf("| - error ||x - s||: %02.04e \n", tol4[niter])
        @printf("| - error ||sh - v||: %02.04e \n", tol5[niter])

        @printf("time for iteration : %f seconds \n", toq())

    end
println("")
@printf("time for ADMM : %f seconds \n", toq())

# figure(2)
# for z = 1:nfreq
#     subplot(5,2,z)
#     imshow(x[:,:,z])
# end
#
# figure(3)
for z = 1:nfreq
    errorrec[:,:,z] = sky[:,:,z] - x[:,:,z]
    errorest[z] =  vecnorm(sky[:,:,z] - x[:,:,z])^2/vecnorm(sky[:,:,z])^2
    errorraw[z] =  vecnorm(mydata[:,:,z] - x[:,:,z])^2/vecnorm(mydata[:,:,z])^2
    #subplot(5,2,z)
    #imshow(errorrec[:,:,z])
end

# figure(4)
# plot([1:nfreq],errorest,color="blue")
# plot([1:nfreq],errorraw,color="red")
