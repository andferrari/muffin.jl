####################################################################
########################## main admm loop ##########################
####################################################################
println("Initialisation...")
include("init.jl");

loop = true

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
        println("nonzeros t","  ",length(nonzeros(t[:,:,1,1])))


        # prox positivity
        tmp = x-taup/rhop
        p = max(0,tmp)

        # prox spec
        tmp = permutedims(tauv + rhov*v,[3,1,2])
        s = estime_s(s,tmp)
        sh = estime_sh(s)

        tmp = sh - tauv/rhov
        v = prox_u(tmp,μv)
        println("nonzeros v","  ",length(nonzeros(v[:,:,1])))

        # update of Lagrange multipliers
        taup = taup + rhop*(p-x)
        taut = taut + rhot*(t-Hx)
        tauv = tauv + rhov*(v-sh)
        taus = taus + rhos*(s-x)

        # computer residues
        push!(tol1,vecnorm(x - xmm, 2)^2/vecnorm(x, 2)^2)
        push!(tol2,vecnorm(x - p, 2)^2/vecnorm(x, 2)^2)
        push!(tol3,vecnorm(Hx - t, 2)^2/vecnorm(t, 2)^2)
        push!(tol4,vecnorm(x - s, 2)^2/vecnorm(x, 2)^2)
        push!(tol5,vecnorm(sh - v, 2)^2/vecnorm(v, 2)^2)


        push!(snr,10*log(vecnorm(x[:,:,1])^2/vecnorm(sky[:,:,1]-x[:,:,1])^2))
        @printf("SNR: %02.04e dB \n", snr[niter+1])

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
        if (niter >= nbitermax) || ((tol1[niter] < 1E-6) && (tol2[niter] < 1E-4))
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

        #include("plot_admm.jl")
    end
println("")
@printf("time for ADMM : %f seconds \n", toq())
####################################################################
####################################################################
####################################################################

for z = 1:nfreq
    errorrec[:,:,z] = sky[:,:,z] - x[:,:,z]
    errorest[z] =  vecnorm(sky[:,:,z] - x[:,:,z])^2/vecnorm(sky[:,:,z])^2
    errorraw[z] =  vecnorm(mydata[:,:,z] - x[:,:,z])^2/vecnorm(mydata[:,:,z])^2
end
