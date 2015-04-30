####################################################################
########################## main admm loop ##########################
####################################################################

loop = true

tic()
    while loop
        tic()
        niter +=1
        println("")
        println("ADMM iteration: $niter")

        ##############################
        ########## update x ##########

        for z = 1:nfreq, b = 1:nspat
            admmst.wlt[:,:,z] = sum(idwt(admmst.taut[:,:,z,b] + rhot*admmst.t[:,:,z,b],wavelet(spatialwlt[b])),4)
        end

        b = admmst.fty + admmst.taup + rhop*admmst.p + admmst.taus + rhos*admmst.s

        admmst.x = forconjgrad(admmst.x, b, psfst.mypsf, psfst.mypsfadj, mu, admmst.wlt, nfreq)


        ##############################
        ######### prox spat ##########

        for z in 1:nfreq, b in 1:nspat
            admmst.Hx[:,:,z,b] = dwt(admmst.x[:,:,z],wavelet(spatialwlt[b]))
        end

        tmp = admmst.Hx - admmst.taut/rhot
        admmst.t = prox_u(tmp,μt/rhot)


        ##############################
        ###### prox positivity #######

        tmp = admmst.x-admmst.taup/rhop
        admmst.p = max(0,tmp)


        ##############################
        ######### prox spec ##########

        tmp = permutedims(admmst.tauv + rhov*admmst.v,[3,1,2])
        admmst.s = st_estime_s(admmst.s,tmp)
        admmst.sh = st_estime_sh(admmst.s)

        tmp = admmst.sh - admmst.tauv/rhov
        admmst.v = prox_u(tmp,μv/rhov)


        ########################################
        #### update of Lagrange multipliers ####

        admmst.taup = admmst.taup + rhop*(admmst.p-admmst.x)
        admmst.taut = admmst.taut + rhot*(admmst.t-admmst.Hx)
        admmst.tauv = admmst.tauv + rhov*(admmst.v-admmst.sh)
        admmst.taus = admmst.taus + rhos*(admmst.s-admmst.x)


        ##############################
        ##### computer residues ######

        push!(toolst.tol1,vecnorm(admmst.x - admmst.xmm, 2)^2)
        push!(toolst.tol2,vecnorm(admmst.x - admmst.p, 2)^2)
        push!(toolst.tol3,vecnorm(admmst.Hx - admmst.t, 2)^2)
        push!(toolst.tol4,vecnorm(admmst.x - admmst.s, 2)^2)
        push!(toolst.tol5,vecnorm(admmst.sh - admmst.v, 2)^2)

        push!(toolst.snr,10*log10(mean(cubefilter(admmst.x,psfst.mypsf).^2)/(skyst.sig)^2))
        @printf("SNR : %02.04e dB \n", toolst.snr[niter+1])


        ##############################
        ############ RMSE ############
        for z = 1:nfreq
            toolst.err[niter,z] = sqrt(sum((admmst.x[:,:,z] - skyst.sky[:,:,z]).^2)/skyst.sumsky2[z])
        end

        ##############################
        ####### stopping rule ########

        if (niter >= nitermax) || ((toolst.tol1[niter] < 1E-6) && (toolst.tol2[niter] < 1E-4))
            loop = false
            lastiter = niter
        end

        admmst.xmm[:] = admmst.x


        @printf("| - error ||x - xm||: %02.04e \n", toolst.tol1[niter])
        @printf("| - error ||x - xp||: %02.04e \n", toolst.tol2[niter])
        @printf("| - error ||Hx - t||: %02.04e \n", toolst.tol3[niter])
        @printf("| - error ||x - s||: %02.04e \n", toolst.tol4[niter])
        @printf("| - error ||sh - v||: %02.04e \n", toolst.tol5[niter])

        @printf("time for iteration : %f seconds \n", toq())

    end
println("")
@printf("time for ADMM : %f seconds \n", toq())
####################################################################
####################################################################
####################################################################

for z = 1:nfreq
    toolst.errorrec[:,:,z] = skyst.sky[:,:,z] - admmst.x[:,:,z]
    toolst.errorest[z] =  vecnorm(skyst.sky[:,:,z] - admmst.x[:,:,z])^2/vecnorm(skyst.sky[:,:,z])^2
    toolst.errorraw[z] =  vecnorm(skyst.mydata[:,:,z] - admmst.x[:,:,z])^2/vecnorm(skyst.mydata[:,:,z])^2
end
