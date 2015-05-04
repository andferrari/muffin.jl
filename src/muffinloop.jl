####################################################################
########################## main admm loop ##########################
####################################################################
function muffinadmm(psfst, skyst, algost, admmst, toolst)

    const rhop = admmst.rhop
    const rhot = admmst.rhot
    const rhov = admmst.rhov
    const rhos = admmst.rhos
    const μt = admmst.μt
    const μv = admmst.μv
    const muesp = admmst.muesp
    const tt = admmst.tt
    const mu = admmst.mu
    const nspat = algost.nspat
    const nfreq = algost.nfreq
    const nspec = algost.nspec
    const nxy = algost.nxy

    spatialwlt  = [WT.db1,WT.db2,WT.db3,WT.db4,WT.db5,WT.db6,WT.db7,WT.db8,WT.haar]

    niter = 0

    loop = true

    tic()
        while loop
            tic()
            niter +=1
            println("")
            println("ADMM iteration: $niter")

            ##############################
            ########## update x ##########
            tic()
            for z in 1:nfreq, b in 1:nspat
                admmst.wlt[:,:,z] = sum(idwt(admmst.taut[:,:,z,b] + rhot*admmst.t[:,:,z,b],wavelet(spatialwlt[b])),4)
            end
            a = toq()
            println("calcul wlt","  ",a)

            tic()
            b = admmst.fty + admmst.taup + rhop*admmst.p + admmst.taus + rhos*admmst.s
            a = toq()
            println("calcul b","  ",a)
            #admmst.x = forconjgrad(admmst.x, b, psfst.mypsf, psfst.mypsfadj, mu, admmst.wlt, nfreq)

            tic()
            @sync @parallel for z in 1:nfreq
                            c = b[:,:,z] + (admmst.wlt)[:,:,z]
                            (admmst.x)[:,:,z] = conjgrad((admmst.x)[:,:,z], c,
                             (psfst.mypsf)[:,:,z], (psfst.mypsfadj)[:,:,z], admmst.mu, tol=1e-4, itermax = 1e3)
                            end
            a = toq()
            println("calcul parallel","  ",a)
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
            tic()
            tmp = permutedims(admmst.tauv + rhov*admmst.v,[3,1,2])
            admmst.s = estime_s(admmst.s,tmp,nxy,nspec,admmst)
            admmst.sh = estime_sh(admmst.s,nxy,admmst)
            a = toq()
            println("calcul s sh","  ",a)

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
            @printf("SNR : %02.04e dB \n", toolst.snr[niter])


            ##############################
            ############ RMSE ############
            for z in 1:nfreq
                toolst.err[niter,z] = sqrt(sum((admmst.x[:,:,z] - skyst.sky[:,:,z]).^2)/skyst.sumsky2[z])
            end

            ##############################
            ####### stopping rule ########

            if (niter >= algost.nitermax) || ((toolst.tol1[niter] < 1E-6) && (toolst.tol2[niter] < 1E-4))
                loop = false
                algost.lastiter = niter
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

    for z in 1:nfreq
        toolst.errorrec[:,:,z] = skyst.sky[:,:,z] - admmst.x[:,:,z]
        toolst.errorest[z] =  vecnorm(skyst.sky[:,:,z] - admmst.x[:,:,z])^2/vecnorm(skyst.sky[:,:,z])^2
        toolst.errorraw[z] =  vecnorm(skyst.mydata[:,:,z] - admmst.x[:,:,z])^2/vecnorm(skyst.mydata[:,:,z])^2
    end

    return psfst, skyst, algost, admmst, toolst

end
