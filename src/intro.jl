using FITSIO
@everywhere using Images
using Wavelets
using GHF

include("structure.jl")
include("tmp.jl")
include("func.jl")
include("prox.jl")

##################################
#### Structure initialisation ####
##################################

##################################
psf = "../data/meerkat_m30_25pix.psf.fits"
obj = "../data/M31.fits"
##################################
psfst = loadpsf(psf,5)
skyst = loadsky(obj,psfst.nu)
##################################
spatialwlt  = [WT.db1,WT.db2,WT.db3,WT.db4,WT.db5,WT.db6,WT.db7,WT.db8,WT.haar]
const nspat = length(spatialwlt)
const nfreq = size(psfst.mypsf)[3]
const nspec = 1
const nxy = size(skyst.mydata)[1]
niter = 0
lastiter = 0
const nitermax = 2000
# algost = loadparam(nspat,nfreq,nspec,nxy,niter,lastiter,nitermax)
##################################
const rhop = 1
const rhot = 5
const rhov = 2
const rhos = 1
const μt = 5e-1
const μv = 1e-0
const muesp = 0.001
const tt = rhot*nspat
const mu = muesp + rhop + tt + rhos

admmst = loadarray()
toolst = loadtools()
##################################

# precompute
snr0 = 10*log10(mean(cubefilter(skyst.sky,psfst.mypsf).^2)/(skyst.sig)^2)
fty = cubefilter(skyst.mydata,psfst.mypsfadj)
push!(toolst.snr,snr0)
admmst.x[:] = skyst.mydata



####################################################################
########################## main admm loop ##########################
####################################################################
println("Initialisation...")

loop = true

#figure(1)
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

        b = fty + admmst.taup + rhop*admmst.p + admmst.taus + rhos*admmst.s


        # xx = admmst.x
        # wwlt = admmst.wlt
        # pmypsf = psfst.mypsf
        # psfadj = psfst.mypsfadj

        # @sync @parallel  for z in 1:nfreq
        #                 b = fty[:,:,z] + admmst.taup[:,:,z] + rhop*admmst.p[:,:,z] + admmst.wlt[:,:,z] + admmst.taus[:,:,z] + rhos*admmst.s[:,:,z]
        #                 admmst.x[:,:,z] = conjgrad(admmst.x[:,:,z],b,psfst.mypsf[:,:,z],psfst.mypsfadj[:,:,z],mu,tol=1e-4,itermax = 1e3)
        #                 #xx[:,:,z] = conjgrad(xx[:,:,z], b[:,:,z] + admmst.wwlt[z], pmypsf[:,:,z], psfadj[:,:,z], mu, tol=1e-4, itermax = 1e3)
        #            end

        admmst.x = forconjgrad(admmst.x, b, psfst.mypsf, psfst.mypsfadj, mu, admmst.wlt, nfreq)

        # admmst.x = xx
        # prox spat
        for z in 1:nfreq, b in 1:nspat
            admmst.Hx[:,:,z,b] = dwt(admmst.x[:,:,z],wavelet(spatialwlt[b]))
        end

        tmp = admmst.Hx - admmst.taut/rhot
        admmst.t = prox_u(tmp,μt/rhot)


        # prox positivity
        tmp = admmst.x-admmst.taup/rhop
        admmst.p = max(0,tmp)

        # prox spec
        tmp = permutedims(admmst.tauv + rhov*admmst.v,[3,1,2])
        admmst.s = estime_s(admmst.s,tmp)
        admmst.sh = estime_sh(admmst.s)

        tmp = admmst.sh - admmst.tauv/rhov
        admmst.v = prox_u(tmp,μv/rhov)

        # update of Lagrange multipliers
        admmst.taup = admmst.taup + rhop*(admmst.p-admmst.x)
        admmst.taut = admmst.taut + rhot*(admmst.t-admmst.Hx)
        admmst.tauv = admmst.tauv + rhov*(admmst.v-admmst.sh)
        admmst.taus = admmst.taus + rhos*(admmst.s-admmst.x)

        # computer residues
        push!(toolst.tol1,vecnorm(admmst.x - admmst.xmm, 2)^2)
        push!(toolst.tol2,vecnorm(admmst.x - admmst.p, 2)^2)
        push!(toolst.tol3,vecnorm(admmst.Hx - admmst.t, 2)^2)
        push!(toolst.tol4,vecnorm(admmst.x - admmst.s, 2)^2)
        push!(toolst.tol5,vecnorm(admmst.sh - admmst.v, 2)^2)

        push!(toolst.snr,10*log10(mean(cubefilter(admmst.x,psfst.mypsf).^2)/(skyst.sig)^2))
        @printf("SNR : %02.04e dB \n", toolst.snr[niter+1])

        # plot
            for z = 1:nfreq
                toolst.err[niter,z] = sqrt(sum((admmst.x[:,:,z] - skyst.sky[:,:,z]).^2)/skyst.sumsky2[z])
            end
            #
            #
            #
            # clf()
            # for z = 1:nfreq
            #     subplot(5,2,z)
            #     plot(toolst.err[1:niter,z])
            # end

        # Spectre

        i = 128
        j = 108

        admmst.spectrex[:,niter] = (squeeze(squeeze(admmst.x[i,j,:],1),1))
        admmst.spectresky[:,niter] = (squeeze(squeeze(skyst.sky[i,j,:],1),1))

        # stopping rule
        if (niter >= nbitermax) || ((toolst.tol1[niter] < 1E-6) && (toolst.tol2[niter] < 1E-4))
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

        #include("plot_admm.jl")
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

include("savedata.jl")
