    ####################################################################
#######                Main Muffin function                  #######
####################################################################

function muffin(;folder="",dataobj="",datapsf="",nitermax = 500, rhop = 1,
                rhot = 5, rhov = 2, rhos = 1, μt = 5e-1, μv = 1e-0, mueps = 1e-3,
                bw = 5, npixpsf = 255, ws="")

println("")
println("MUFFIN initialisation")

                 ##################################
    ################### data initialisation #################
                 ##################################


    if typeof(dataobj) == ASCIIString
        if dataobj == "m31"
            psf = "data/meerkat_m30_25pix.psf.fits"
            obj = "data/M31.fits"
        elseif dataobj == "andro"
            psf = "data/meerkat_m30_25pix.psf.fits"
            obj = "data/andro.fits"
        elseif dataobj == "2gauss"
            psf = "data/meerkat_m30_25pix.psf.fits"
            obj = "data/2gauss.fits"
        elseif dataobj == "chiara"
            psf = "/home/deguignet/Julia/example_sim_psf.fits"
            obj = "/home/deguignet/Julia/example_sim_dirty.fits"

        elseif isempty(folder)
            tmp = pwd()
            psf = string(tmp,tmp[1],datapsf)
            obj = string(tmp,tmp[1],dataobj)
        elseif typeof(folder) == ASCIIString
                 psf = string(folder,folder[1],datapsf)
                 obj = string(folder,folder[1],dataobj)
        else
                 error("data folder is not correct")
        end

    elseif isempty(dataobj)
        psf = "data/meerkat_m30_25pix.psf.fits"
        obj = "data/M31.fits"
    end

println("psf :"," ",psf)
println("obj :"," ",obj)

                ##################################
    ################# Structure initialisation #################
                ##################################

if dataobj == "chiara"
    ##################################
    println("mode 1")
    println("loading psf...")
    psfst = loadpsf_dirty(psf)
    println("loading sky...")
    skyst = loadsky_dirty(obj,psfst.mypsf,psfst.nu)
    ##################################
else
    ##################################
    println("mode 2")
    println("loading psf...")
    psfst = loadpsf(psf,bw)
    println("loading sky...")
    skyst = loadsky(obj,psfst.mypsf,psfst.nu)
    ##################################
end


    ##################################
    spatialwlt  = [WT.db1,WT.db2,WT.db3,WT.db4,WT.db5,WT.db6,WT.db7,WT.db8,WT.haar]
    const nspat = length(spatialwlt)
    const nfreq = size(psfst.mypsf)[3]
    const nspec = 1
    const nxy = size(skyst.mydata)[1]
    niter = 0
    lastiter = 0
    #################################

    #################################
    println("loading param...")
    algost = loadparam(nspat,nfreq,nspec,nxy,niter,lastiter,nitermax)

    println("loading arrays...")
    if ws == "true"
        file = string("/home/deguignet/Julia/resultats_100iter_woshaar.jld")
        println(size(skyst.mydata),"  ",size(psfst.mypsfadj))
        admmst = loadarray_ws(rhop,rhot,rhov,rhos,μt,μv,mueps,nspat,nfreq,nxy,
                            skyst.mydata,psfst.mypsfadj,file)
    else
        println(size(skyst.mydata),"  ",size(psfst.mypsfadj))
        admmst = loadarray(rhop,rhot,rhov,rhos,μt,μv,mueps,nspat,nfreq,nxy,
                            skyst.mydata,psfst.mypsfadj)
    end

    toolst = loadtools(nitermax,nfreq,nxy)
    ##################################

                 ##################################
    #####################  Main Admm Loop  #####################
                 ##################################
    println("Starting ADMM")
    #################################
    psfst, skyst, algost, admmst, toolst = muffinadmm(psfst, skyst, algost, admmst, toolst)
    #################################


    return psfst, skyst, algost, admmst, toolst
end



####################################################################
#######                   Main Admm Loop                     #######
####################################################################

function muffinadmm(psfst, skyst, algost, admmst, toolst)

    const rhop = admmst.rhop
    const rhot = admmst.rhot
    const rhov = admmst.rhov
    const rhos = admmst.rhos
    const μt = admmst.μt
    const μv = admmst.μv
    const mueps = admmst.mueps
    const tt = admmst.tt
    const mu = admmst.mu
    const nspat = algost.nspat
    const nfreq = algost.nfreq
    const nspec = algost.nspec
    const nxy = algost.nxy
    const fty = admmst.fty
    const nitermax = algost.nitermax

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
            for z in 1:nfreq
                admmst.wlt[:,:,z] = myidwt((admmst.wlt)[:,:,z], nspat, (admmst.taut)[:,:,z,:], rhot,
                                    (admmst.t)[:,:,z,:], spatialwlt)
            end
            a = toq()
            println("calcul wlt","  ",a)

            # tic()
            # @sync @parallel for z in 1:nfreq
            #                     admmst.wlttmp[:,:,z,:] = myidwt(admmst.wlttmp[:,:,z,:],nspat,(admmst.taut)[:,:,z,:],rhot,
            #                                             (admmst.t)[:,:,z,:],spatialwlt)
            #                 end
            # admmst.wlt = convert(Array,squeeze(sum(admmst.wlttmp,4),4))
            # a = toq()
            # println("calcul wlt","  ",a)

            tic()
            b = fty + admmst.taup + rhop*admmst.p + admmst.taus + rhos*admmst.s
            a = toq()
            println("calcul b","  ",a)

            tic()
            admmst.x = estime_x_par(admmst.x,psfst.mypsf,psfst.mypsfadj,admmst.wlt + b,mu,nfreq)
            a = toq()
            println("calcul parallel","  ",a)
            # b = 0
            ##############################
            ######### prox spat ##########

            ############################################################
            ############################################################
            # tic()
            # for z in 1:nfreq, b in 1:nspat
            #     admmst.Hx[:,:,z,b] = dwt(admmst.x[:,:,z],wavelet(spatialwlt[b]))
            # end
            # a = toq()
            # println("calcul HX", "  ", a)
            #
            # tic()
            # @time tmp = admmst.Hx - (admmst.taut)/rhot
            #
            # @time admmst.t = prox_u(tmp,μt/rhot)
            # # # tmp = 0
            ##############################
            ##############################
            tic()
            tmp1 = 0.0
            tmp2 = zeros(Float64,nxy,nxy)
            for z in 1:nfreq
                for b in 1:nspat
                        hx = dwt(admmst.x[:,:,z],wavelet(spatialwlt[b]))
                        tmp = hx - admmst.taut[:,:,z,b]/rhot
                        admmst.t[:,:,z,b] = prox_u(tmp,μt/rhot)
                        admmst.taut[:,:,z,b] = admmst.taut[:,:,z,b] + rhot*(admmst.t[:,:,z,b]-hx)
                        tmp1 = vecnorm([tmp2 (hx-(admmst.t)[:,:,z,b])],2)
                        tmp2 = (hx-(admmst.t)[:,:,z,b])
                end
            end

            # for z in 1:nfreq
            #     admmst.t, admmst.taut, tmp1 = myspat(admmst.x[:,:,z], admmst.t[:,:,z,:], admmst.taut[:,:,z,:],
            #                                         tmp1, tmp2, nspat, spatialwlt, rhot, μt)
            # end


            tmp2[:] = 0
            a = toq()
            println("new", "  ", a)
            tic()
            ############################################################
            ############################################################





            ##############################
            ###### prox positivity #######

            @time tmp = admmst.x-admmst.taup/rhop

            @time admmst.p = max(0,tmp)

            a = toq()
            println("calcul prox","  ", a)
            # tmp = 0

            ##############################
            ######### prox spec ##########
            tic()
            tmp = admmst.tauv + rhov*admmst.v

            admmst.s, admmst.sh = estime_ssh(admmst.s,admmst.sh,tmp,nxy,nspec,admmst.spectralwlt,
                                              admmst.x,admmst.taus,rhov,rhos)

            a = toq()
            println("calcul s sh","  ",a)
            # tmp = 0

            tmp = admmst.sh - admmst.tauv/rhov
            admmst.v = prox_u(tmp,μv/rhov)
            # tmp = 0

            ########################################
            #### update of Lagrange multipliers ####

            admmst.taup = admmst.taup + rhop*(admmst.p-admmst.x)
            # admmst.taut = admmst.taut + rhot*(admmst.t-admmst.Hx)
            admmst.tauv = admmst.tauv + rhov*(admmst.v-admmst.sh)
            admmst.taus = admmst.taus + rhos*(admmst.s-admmst.x)


            ##############################
            ##### computer residues ######

            push!(toolst.tol1,vecnorm(admmst.x - admmst.xmm, 2)^2)
            push!(toolst.tol2,vecnorm(admmst.x - admmst.p, 2)^2)
            # push!(toolst.tol3,vecnorm(admmst.Hx - admmst.t, 2)^2)
            push!(toolst.tol3,tmp1^2)
            push!(toolst.tol4,vecnorm(admmst.x - admmst.s, 2)^2)
            push!(toolst.tol5,vecnorm(admmst.sh - admmst.v, 2)^2)

            # push!(toolst.snr,10*log10(mean(cubefilter(admmst.x,psfst.mypsf).^2)/(skyst.sig)^2))
            # @printf("SNR : %02.04e dB \n", toolst.snr[niter])


            # ##############################
            # ############ RMSE ############
            # println(size(toolst.err)," ",size(admmst.x)," ",size(skyst.sky)," ",size(skyst.sumsky2))
            # for z in 1:nfreq
            #     toolst.err[niter,z] = sqrt(sum((admmst.x[:,:,z] - skyst.sky[:,:,z]).^2)/skyst.sumsky2[z])
            # end

            # ##############################
            # ####### stopping rule ########

            if (niter >= nitermax) #|| ((toolst.tol1[niter] < 1E-6) && (toolst.tol2[niter] < 1E-4))
                loop = false
                algost.lastiter = niter
            end

            admmst.xmm[:] = admmst.x
            # println(run(`free -g`))

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

    # for z in 1:nfreq
    #     toolst.errorrec[:,:,z] = skyst.sky[:,:,z] - admmst.x[:,:,z]
    #     toolst.errorest[z] =  vecnorm(skyst.sky[:,:,z] - admmst.x[:,:,z])^2/skyst.sumsky2[z]
    #     toolst.errorraw[z] =  vecnorm(skyst.mydata[:,:,z] - admmst.x[:,:,z])^2/vecnorm(skyst.mydata[:,:,z])^2
    # end
    # savedata("result_2048pix_100ite.jld", psfst, skyst, algost, admmst, toolst)
    return psfst, skyst, algost, admmst, toolst

end


####################################################################
#######                  Data functions                      #######
####################################################################


##################################
function cubefilter{T<:FloatingPoint}(imagecube::Array{T,3},psfcube::Array{T,3})
    nximg, nyimg, nfreq = size(imagecube)
    nxpsf, nypsf, nfreq = size(psfcube)
    rescube = Array(Float64,nximg,nyimg,nfreq)
    for k in 1:nfreq
        rescube[:,:,k] = imfilter_fft(imagecube[:,:,k],psfcube[:,:,k],"circular")
    end
    return rescube
end

function cubefilter{T<:FloatingPoint}(imagecube::SharedArray{T,3},psfcube::Array{T,3})
    nximg, nyimg, nfreq = size(imagecube)
    nxpsf, nypsf, nfreq = size(psfcube)
    rescube = Array(Float64,nximg,nyimg,nfreq)
    for k in 1:nfreq
        rescube[:,:,k] = imfilter_fft(imagecube[:,:,k],psfcube[:,:,k],"circular")
    end
    return rescube
end

function cubefilter{T<:FloatingPoint}(imagecube::Array{T,2},psfcube::Array{T,3})
    nximg, nyimg = size(imagecube)
    nxpsf, nypsf, nfreq = size(psfcube)
    rescube = Array(Float64,nximg,nyimg,nfreq)
    for k in 1:nfreq
        rescube[:,:,k] = imfilter_fft(imagecube,psfcube[:,:,k],"circular")
    end
    return rescube
end
##################################

##################################
function cubeaverage{T<:FloatingPoint}(imagecube::Array{T,3},M::Int)
    nxpsf, nypsf, nfreq = size(imagecube)
    if nfreq < M
        error("Not enough channels to average!")
    end
    nfreqavg = itrunc(nfreq/M)
    rescube = Array(Float64, nxpsf, nypsf, nfreqavg)
    for k in 1:nfreqavg
        rescube[:,:,k] = sum(imagecube[:,:,(k-1)*M+1:k*M], 3)/M
    end

    return rescube
end
##################################

##################################
function cubefreq(psf::ASCIIString,imagecube::Array,M::Int)
    nxpsf, nypsf, nfreq = size(imagecube)
    if nfreq < M
        error("Not enough channels to average!")
    end
    nfreqavg = itrunc(nfreq/M)
    file = FITS(psf)
    header = read_header(file[1])
    nustart = header["CRVAL4"]
    nustep = header["CDELT4"]
    nu0 = header["RESTFRQ"]
    nu = zeros(Float64,nfreqavg)
    for i in 1:nfreqavg
        nu[i] = nustart + (M-1)*nustep + (i-1)*M*nustep
    end
    close(file)
    return nu,nu0
end
##################################

##################################
function cropcubexy{T<:FloatingPoint}(imagecube::Array{T,3},M::Int)
    nxpsf, nypsf, nfreq = size(imagecube)
    if ~(nxpsf == nypsf)
        error("Image must be square!")
    end
    if nxpsf < M
        println("no crop")
        return imagecube
    end

    rescube = Array(Float64, M, M, nfreq)
    if (iseven(nxpsf) && iseven(M)) | (isodd(nxpsf) && isodd(M))
        nstart = (nxpsf-M)/2 +1
        nend = (nxpsf+M)/2
        rescube = imagecube[nstart:nend,nstart:nend,:]
    else
        error("size of image and crop must have same parity!")
    end
    return rescube
end
##################################

##################################
function sky2cube{T<:FloatingPoint}(sky::Array{T,2},nu::Array{T,1})
    nx, ny = size(sky)
    nbands = length(nu)
    corl = minimum([nx,ny])/5.0
    rho(x,y) = exp(-norm([x, y])/corl)
    tx = linspace(0,nx-1,nx)
    ty = linspace(0,ny-1,ny)
    field = genghf(tx, ty, rho)
    skycube = Array(Float64, nx, ny, nbands)

    sm, sM = extrema(sky)
    fm, fM = extrema(field)

    alpha = 0.8 + 0.4*(sky-sm)/(sM-sm) + 0.4*(field-fm)/(fM-fm)

    nu0 = (nu[end]+nu[1])/2
    for k in 1:nbands
        skycube[:,:,k] = sky.* (nu[k]/nu0) .^(-alpha)
    end
    return skycube,alpha
end
##################################

##################################
function lecture(directory::ASCIIString)
    file = FITS(directory)
    data = float64(read(file[1]))
    close(file)
    println("taille data lecture"," ",size(data))
    data = squeeze(data,find(([size(data)...].==1)))

    data = data[:,:,:]
    println("taille tronquee","  ",size(data))
    return data
end
##################################



####################################################################
#######                  Admm functions                      #######
####################################################################

#################################
######### x estimation ##########
function estime_x_par(x::Array{Float64,3},mypsf::Array{Float64,3},mypsfadj::Array{Float64,3},
                        wlt_b::Array{Float64,3},mu::Float64,nfreq::Int64)

    nxy = (size(x))[1]
    nxypsf = (size(mypsf))[1]
    fftpsf = zeros(Complex64,nxy,nxy,nfreq)
    psfcbe = zeros(Complex64,nxy,nxy,nfreq)
    psfpad = zeros(Float64,nxy,nxy,nfreq)

    for z in 1:nfreq
        psfpad[1:nxypsf,1:nxypsf,z] = mypsf[:,:,z]
        psfcbe[:,:,z] = 1./(abs(fft(psfpad[:,:,z])).^2+mu)

        ########################
        x[:,:,z] = real(ifft(psfcbe[:,:,z].*fft(wlt_b[:,:,z])))
        ########################

        ########################
        # x[:,:,z] = imfilter_fft(wlt_b[:,:,z],ifftshift(ifft(psfcbe[:,:,z])),"circular")
        # xtmp = copy(x[:,:,z])
        # x[:,2:nxy,z] = xtmp[:,1:nxy-1]
        # x[:,1,z] = xtmp[:,nxy]
        # xtmp = copy(x[:,:,z])
        # x[2:nxy,:,z] = xtmp[1:nxy-1,:]
        # x[1,:,z] = xtmp[nxy,:]
        ########################
    end



    return x
end

#################################
###### proximity operators ######
# function prox_u(u::SharedArray,μ::Float64)
#     return (max(1-μ./abs(u),0).*u)
# end

function prox_u(u::Array,μ::Float64)
    return max(0, 1-μ./abs(u)).*u
end
#################################
####### s / sh estimation #######
function estime_ssh(s::Array{Float64,3},sh::Array{Float64,3},tmp::Array{Float64,3},
                    nxy::Int64,nspec::Int64,spectralwlt::Array{Float64,3},
                    x::Array{Float64,3},taus::Array{Float64,3},rhov::Float64,rhos::Float64)

    spectralwlt = idct(tmp,3)
    s = (spectralwlt + rhos*x - taus)/(rhov*nspec + rhos)
    sh = dct(s,3)

    return s,sh
end

function myidwt(wlt,nspat,taut,rhot,t,spatialwlt)
        wlt = idwt(taut[:,:,1,1] + rhot*t[:,:,1,1],wavelet(spatialwlt[1]))
            for b in 2:nspat
                wlt = wlt + idwt(taut[:,:,1,b] + rhot*t[:,:,1,b],wavelet(spatialwlt[b]))
            end
        return wlt
end
##################################

##################################
function cubefreqchiara(nfrequencies::Int)
    nfreq = nfrequencies
    nustart = 9.85e8
    nustep = 2e6

    nu = zeros(Float64,nfreq)
    for i in 1:nfreq
        nu[i] = nustart + (i-1)*nustep
    end
    nu0 = (nu[1] + nu[nfreq])/2

    return nu,nu0
end


function myspat(x,t,taut,tmp1,tmp2, nspat, spatialwlt, rhot, μt)
    println(size(x)," ",size(t)," ",size(taut))
    for b in 1:nspat
        hx = dwt(x[:,:], wavelet(spatialwlt[b]))
        tmp = hx - taut[:,:,1,b]/rhot
        t[:,:,1,b] = prox_u(tmp,μt/rhot)
        taut[:,:,1,b] = taut[:,:,1,b] + rhot*(t[:,:,1,b]-hx)
        tmp1 = vecnorm([tmp2 (hx-t[:,:,1,b])],2)
        tmp2 = (hx-t[:,:,1,b])
    end
    return t,taut,tmp1
end
