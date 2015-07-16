    ####################################################################
#######                Main Muffin function                  #######
####################################################################

function muffin(;folder="",dataobj="",datapsf="",nitermax = 500, rhop = 1,
                rhot = 5, rhov = 2, rhos = 1, μt = 5e-1, μv = 1e-0, mueps = 1e-3,
                bw = 5, ws="",parallel="")

println("")
println("MUFFIN initialisation")

                 ##################################
    ################### data initialisation #################
                 ##################################


    if typeof(dataobj) == ASCIIString
        if dataobj == "m31"
            psf = "data/meerkat_m30_25pix.psf.fits"
            obj = "data/M31.fits"
            tmp = string(Pkg.dir("Muffin"))
            psf = string(tmp,tmp[1],psf)
            obj = string(tmp,tmp[1],obj)
        elseif dataobj == "andro"
            psf = "data/meerkat_m30_25pix.psf.fits"
            obj = "data/andro.fits"
            tmp = string(Pkg.dir("Muffin"))
            psf = string(tmp,tmp[1],psf)
            obj = string(tmp,tmp[1],obj)
        elseif dataobj == "2gauss"
            psf = "data/meerkat_m30_25pix.psf.fits"
            obj = "data/2gauss.fits"
            tmp = string(Pkg.dir("Muffin"))
            psf = string(tmp,tmp[1],psf)
            obj = string(tmp,tmp[1],obj)
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
    println("loading psf...")
    psfst = loadpsf_dirty(psf)
    println("loading sky...")
    skyst = loadsky_dirty(obj,psfst.mypsf,psfst.nu)
    ##################################
else
    ##################################
    println("loading psf...")
    psfst = loadpsf(psf,bw)
    println("loading sky...")
    skyst = loadsky(obj,psfst.mypsf,psfst.nu)
    ##################################
end


    ##################################
    spatialwlt  = [WT.db1,WT.db2,WT.db3,WT.db4,WT.db5,WT.db6,WT.db7,WT.db8]
    # spatialwlt  = [WT.db1,WT.db2,WT.db3,WT.db4,WT.db5,WT.db6,WT.db7,WT.db8,WT.haar]

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
        file = string("/home/deguignet/resultats_400ite_woshar.jld")
        # file = string("/Users/deguignet/test.jld")
        admmst = loadarray_ws(rhop,rhot,rhov,rhos,μt,μv,mueps,nspat,nfreq,nxy,
                            skyst.mydata,psfst.mypsfadj,file)
    algost.niter = load(file,"algost.lastiter")
    toolst = loadtools(nitermax,nfreq,nxy)
    toolst.tol1 = load(file,"toolst.tol1")
    toolst.tol2 = load(file,"toolst.tol2")
    toolst.tol3 = load(file,"toolst.tol3")
    toolst.tol4 = load(file,"toolst.tol4")
    toolst.tol5 = load(file,"toolst.tol5")
    skyst.mydata = load(file,"skyst.mydata")


    else
        admmst = loadarray(rhop,rhot,rhov,rhos,μt,μv,mueps,nspat,nfreq,nxy,
                            skyst.mydata,psfst.mypsfadj,parallel)
        toolst = loadtools(nitermax,nfreq,nxy)
    end


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

    niter = algost.niter

    loop = true

    println("")
    println("                                =======================                                ")
    println("                               | MUFFIN RECONSTRUCTION |                               ")

    println("========================================================================================")
    println("ADMM"," ","|"," ","||x - xm||"," ","|"," ","||x - xp||"," ","|"," ",
            "||Hx - t||"," ","|"," ","||x - s||"," "," |"," ","||sh - v||"," ","|"," ","     time      |")
    println("========================================================================================")
    tic()
        while loop
            tic()
            niter +=1
            println("    "," ","|"," ","          "," ","|"," ","          "," ","|"," ",
                    "          "," ","|"," ","         "," "," |"," ","          "," ","|"," ")
            # println("ADMM iteration: $niter")

    ######################################
    ######################################
            # tic()
            # @parallel for i in 2:nfreq+1
            #     z = i-1
            #     admmst.wlt[:,:,z],admmst.x[:,:,z],admmst.t[:,:,z,:],admmst.taut[:,:,z,:],admmst.p[:,:,z],admmst.taup[:,:,z] =
            #                                 @fetchfrom(i,parallelmuffin(admmst.wlt[:,:,z], admmst.taut[:,:,z,:], admmst.t[:,:,z,:], rhot, admmst.x[:,:,z],
            #                                 psfst.mypsf[:,:,z], admmst.p[:,:,z], admmst.taup[:,:,z],
            #                                 fty[:,:,z], rhop, admmst.taus[:,:,z], admmst.s[:,:,z], rhos, admmst.mu, spatialwlt,
            #                                 μt, nspat))
########################

            @sync @parallel for z in 1:nfreq
                admmst.wlt[:,:,z],admmst.x[:,:,z],admmst.t[:,:,z,:],admmst.taut[:,:,z,:],admmst.p[:,:,z],admmst.taup[:,:,z] =
                                            parallelmuffin(admmst.wlt[:,:,z], admmst.taut[:,:,z,:], admmst.t[:,:,z,:], rhot, admmst.x[:,:,z],
                                            psfst.mypsf[:,:,z], admmst.p[:,:,z], admmst.taup[:,:,z],
                                            fty[:,:,z], rhop, admmst.taus[:,:,z], admmst.s[:,:,z], rhos, admmst.mu, spatialwlt,
                                            μt, nspat)

            end






            ##############################
            ######### prox spec ##########

            tmp = admmst.tauv + rhov*admmst.v

            admmst.s, admmst.sh = estime_ssh(admmst.s,admmst.sh,tmp,nxy,nspec,admmst.spectralwlt,
                                              admmst.x,admmst.taus,rhov,rhos)

            tmp = admmst.sh - admmst.tauv/rhov
            admmst.v = prox_u(tmp,μv/rhov)


            ########################################
            #### update of Lagrange multipliers ####

            admmst.tauv = admmst.tauv + rhov*(admmst.v-admmst.sh)
            admmst.taus = admmst.taus + rhos*(admmst.s-admmst.x)

    ######################################
    ######################################

            #
            # ##############################
            # ########## update x ##########
            # for z in 1:nfreq
            #     admmst.wlt[:,:,z] = myidwt((admmst.wlt)[:,:,z], nspat, (admmst.taut)[:,:,z,:], rhot,
            #                         (admmst.t)[:,:,z,:], spatialwlt)
            # end
            #
            # b = fty + admmst.taup + rhop*admmst.p + admmst.taus + rhos*admmst.s
            #
            # admmst.x = estime_x_par(admmst.x,psfst.mypsf,admmst.wlt + b,mu,nfreq)
            # ##############################
            # ######### prox spat ##########
            # tmp1 = 0.0
            # tmp2 = zeros(Float64,nxy,nxy)
            # for z in 1:nfreq
            #     for b in 1:nspat
            #             hx = dwt(admmst.x[:,:,z],wavelet(spatialwlt[b]))
            #             tmp = hx - admmst.taut[:,:,z,b]/rhot
            #             admmst.t[:,:,z,b] = prox_u(tmp,μt/rhot)
            #             admmst.taut[:,:,z,b] = admmst.taut[:,:,z,b] + rhot*(admmst.t[:,:,z,b]-hx)
            #             tmp1 = vecnorm([tmp2 (hx-(admmst.t)[:,:,z,b])],2)
            #             tmp2 = (hx-(admmst.t)[:,:,z,b])
            #     end
            # end
            # tmp2[:] = 0
            # ##############################
            # ###### prox positivity #######
            #
            # tmp = admmst.x-admmst.taup/rhop
            #
            # admmst.p = max(0,tmp)
            #
            # ##############################
            # ######### prox spec ##########
            #
            # tmp = admmst.tauv + rhov*admmst.v
            #
            # admmst.s, admmst.sh = estime_ssh(admmst.s,admmst.sh,tmp,nxy,nspec,admmst.spectralwlt,
            #                                   admmst.x,admmst.taus,rhov,rhos)
            #
            # tmp = admmst.sh - admmst.tauv/rhov
            # admmst.v = prox_u(tmp,μv/rhov)
            #
            #
            # ########################################
            # #### update of Lagrange multipliers ####
            #
            # admmst.taup = admmst.taup + rhop*(admmst.p-admmst.x)
            # admmst.tauv = admmst.tauv + rhov*(admmst.v-admmst.sh)
            # admmst.taus = admmst.taus + rhos*(admmst.s-admmst.x)

    ######################################
    ######################################


            ##############################
            ##### computer residues ######

            push!(toolst.tol1,vecnorm(admmst.x - admmst.xmm, 2)^2)
            push!(toolst.tol2,vecnorm(admmst.x - admmst.p, 2)^2)
            # push!(toolst.tol3,tmp1^2)
            push!(toolst.tol4,vecnorm(admmst.x - admmst.s, 2)^2)
            push!(toolst.tol5,vecnorm(admmst.sh - admmst.v, 2)^2)


            # ##############################
            # ####### stopping rule ########

            if (niter >= nitermax) #|| ((toolst.tol1[niter] < 1E-6) && (toolst.tol2[niter] < 1E-4))
                loop = false
                algost.lastiter = niter
            end

            admmst.xmm[:] = admmst.x

            #
            # println("ADMM"," ","|"," ","||x - xm||"," ","|"," ","||x - xp||"," ","|"," ",
            #         "||Hx - t||"," ","|"," ","||x - s||"," ","|"," ","||sh - v||"," ","|"," ","time")
            # println(niter," ","|"," ",toolst.tol1[niter]," ","|"," ",toolst.tol2[niter]," ","|"," ",
            #         toolst.tol3[niter]," ","|"," ",toolst.tol4[niter]," ","|"," ",toolst.tol5[niter]," ","|"," ",toq())



            @printf("%03d  | ", niter)
            @printf("%02.04e | ", toolst.tol1[niter])
            @printf("%02.04e | ", toolst.tol2[niter])
            # @printf("%02.04e | ", toolst.tol3[niter])
            @printf("%02.04e | ", 0)
            @printf("%02.04e | ", toolst.tol4[niter])
            @printf("%02.04e | ", toolst.tol5[niter])
            @printf("%f seconds  \n", toq())

        end
    println("")
    @printf("time for ADMM : %f seconds \n", toq())
    ####################################################################
    ####################################################################
    ####################################################################

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
    data = squeeze(data,find(([size(data)...].==1)))
    data = data[:,:,:]

    return data
end
##################################



####################################################################
#######                  Admm functions                      #######
####################################################################

#################################
######### x estimation ##########
# function estime_x_par(x::Array{Float64,3},mypsf::Array{Float64,3},
#                         wlt_b::Array{Float64,3},mu::Float64,nfreq::Int64)
#
#     nxy = (size(x))[1]
#     nxypsf = (size(mypsf))[1]
#     fftpsf = zeros(Complex64,nxy,nxy,nfreq)
#     psfcbe = zeros(Complex64,nxy,nxy,nfreq)
#     psfpad = zeros(Float64,nxy,nxy,nfreq)
#
#     for z in 1:nfreq
#         psfpad[1:nxypsf,1:nxypsf,z] = fftshift(mypsf[:,:,z])
#         psfcbe[:,:,z] = 1./(abs(fft(psfpad[:,:,z])).^2+mu)
#         x[:,:,z] = real(ifft(psfcbe[:,:,z].*fft(wlt_b[:,:,z])))
#     end
#
#     return x
# end
#
# function estime_x_par(x::SharedArray{Float64,3},mypsf::Array{Float64,3},
#                         wlt_b::SharedArray{Float64,3},mu::Float64,nfreq::Int64)
#
#     nxy = (size(x))[1]
#     nxypsf = (size(mypsf))[1]
#     fftpsf = zeros(Complex64,nxy,nxy,nfreq)
#     psfcbe = zeros(Complex64,nxy,nxy,nfreq)
#     psfpad = zeros(Float64,nxy,nxy,nfreq)
#
#     for z in 1:nfreq
#         psfpad[1:nxypsf,1:nxypsf,z] = fftshift(mypsf[:,:,z])
#         psfcbe[:,:,z] = 1./(abs(fft(psfpad[:,:,z])).^2+mu)
#         x[:,:,z] = real(ifft(psfcbe[:,:,z].*fft(wlt_b[:,:,z])))
#     end
#
#     return x
# end

#################################
###### proximity operators ######
function prox_u(u::Array,μ::Float64)
    return max(0, 1-μ./abs(u)).*u
end

function prox_u(u::SharedArray,μ::Float64)
    return max(0, 1-μ./abs(u)).*u
end

# function prox_u(u::Array,μ::Float64)
# 	res=zeros(size(u));
# 	x=u;
# 	ind = (abs(x).>(.75*μ^(2/3))) ;
# 	res[ind]= 2/3*x[ind].*(1+cos(2*pi/3-2/3*acos(μ/8*(abs(x[ind])/3).^(-3/2)))) ;
#    return res
# end

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

function estime_ssh(s::Array{Float64,3},sh::Array{Float64,3},tmp::Array{Float64,3},
                    nxy::Int64,nspec::Int64,spectralwlt::Array{Float64,3},
                    x::SharedArray{Float64,3},taus::Array{Float64,3},rhov::Float64,rhos::Float64)

    spectralwlt = idct(tmp,3)
    s = (spectralwlt + rhos*x - taus)/(rhov*nspec + rhos)
    sh = dct(s,3)

    return s,sh
end

function myidwt(wlt::Array{Float64,2},nspat::Int,taut::Array{Float64,4},rhot::Float64,t::Array{Float64,4},spatialwlt)
        wlt = idwt(taut[:,:,1,1] + rhot*t[:,:,1,1],wavelet(spatialwlt[1]))
            for b in 2:nspat
                wlt = wlt + idwt(taut[:,:,1,b] + rhot*t[:,:,1,b],wavelet(spatialwlt[b]))
            end
        return wlt
end
function myidwt(wlt,nspat::Int,taut,rhot::Float64,t,spatialwlt)
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

function parallelmuffin(wlt::Array{Float64,2},taut::Array{Float64,4},t::Array{Float64,4},rhot::Float64,
                        x::Array{Float64,2},psf::Array{Float64,2},p::Array{Float64,2},taup::Array{Float64,2},
                        fty::Array{Float64,2},rhop::Float64,taus::Array{Float64,2},s::Array{Float64,2},rhos::Float64,
                        mu::Float64,spatialwlt,μt::Float64,nspat::Int)

        wlt = myidwt(wlt, nspat, taut[:,:,1,:], rhot, t[:,:,1,:], spatialwlt)
        b = fty + taup + rhop*p + taus + rhos*s
        wlt_b = wlt + b

        nxy = (size(x))[1]
        nxypsf = (size(psf))[1]
        psfcbe = zeros(Complex64,nxy,nxy)
        psfpad = zeros(Float64,nxy,nxy)
        psfpad[1:nxypsf,1:nxypsf] = psf[:,:]
        psfcbe = 1./(abs(fft(psfpad)).^2+mu)
        x = real(ifft(psfcbe.*fft(wlt_b)))

        # tmp1 = 0.0
        # tmp2 = zeros(Float64,nxy,nxy)
        for b in 1:nspat
                    hx = dwt(x,wavelet(spatialwlt[b]))
                    tmp = hx - taut[:,:,1,b]/rhot
                    t[:,:,1,b] = prox_u(tmp,μt/rhot)
                    taut[:,:,1,b] = taut[:,:,1,b] + rhot*(t[:,:,1,b]-hx)
                    # tmp1 = vecnorm([tmp2 (hx-(t)[:,:,z,b])],2)
                    # tmp2 = (hx-(t)[:,:,z,b])
        end
        # tmp2[:] = 0


        tmp = x-taup/rhop
        p = max(0,tmp)
        taup = taup + rhop*(p-x)

    return wlt,x,t,taut,p,taup

end

function parallelmuffin(wlt::SharedArray{Float64,2},taut::SharedArray{Float64,4},t::SharedArray{Float64,4},rhot::Float64,
                        x::SharedArray{Float64,2},psf::Array{Float64,2},p::SharedArray{Float64,2},taup::SharedArray{Float64,2},
                        fty::Array{Float64,2},rhop::Float64,taus::Array{Float64,2},s::Array{Float64,2},rhos::Float64,
                        mu::Float64,spatialwlt,μt::Float64,nspat::Int)

        wlt = myidwt(wlt, nspat, taut[:,:,1,:], rhot, t[:,:,1,:], spatialwlt)
        b = fty + taup + rhop*p + taus + rhos*s
        wlt_b = wlt + b

        nxy = (size(x))[1]
        nxypsf = (size(psf))[1]
        psfcbe = zeros(Complex64,nxy,nxy)
        psfpad = zeros(Float64,nxy,nxy)
        psfpad[1:nxypsf,1:nxypsf] = psf[:,:]
        psfcbe = 1./(abs(fft(psfpad)).^2+mu)
        x = real(ifft(psfcbe.*fft(wlt_b)))

        # tmp1 = 0.0
        # tmp2 = zeros(Float64,nxy,nxy)
        for b in 1:nspat
                    hx = dwt(x,wavelet(spatialwlt[b]))
                    tmp = hx - taut[:,:,1,b]/rhot
                    t[:,:,1,b] = prox_u(tmp,μt/rhot)
                    taut[:,:,1,b] = taut[:,:,1,b] + rhot*(t[:,:,1,b]-hx)
                    # tmp1 = vecnorm([tmp2 (hx-(t)[:,:,z,b])],2)
                    # tmp2 = (hx-(t)[:,:,z,b])
        end
        # tmp2[:] = 0


        tmp = x-taup/rhop
        p = max(0,tmp)
        taup = taup + rhop*(p-x)

    return wlt,x,t,taut,p,taup

end
