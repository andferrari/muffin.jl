function loadsky(obj::ASCIIString,mypsf,nu::Array)

    skyst = init_SKY()
    sky0 = lecture(obj)
    if size(sky0)[3] == 1
        println("methode 1")
        sky0 = squeeze(sky0/maximum(sky0),3)
        skyst.sky,skyst.alpha = sky2cube(sky0,nu)
        skyconv = cubefilter(skyst.sky,mypsf)
        skyst.sig = sqrt(mean(skyconv.^2)/100)
        skyst.var = skyst.sig*skyst.sig
        skyst.noise = skyst.sig*randn( size(skyst.sky)[1], size(skyst.sky)[1], size(mypsf)[3])
        skyst.mydata = skyconv + skyst.noise
        for z in 1:length(nu)
            push!(skyst.sumsky2, sum(skyst.sky[:,:,z].^2))
        end
    elseif size(sky0)[3] != 1
        println("methode 2")

        skyst = init_SKY()
        skyst.sky = lecture(obj)
        skyst.mydata = cubefilter(skyst.sky,mypsf)
    end
    return skyst
end

function loadsky_dirty(obj::ASCIIString,mypsf,nu::Array)
    skyst = init_SKY_dirty()
    sky0 = "/home/deguignet/Julia/I_HALO_CL_RS_SKY.FITS"
    skyst.sky = lecture(obj)
    skyst.mydata = lecture(obj)
    for z in 1:length(nu)
        push!(skyst.sumsky2, sum(skyst.sky[:,:,z].^2))
    end
    return skyst
end


function loadpsf(psf::ASCIIString,M::Int)
    psfst = init_PSF()
    psfcube = lecture(psf)
    d = round((size(psfcube)[1])/2)
    # psfcube = psfcube[d-128:d+127,d-128:d+127,:]
    println(size(psfcube))
    psfavg = cubeaverage(psfcube,M)
    psfst.nu, psfst.nu0 = cubefreq(psf,psfcube,M)
    psfst.mypsf = cropcubexy(psfavg,size(psfcube)[1])
    psfst.mypsfadj = flipdim(flipdim(psfst.mypsf,1),2)
    return psfst
end

function loadpsf_dirty(psf::ASCIIString)
    psfst = init_PSF_dirty()
    psfcube = lecture(psf)
    psfst.nu, psfst.nu0 = cubefreqchiara(size(psfcube)[3])
    psfst.mypsf = cropcubexy(psfcube,size(psfcube)[1])
    psfst.mypsfadj = flipdim(flipdim(psfst.mypsf,1),2)
    println("size toto"," ",size(psfst.mypsf))
    return psfst
end

function loadparam(nspat,nfreq,nspec,nxy,niter,lastiter,nitermax)
    algost = init_Algoparam()
    algost.nspat = nspat
    algost.nfreq = nfreq
    algost.nspec = nspec
    algost.nxy = nxy
    algost.niter = niter
    algost.lastiter = lastiter
    algost.nitermax = nitermax
    return algost
end

function loadarray(rhop,rhot,rhov,rhos,μt,μv,mueps,nspat,nfreq,nxy,mydata,mypsfadj)
    println("init")
    admmst = init_Admmarray()
    println("s")
    admmst.s = zeros(Float64,nxy,nxy,nfreq)
    # admmst.s = SharedArray(Float64,nxy,nxy,nfreq)
    println("taus")
    admmst.taus = zeros(Float64,nxy,nxy,nfreq)
    println("sh")
    admmst.sh = zeros(Float64,nxy,nxy,nfreq)
    # admmst.sh = SharedArray(Float64,nxy,nxy,nfreq)
    println("taup")
    admmst.taup = zeros(Float64,nxy,nxy,nfreq)
    println("p")
    admmst.p = zeros(Float64,nxy,nxy,nfreq)
    println("tauv")
    admmst.tauv = zeros(Float64,nxy,nxy,nfreq)
    println("v")
    admmst.v = zeros(Float64,nxy,nxy,nfreq)
    println("t")
    admmst.t = zeros(Float64,nxy,nxy,nfreq,nspat)
    println("taut")
    admmst.taut = zeros(Float64,nxy,nxy,nfreq,nspat)

    println("wlt")
    admmst.wlt = zeros(Float64,nxy,nxy,nfreq)
    # admmst.wlt = SharedArray(Float64,nxy,nxy,nfreq)
    println("x")
    admmst.x = copy(mydata)
    # admmst.Hx = SharedArray(Float64,nxy,nxy,nfreq,nspat)
    println("Hx")
    admmst.Hx = zeros(Float64,nxy,nxy,nfreq,nspat)
    println("xmm")
    admmst.xmm = zeros(Float64,nxy,nxy,nfreq)
    println("spectralwlt")
    admmst.spectralwlt = zeros(Float64,nxy,nxy,nfreq)
    # admmst.spectralwlt = SharedArray(Float64,nxy,nxy,nfreq)
    println("fty")
    admmst.fty = cubefilter(mydata,mypsfadj)
    admmst.rhop = rhop
    admmst.rhot = rhot
    admmst.rhov = rhov
    admmst.rhos = rhos
    admmst.μt = μt
    admmst.μv = μv
    admmst.mueps = mueps
    admmst.tt = rhot*nspat
    admmst.mu = mueps + rhop + admmst.tt + rhos
    return admmst
end

function loadtools(nitermax,nfreq,nxy)
    toolst = init_TOOLS()
    toolst.snr = Float64[]
    toolst.tol1 = Float64[]
    toolst.tol2 = Float64[]
    toolst.tol3 = Float64[]
    toolst.tol4 = Float64[]
    toolst.tol5 = Float64[]
    toolst.err = zeros(Float64,nitermax,nfreq)
    toolst.errorrec = zeros(Float64,nxy,nxy,nfreq)
    toolst.errorest = zeros(Float64,nfreq)
    toolst.errorraw = zeros(Float64,nfreq)
    return toolst
end
