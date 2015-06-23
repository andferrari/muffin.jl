function loadsky(obj::ASCIIString,mypsf,nu::Array)

    skyst = init_SKY()
    skyst.obj = obj
    skyst.sky0 = lecture(obj)
    if length(size(skyst.sky0)) == 2
        skyst.sky0 = (lecture(obj)/maximum(lecture(obj)))
        skyst.sky,skyst.alpha = sky2cube(skyst.sky0,nu)
        skyst.skyconv = cubefilter(skyst.sky,mypsf)
        skyst.sig = sqrt(mean(skyst.skyconv.^2)/100)
        skyst.var = skyst.sig*skyst.sig
        skyst.noise = skyst.sig*randn(size(skyst.sky)[1],size(skyst.sky)[1],size(mypsf)[3])
        skyst.mydata = skyst.skyconv + skyst.noise
        for z in 1:length(nu)
            push!(skyst.sumsky2, sum(skyst.sky[:,:,z].^2))
        end
    elseif length(size(skyst.sky0)) == 3
        skyst = init_SKY()
        skyst.mydata = lecture(obj)
    end
    return skyst
end

function loadpsf(psf::ASCIIString,M::Int,npixpsf::Int)
    psfst = init_PSF()
    psfst.psf = psf
    psfst.psfcube = lecture(psfst.psf)
    psfst.psfavg = cubeaverage(psfst.psfcube,M)
    # psfst.nu, psfst.nu0 = cubefreq(psfst.psf,psfst.psfcube,M)
    psfst.nu, psfst.nu0 = cubefreqchiara(psfst.psf,psfst.psfcube,M)
    psfst.mypsf = cropcubexy(psfst.psfavg,npixpsf)
    psfst.mypsfadj = flipdim(flipdim(psfst.mypsf,1),2)
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
    # admmst.wlt = SharedArray(Float64,nxy,nxy,nfreq)
    println("wlt")
    admmst.wlt = zeros(Float64,nxy,nxy,nfreq)
    println("wlttmp")
    admmst.wlttmp = zeros(Float64,nxy,nxy,nfreq,nspat)
    # admmst.wlttmp = SharedArray(Float64,nxy,nxy,nfreq,nspat)
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
