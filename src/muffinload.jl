function loadsky(obj::ASCIIString,mypsf,nu::Array)
    skyst = init_SKY()
    skyst.obj = obj
    skyst.sky0 = (lecture(obj)/maximum(lecture(obj)))[112:143,112:143]
    skyst.sky,skyst.alpha = sky2cube(skyst.sky0,nu)
    skyst.skyconv = cubefilter(skyst.sky,mypsf)
    skyst.sig = sqrt(mean(skyst.skyconv.^2)/100)
    skyst.var = skyst.sig*skyst.sig
    skyst.noise = skyst.sig*randn(size(skyst.sky)[1],size(skyst.sky)[1],size(mypsf)[3])
    skyst.mydata = skyst.skyconv + skyst.noise
    for z in 1:length(nu)
        push!(skyst.sumsky2, sum(skyst.sky[:,:,z].^2))
    end
    return skyst
end

function loadpsf(psf::ASCIIString,M::Int)
    psfst = init_PSF()
    psfst.psf = psf
    psfst.psfcube = lecture(psfst.psf)
    psfst.psfavg = cubeaverage(psfst.psfcube,M)
    psfst.nu, psfst.nu0 = cubefreq(psfst.psf,psfst.psfcube,M)
    psfst.mypsf = cropcubexy(psfst.psfavg,31)
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

function loadarray(rhop,rhot,rhov,rhos,μt,μv,muesp,nspat,nfreq,nxy,mydata,mypsfadj)
    admmst = init_Admmarray()
    admmst.s = zeros(Float64,nxy,nxy,nfreq)
    admmst.taus = zeros(Float64,nxy,nxy,nfreq)
    admmst.sh = zeros(Float64,nxy,nxy,nfreq)
    admmst.taup = zeros(Float64,nxy,nxy,nfreq)
    admmst.p = zeros(Float64,nxy,nxy,nfreq)
    admmst.tauv = zeros(Float64,nxy,nxy,nfreq)
    admmst.v = zeros(Float64,nxy,nxy,nfreq)
    admmst.t = zeros(Float64,nxy,nxy,nfreq,nspat)
    admmst.taut = zeros(Float64,nxy,nxy,nfreq,nspat)
    admmst.wlt = SharedArray(Float64,nxy,nxy,nfreq)
    admmst.x = copy(mydata)
    admmst.Hx = SharedArray(Float64,nxy,nxy,nfreq,nspat)
    admmst.xmm = zeros(Float64,nxy,nxy,nfreq)
    admmst.spectralwlt = zeros(Float64,nxy,nxy,nfreq)
    admmst.fty = cubefilter(mydata,mypsfadj)
    admmst.rhop = rhop
    admmst.rhot = rhot
    admmst.rhov = rhov
    admmst.rhos = rhos
    admmst.μt = 5e-1
    admmst.μv = 1e-0
    admmst.muesp = 0.001
    admmst.tt = rhot*nspat
    admmst.mu = muesp + rhop + admmst.tt + rhos
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
