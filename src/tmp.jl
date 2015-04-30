function loadsky(obj::ASCIIString,nu::Array)
    skyst = init_SKY()
    skyst.obj = obj
    skyst.sky0 = lecture(obj)/maximum(lecture(obj))
    skyst.sky,skyst.alpha = sky2cube(skyst.sky0,nu)
    skyst.skyconv = cubefilter(skyst.sky,psfst.mypsf)
    skyst.sig = sqrt(mean(skyst.skyconv.^2)/100)
    skyst.var = skyst.sig*skyst.sig
    skyst.noise = skyst.sig*randn(size(skyst.sky)[1],size(skyst.sky)[1],size(psfst.mypsf)[3])
    skyst.mydata = skyst.skyconv + skyst.noise
    for z = 1:length(psfst.nu)
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
    psfst.mypsf = cropcubexy(psfst.psfavg,255)
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

function loadarray()
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
    admmst.x = SharedArray(Float64,nxy,nxy,nfreq)
    admmst.Hx = SharedArray(Float64,nxy,nxy,nfreq,nspat)
    admmst.xmm = zeros(Float64,nxy,nxy,nfreq)
    admmst.spectrex = zeros(Float64,nfreq,nitermax)
    admmst.spectresky = zeros(Float64,nfreq,nitermax)
    admmst.spectralwlt = zeros(Float64,nxy,nxy,nfreq)
    admmst.fty = zeros(Float64,nxy,nxy,nfreq)
    return admmst
end

function loadtools()
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
