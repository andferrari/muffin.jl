


##################################
#### Structure initialisation ####
##################################
psfst = init_PSF()
skyst = init_SKY()
algost = init_Algoparam()
admmst = init_Admmarray()
toolst = init_TOOLS()
##################################

psfst.psf = "../data/meerkat_m30_25pix.psf.fits"
skyst.obj = "../data/M31.fits"

function loadsky(skyst.obj::ASCIIString,nu::Array)
    skyst.sky0 = lecture(skyst.obj)/maximum(lecture(skyst.obj))
    skyst.sky,skyst.alpha = sky2cube(skyst.sky0,nu)
    skyst.skyconv = cubefilter(skyst.sky,skyst.mypsf)
    skyst.sig = sqrt(mean(skyst.skyconv.^2)/100)
    skyst.var = skyst.sig*skyst.sig
    skyst.noise = sig*randn(size(skyst.sky)[1],size(skyst.sky)[1],size(skyst.mypsf)[3])
    skyst.mydata = skyst.skyconv + skyst.noise
    skyst.sumsky2 = for z = 1:nw
                        sumsky2[z] = sum(sky[:,:,z].^2)
                    end
end

function loadpsf(psfst.psf::ASCIIString)
    psfst.psfcube = lecture(psfst.psf)
    psfst.psfabg = cubeaverage(psfst.psfcube,5)
    psfst.mypsf = cropcubexy(psfst.psfavg,255)
    psfst.mypsfadj = flipdim(flipdim(psfst.mypsf,1),2)
end

function loadparam(nspat,nfreq,nspec,nxy,niter,lastiter,nitermax)
    algost.nspat = nspat
    algost.nfreq = nfreq
    algost.nspec = nspec
    algost.nxy = nxy
    algost.niter = niter
    algost.lastiter = lastiter
    algost.nitermax = nitermax
end
