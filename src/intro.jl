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
admmst.fty = cubefilter(skyst.mydata,psfst.mypsfadj)
push!(toolst.snr,snr0)
admmst.x[:] = skyst.mydata
