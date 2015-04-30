##################################
#### Structure initialisation ####
##################################

function muffin(;nitermax = 2000, rhop = 1, rhot = 2, rhov = 2, rhos = 1,
                 μt = 5e-1, μv = 1e-0, muesp = 1e-3)

    ##################################
    psf = "../data/meerkat_m30_25pix.psf.fits"
    obj = "../data/M31.fits"
    ##################################
    psfst = loadpsf(psf,5)
    skyst = loadsky(obj,psfst.mypsf,psfst.nu)
    ##################################


    ##################################
    spatialwlt  = [WT.db1,WT.db2,WT.db3,WT.db4,WT.db5,WT.db6,WT.db7,WT.db8,WT.haar]
    const nspat = length(spatialwlt)
    const nfreq = size(psfst.mypsf)[3]
    const nspec = 1
    const nxy = size(skyst.mydata)[1]

    niter = 0
    lastiter = 0
    #const nitermax = 2000

    algost = loadparam(nspat,nfreq,nspec,nxy,niter,lastiter,nitermax)
    ##################################

    ##################################


    admmst = loadarray(rhop,rhot,rhov,rhos,μt,μv,muesp,nspat,nfreq,nxy,
                        skyst.mydata,psfst.mypsfadj)
    toolst = loadtools(nitermax,nfreq,nxy)
    ##################################


    psfst, skyst, algost, admmst, toolst = muffinadmm(psfst, skyst, algost, admmst, toolst)



    return psfst, skyst, algost, admmst, toolst
end
