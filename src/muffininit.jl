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
    skyst = loadsky(obj,psfst.nu)
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
    # const rhop = 1
    # const rhot = 5
    # const rhov = 2
    # const rhos = 1
    # const μt = 5e-1
    # const μv = 1e-0
    # const muesp = 0.001
    # const tt = rhot*nspat
    # const mu = muesp + rhop + tt + rhos

    admmst = loadarray(rhop,rhot,rhov,rhos,μt,μv,muesp)
    toolst = loadtools()
    ##################################


    include("muffinloop.jl")


    return psfst, skyst, algost, admmst, toolst
end
