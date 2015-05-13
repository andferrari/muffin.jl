
                 ##################################
    ###################     Muffin Demo     #################
                 ##################################

################################

folder = string(Pkg.dir("Muffin"))
dataobj =  "data/M31.fits"
datapsf =  "data/meerkat_m30_25pix.psf.fits"
nitermax = 10
rhop = 1
rhot = 5
rhov = 2
rhos = 1
μt = 5e-1
μv = 1
mueps = 1e-3
savepath = string(folder,folder[1],"data/results/demo_results.jld")

################################


psfst, skyst, algost, admmst, toolst = muffin(folder,dataobj,datapsf,nitermax,
                                              rhop,rhot,rhov,rhos,μt,μv,mueps)


################################

savedata(savepath,psfst, skyst, algost, admmst, toolst)

################################
