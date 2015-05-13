
                 ##################################
    ###################     Muffin Demo     #################
                 ##################################

################################

myfolder = string(Pkg.dir("Muffin"))
mydataobj =  "data/M31.fits"
mydatapsf =  "data/meerkat_m30_25pix.psf.fits"
mynitermax = 1
myrhop = 1
myrhot = 5
myrhov = 2
myrhos = 1
myμt = 5e-1
myμv = 1
mymueps = 1e-3
mysavepath = string(myfolder,myfolder[1],"data/results/demo_results.jld")

################################


psfst, skyst, algost, admmst, toolst = muffin(folder = myfolder,
                                              dataobj = mydataobj,
                                              datapsf = mydatapsf,
                                              nitermax = mynitermax,
                                              rhop = myrhop,
                                              rhot = myrhot,
                                              rhov = myrhov,
                                              rhos = myrhos,
                                              μt = myμt,
                                              μv = myμv,
                                              mueps = mymueps)


################################

savedata(mysavepath, psfst, skyst, algost, admmst, toolst)

################################
