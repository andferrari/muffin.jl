using HDF5, JLD

function savedata(savepath, psfst, skyst, algost, admmst, toolst)

    JLD.save(savepath,"admmst.x",admmst.x,
                      "algost.lastiter",algost.lastiter,
                      "algost.nfreq",algost.nfreq,
                    #   "toolst.errorrec",toolst.errorrec,
                    #   "toolst.errorest",toolst.errorest,
                    #   "toolst.errorraw",toolst.errorraw,
                    #   "toolst.err",toolst.err,
                      "toolst.tol1",toolst.tol1,
                      "toolst.tol2",toolst.tol2,
                    #   "toolst.tol3",toolst.tol3,
                      "toolst.tol4",toolst.tol4,
                      "toolst.tol5",toolst.tol5,
                      "skyst.mydata",skyst.mydata,
                    #   "skyst.sky",skyst.sky,
                    #   "toolst.snr",toolst.snr,
                      "psfst.nu",psfst.nu,
                      "psfst.nu0",psfst.nu0)
                    #   "skyst.noise",skyst.noise)
end
