using HDF5, JLD

function savedata(savepath, psfst, skyst, algost, admmst, toolst)

    JLD.save("imagerec.jld","admmst.x",admmst.x)

    JLD.save(savepath,"algost.lastiter",algost.lastiter,
                      "algost.nfreq",algost.nfreq,
                      "toolst.err",toolst.err,
                      "toolst.tol1",toolst.tol1,
                      "toolst.tol2",toolst.tol2,
                      "toolst.tol3",toolst.tol3,
                      "toolst.tol4",toolst.tol4,
                      "toolst.tol5",toolst.tol5,
                      "skyst.mydata",skyst.mydata,
                      "skyst.sky",skyst.sky,
                      "psfst.nu",psfst.nu,
                      "psfst.nu0",psfst.nu0,
                      "admmst",admmst)



                    #   "skyst.noise",skyst.noise)
end

#   "toolst.snr",toolst.snr,
#   "toolst.errorrec",toolst.errorrec,
#   "toolst.errorest",toolst.errorest,
#   "toolst.errorraw",toolst.errorraw,
# "admmst.x",admmst.x
# "admmst.t",admmst.t,
# "admmst.taut",admmst.taut,
# "admmst.v",admmst.v,
# "admmst.tauv",admmst.tauv,
# "admmst.s",admmst.s,
# "admmst.taus",admmst.taus,
# "admmst.p",admmst.p,
# "admmst.taup",admmst.taup
