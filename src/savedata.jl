using HDF5, JLD

save("../data/results/data.jld","nfreq",nfreq,
                                "lastiter",lastiter,
                                "x", x,
                                "errorest",errorest,
                                "errorraw",errorraw,
                                "errorrec",errorrec,
                                "err",err,
                                "tol1",tol1,
                                "tol2",tol2,
                                "tol3",tol3,
                                "tol4",tol4,
                                "tol5",tol5)
