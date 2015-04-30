using HDF5, JLD

save("../data/results/struct_2g_2000iter.jld","nfreq",nfreq,
                                "lastiter",lastiter,
                                "psfst",psfst,
                                "skyst",skyst,
                                "admmst",admmst,
                                "toolst",toolst)
