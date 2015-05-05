using HDF5, JLD

function savedata(psfst, skyst, algost, admmst, toolst)


save("../data/results/savetest_structure.jld", "psfst",psfst,
                                             "skyst",skyst,
                                             "algost",algost,
                                             "admmst",admmst,
                                             "toolst",toolst)
