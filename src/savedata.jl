using HDF5, JLD

function savedata(psfst, skyst, algost, admmst, toolst)

    JLD.save("/home/deguignet/savetest_structure.jld","psfst",psfst,
                                                    "skyst",skyst,
                                                    "algost",algost,
                                                    "toolst",toolst,
                                                    "admmst.x",admmst.x
                                                    )






end
