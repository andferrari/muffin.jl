
####################################################################
#######                     MUFFIN                           #######
####### Multi-Frequency Sparse Radio Interferometric imaging #######
####################################################################

################################

# psfst, skyst, algost, admmst, toolst = muffin(nitermax = 30);

# demo = string(Pkg.dir("Muffin"),"/src/muffindemo.jl")

################################


module Muffin
    using FITSIO
    using Images
    using Wavelets
    using GHF


    export muffin,savedata

    include("muffintype.jl")

    #export admmst


    include("muffinfunc.jl")
    include("muffinload.jl")
    include("savedata.jl")

end


################################

# Pkg.clone("https://github.com/andferrari/GHF.jl.git")
# Pkg.add("PyPlot")
# Pkg.add("FITSIO")
# Pkg.add("Images")
# Pkg.add("Wavelets")
# Pkg.add("HDF5")

################################
