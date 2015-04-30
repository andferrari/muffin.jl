####################################################################
############################## MUFFIN ##############################
####### Multi-Frequency Sparse Radio Interferometric imaging #######
####################################################################

module Muffin
    using FITSIO
    @everywhere using Images
    using Wavelets
    using GHF


    export muffin

    include("muffintype.jl")
    include("muffinload.jl")
    include("muffinfunc.jl")
    include("muffinprox.jl")
    include("muffinloop.jl")
    include("muffininit.jl")
    
end


# psfst, skyst, algost, admmst, toolst = muffin()




################################

# Pkg.clone("https://github.com/andferrari/GHF.jl.git")
# Pkg.add("PyPlot")
# Pkg.add("FITSIO")
# Pkg.add("Images")
# Pkg.add("Wavelets")
# Pkg.add("HDF5")

################################
