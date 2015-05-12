
####################################################################
#######                     MUFFIN                           #######
####### Multi-Frequency Sparse Radio Interferometric imaging #######
####################################################################

################################

# psfst, skyst, algost, admmst, toolst = muffin()

# psfst, skyst, algost, admmst, toolst = muffin(folder = "/home/deguignet/.julia/v0.3/Muffin/data",datapsf="meerkat_m30_25pix.psf.fits",dataobj="testobject.fits",nitermax=1 )
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
