####################################################################
############################### Main ###############################
####################################################################

########
println("Initialisation...")
using FITSIO
@everywhere using Images
using Wavelets
using GHF
########

########
include("muffintype.jl")
include("muffinload.jl")
include("muffinfunc.jl")
include("muffinprox.jl")
########

########
include("muffininit.jl")
include("muffinloop.jl")
########


################################

# Pkg.clone("https://github.com/andferrari/GHF.jl.git")
# Pkg.add("PyPlot")
# Pkg.add("FITSIO")
# Pkg.add("Images")
# Pkg.add("Wavelets")
# Pkg.add("HDF5")

################################
