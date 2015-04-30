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
include("structure.jl")
include("tmp.jl")
include("func.jl")
include("prox.jl")
########

########
include("intro.jl")
include("admmloop.jl")
########







################################

# Pkg.clone("https://github.com/andferrari/GHF.jl.git")
# Pkg.add("PyPlot")
# Pkg.add("FITSIO")
# Pkg.add("Images")
# Pkg.add("Wavelets")
# Pkg.add("HDF5")

################################
