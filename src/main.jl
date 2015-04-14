####################################################################
############################### Main ###############################
####################################################################



########
include("func.jl")
include("prox.jl")
include("testobj.jl")
########


########
include("init.jl")
include("admm.jl")  # include("plot_admm.jl")
include("savedata.jl")
########


########
include("readdata.jl")
include("alpha_rec.jl")
########



################################

# Pkg.clone("https://github.com/andferrari/GHF.jl.git")
# Pkg.add("PyPlot")
# Pkg.add("FITSIO")
# Pkg.add("Images")
# Pkg.add("Wavelets")
# Pkg.add("HDF5")

################################
