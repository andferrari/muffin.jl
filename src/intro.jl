using FITSIO
@everywhere using Images
using Wavelets
using GHF

include("structure.jl")
include("tmp.jl")
include("func.jl")

nw = 15
nu = zeros(Float64,nw)
for i = 1:nw
    nu[i] = 1.025e9 + (i-1)*50e6
end
nu0 = (nu[end]+nu[1])/2

file = FITS("")
header = readheader(file[1])
nustart = header["CRVAL4"]
nustep = header["CDELT4"]
nu0 = header["RESTFRQ"]

##################################
#### Structure initialisation ####
##################################



##################################
psf = "../data/meerkat_m30_25pix.psf.fits"
obj = "../data/M31.fits"
##################################
psfst = loadpsf(psf)
skyst = loadsky(obj,nu)
##################################
spatialwlt  = [WT.db1,WT.db2,WT.db3,WT.db4,WT.db5,WT.db6,WT.db7,WT.db8,WT.haar]
const nspat = length(spatialwlt)
const nfreq = size(psfst.mypsf)[3]
const nspec = 1
const nxy = size(skyst.mydata)[1]
niter = 0
lastiter = 0
const nitermax = 2000
##################################
algost = loadparam(nspat,nfreq,nspec,nxy,niter,lastiter,nitermax)
admmst = loadarray()
toolst = loadtools()
##################################
