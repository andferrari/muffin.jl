####################################################################
########################## Initialisation ##########################
####################################################################
include("func.jl")
include("prox.jl")

using FITSIO
@everywhere using Images
using Wavelets
using GHF

# freq

nw = 15
nu = zeros(Float64,nw)
for i = 1:nw
    nu[i] = 1.025e9 + (i-1)*50e6
end



# load psf fits file created by meqtrees

file = FITS("../data/meerkat_m30_25pix.psf.fits")
data = float64(read(file[1]))
close(file)
psfcube = squeeze(data,3)
psfavg = cubeaverage(psfcube,5)
mypsf = cropcubexy(psfavg,255)
mypsfadj = float64(flipdim(flipdim(mypsf,1),2))


# load gray sky model fits file

file = FITS("../data/M31.fits")
data = float64(read(file[1]))
close(file)
sky0 = squeeze(data,3)
sky,alpha = sky2cube(sky0,nu)
mydata = cubefilter(sky,mypsf)


spatialwlt  = [WT.db1,WT.db2,WT.db3,WT.db4,WT.db5,WT.db6,WT.db7,WT.db8,WT.haar]

const nspat = length(spatialwlt)
const nfreq = size(mypsf)[3]
const nspec = 1
const nxy = size(mydata)[1]

niter = 0
lastiter = 0
const nbitermax = 1000

const rhop = 2
const rhot = 1
const rhov = 5
const rhos = 5
const μt = 1.0
const μv = 1.0
const muesp = 1.0
const tt = rhot*nspat
const mu = muesp + rhop + tt + rhos

spectralwlt = zeros(Float64,nxy,nxy,nfreq)

tol1 = Float64[]
tol2 = Float64[]
tol3 = Float64[]
tol4 = Float64[]
tol5 = Float64[]
err = Array(Float64,nbitermax,nfreq)
errorrec = zeros(Float64,nxy,nxy,nfreq)
errorest = zeros(Float64,nfreq)
errorraw = zeros(Float64,nfreq)

s = zeros(Float64,nxy,nxy,nfreq)
taus = zeros(Float64,nxy,nxy,nfreq)
sh = zeros(Float64,nxy,nxy,nfreq)

taup = zeros(Float64,nxy,nxy,nfreq)
p = zeros(Float64,nxy,nxy,nfreq)

tauv = zeros(Float64,nxy,nxy,nfreq)
v = zeros(Float64,nxy,nxy,nfreq)

t = zeros(Float64,nxy,nxy,nfreq,nspat)
taut = zeros(Float64,nxy,nxy,nfreq,nspat)

wlt = SharedArray(Float64,nxy,nxy,nfreq)

x = SharedArray(Float64,nxy,nxy,nfreq)
Hx = SharedArray(Float64,nxy,nxy,nfreq,nspat)
xmm = zeros(Float64,nxy,nxy,nfreq)

# precompute

fty = cubefilter(mydata,mypsfadj)
