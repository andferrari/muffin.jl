####################################################################
########################## Initialisation ##########################
####################################################################
include("func.jl")
include("prox.jl")
include("testobj.jl")

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

const fov = 45
const arcminrad = 2*pi*60/360
const k = fov*arcminrad


# load psf fits file created by meqtrees
psf = "../data/meerkat_m30_25pix.psf.fits"
psfcube = lecture(psf)
psfavg = cubeaverage(psfcube,5)
mypsf = cropcubexy(psfavg,255)
# for z = 1:nw
#     mypsf[:,:,z] = mypsf[:,:,z]/abs(sum(mypsf[:,:,z]))
# end
mypsfadj = flipdim(flipdim(mypsf,1),2)

# load gray sky model fits file
# obj = "../data/M31.fits"
# sky0 = lecture(obj)/maximum(lecture(obj))
# sky,alpha = sky2cube(sky0,nu)
# noise = randn(size(sky)[1],size(sky)[1],size(mypsf)[3])/k
# mydata = cubefilter(sky,mypsf) + 10*noise

objdum = zeros(Float64,256,1)
sky = createobj(objdum)
noise = randn(size(sky)[1],size(sky)[1],size(mypsf)[3])#/k
mydata = cubefilter(sky,mypsf) #+ noise


spatialwlt  = [WT.db1,WT.db2,WT.db3,WT.db4,WT.db5,WT.db6,WT.db7,WT.db8,WT.haar]

const nspat = length(spatialwlt)
const nfreq = size(mypsf)[3]
#nfreq = 1
const nspec = 1
const nxy = size(mydata)[1]

niter = 0
lastiter = 0
const nbitermax = 350

const rhop = 1
const rhot = 5
const rhov = 2
const rhos = 1
const μt = 1e-1
const μv = 5e-1
const muesp = 0.1
const tt = rhot*nspat
const mu = muesp + rhop + tt + rhos


spectralwlt = zeros(Float64,nxy,nxy,nfreq)

snr = Float64[]
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
snr0 = 10*log10(vecnorm(cubefilter(sky,mypsf))^2/(vecnorm(noise)^2))
fty = cubefilter(mydata,mypsfadj)
push!(snr,10*log10(vecnorm(sky[:,:,1])^2/vecnorm(sky[:,:,1]-mydata[:,:,1])^2))
x[:] = mydata
