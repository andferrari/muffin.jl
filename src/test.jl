#using PyPlot
using NFFT
using Wavelets

include("test_tools.jl")

N = 64
Nb = 128
nw = 14

ind = 0
nbitermax = 100
const μ__tv = 1
const μ__ss = 1e-1
const μ_ϵ = 1E-9

const rho_y = 0.00001
const rho_T = 10
const rho_P = 1
const rho_S = 1
const rho_V = 1

FOV = 0.05

DegRad = 2.0*pi/360.0
RadArcsec = 3600.0/DegRad
coef = 1#FOV/RadArcsec

SpatialWaveletsBasisList  = ["db1","db2","db3","db4","db5","db6","db7","db8","haar"]
nspat = length(SpatialWaveletsBasisList)

Hx = zeros(Float64,N,N,nw,nspat)
sH = zeros(Float64,N,N,nw)
Fx = SharedArray(Complex{Float64},(Nb,nw))
SpcWvlt = zeros(Float64,N,N,nw)
# initialization Lagrange multipliers
tau_T = zeros(Float64,N,N,nw,nspat)
tau_P = zeros(Float64,N,N,nw)
tau_S = zeros(Float64,N,N,nw)
tau_V = zeros(Float64,N,N,nw)
tau_y = zeros(Complex128,Nb,nw)
tau_yc = zeros(Complex128,Nb,nw)

T = zeros(Float64,N,N,nw,nspat)
P = zeros(Float64,N,N,nw)
S = zeros(Float64,N,N,nw)
V = zeros(Float64,N,N,nw)

visi = Array(Complex128,Nb,nw)

yp = zeros(Complex128,Nb,nw)
y_est = zeros(Complex128,Nb,nw)
yc = zeros(Complex128,Nb,nw)


x = SharedArray(Float64,(N,N,nw))
xtt = Array(Complex128,N,N,nw)
x1 = Array(Complex128,N,N,nw)
x2 = Array(Float64,N,N,nw,nspat)

wvl = linspace(0.12,0.37,nw)

########################################
########### Creation objet ############
########################################
rd = 16
mid = N/2 +1

Obj = zeros(Float64,N,N,nw)
for i=1:N, j=1:N
 Obj[i,j,:] = exp(-((i-mid)^2+(j-mid)^2)/100)
end

pos = zeros(Float64,N,N,nw)

u = rand(Nb) - 0.5
v = rand(Nb) - 0.5

umat = Array(Float64,Nb,nw)
vmat = Array(Float64,Nb,nw)

for z = 1:nw, i = 1:Nb
 umat[i,z] = u[i]/wvl[z]/N
 vmat[i,z] = v[i]/wvl[z]/N
end

########################################

Mtmp = {((hcat(umat[:,n]*coef,vmat[:,n]*coef))',N) for n       = 1 : nw}

Plan = pmap(MyPlan,Mtmp)
const F3D  = nudft3d(umat*coef,vmat*coef)
K=1
matinv = Array(Complex128,Nb,Nb,nw)
for z = 1:nw
 visi[:,z] = nfft(Plan[z],Obj[:,:,z]+0.*im)
 matinv[:,:,z] = inv(speye(Nb)+rho_y/K*F3D[z]*F3D[z]')
end

include("test_loop.jl")
