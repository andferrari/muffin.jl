include("func.jl")
include("prox.jl")

function minisky(sky,nu,alpha)
    nx, ny = size(sky)
    nbands = length(nu)
    skycube = Array(Float64, nx, ny, nbands)
    nu0 = (nu[end]+nu[1])/2

    for k =1:nbands
        skyc = sky.* ((nu[k]/nu0).^(-alpha))
        skycube[:,:,k] = skyc
        end
    return skycube
end

function createobj(toto;l=20,sig=400)

    nxy = size(toto)[1]
    obj = zeros(Float64,nxy,nxy,nw)

    aa = nxy/2-l
    bb = nxy/2+l

    gauss1 = zeros(Float64,nxy,nxy)
    gauss2 = zeros(Float64,nxy,nxy)

    for i = 1:nxy, j = 1:nxy
        gauss1[i,j] = exp(-((i-nxy/2)^2+(j-aa)^2)/sig)
        gauss2[i,j] = exp(-((i-nxy/2)^2+(j-bb)^2)/sig)
    end

    alpha1 = 1
    alpha2 = -1

    obj = minisky(gauss1,nu,alpha1) + minisky(gauss2,nu,alpha2)
    obj = prox_u(obj,0.01)

    return obj
end

# figure(1)
# clf();colorbar(imshow(obj[:,:,1]))

# x = obj

######################################
#include("alpha_rec.jl")

# nw = 15
# nu = zeros(Float64,nw)
# for i = 1:nw
#     nu[i] = 1.025e9 + (i-1)*50e6
# end
# nu0 = (nu[end]-nu[1])/2
#
# nxy = 256
# l = 20
# obj = zeros(Float64,nxy,nxy,nw)
#
# aa = nxy/2-l
# bb = nxy/2+l
# sig = 400
#
# gauss1 = zeros(Float64,nxy,nxy)
# gauss2 = zeros(Float64,nxy,nxy)
#
# for i = 1:nxy, j = 1:nxy
#     gauss1[i,j] = exp(-((i-nxy/2)^2+(j-aa)^2)/sig)
#     gauss2[i,j] = exp(-((i-nxy/2)^2+(j-bb)^2)/sig)
# end
#
# alpha1 = 1
# alpha2 = -1
#
#
#
# obj = minisky(gauss1,nu,alpha1) + minisky(gauss2,nu,alpha2)
# obj = prox_u(obj,0.01)
