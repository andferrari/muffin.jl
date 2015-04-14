using HDF5, JLD
using PyPlot
include("testobj.jl");

x = load("../data/results/data.jld", "x")
lastiter = load("../data/results/data.jld", "lastiter")
nfreq = load("../data/results/data.jld", "nfreq")
errorrec = load("../data/results/data.jld", "errorrec")
errorest = load("../data/results/data.jld", "errorest")
errorraw = load("../data/results/data.jld", "errorraw")
err = load("../data/results/data.jld", "err")
tol1 = load("../data/results/data.jld", "tol1")
tol2 = load("../data/results/data.jld", "tol2")
tol3 = load("../data/results/data.jld", "tol3")
tol4 = load("../data/results/data.jld", "tol4")
tol5 = load("../data/results/data.jld", "tol5")
mydata = load("../data/results/data.jld", "mydata")
sky = load("../data/results/data.jld", "sky")
#snr = load("../data/results/data.jld", "snr")
# nu = load("../data/results/data.jld", "nu")
# nu0 = load("../data/results/data.jld", "nu0")

#############################
nw = 15
nu = zeros(Float64,nw)
for i = 1:nw
    nu[i] = 1.025e9 + (i-1)*50e6
end
nxy = size(x)[1]
##############################

figure(1)
for z = 1:nfreq
    subplot(5,3,z)
    colorbar(imshow(x[:,:,z]))
end

figure(2)
for z = 1:nfreq
    subplot(5,3,z)
    colorbar(imshow(errorrec[:,:,z]))
end

figure(3)
plot([1:nfreq],errorest,color="blue")
#plot([1:nfreq],errorraw,color="red")

# figure(4)
# for niter = 1:5
#
#     clf()
#     println(10*niter)
#     subplot(5,1,1)
#     plot(tol1[1:10*niter])
#     subplot(5,1,2)
#     plot(tol2[1:10*niter])
#     subplot(5,1,3)
#     plot(tol3[1:10*niter])
#     subplot(5,1,4)
#     plot(tol4[1:10*niter])
#     subplot(5,1,5)
#     plot(tol5[1:10*niter])
# end
#
# figure(5)
# for niter = 1:2
#     println(25*niter)
#     clf()
#     for z = 1:nfreq
#         subplot(5,3,z)
#         plot(err[1:25*niter,z])
#     end
# end

include("alpha_rec.jl")
