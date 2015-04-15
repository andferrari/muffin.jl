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
#sky = load("../data/results/data.jld", "sky")
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
clf()
title("Reconstructed object")
e = 0
for z in [1 5 10 15]
    e += 1
    subplot(2,2,e)
    a = nu[z]/1e9
    axis("off")
    imshow(mypsf[108:148,108:148,z])
    title("v = $a GHz")
end
subplots_adjust(bottom=0.1, right=0.8, top=0.9)
cax = axes([0.85, 0.1, 0.025, 0.8])
colorbar(cax=cax)
# savefig("psf.pdf")

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


# figure(10)
# clf()
# title("multiplot")
# e = 1
# for z in [1 5 10 15]
#     subplot(4,3,e)
#     a = nu[z]/1e9
#     axis("off")
#     imshow(sky[:,:,z])
#     title("v = $a GHz")
#     e += 3
# end
# e = 2
# for z in [1 5 10 15]
#     subplot(4,3,e)
#     a = nu[z]/1e9
#     axis("off")
#     imshow(mydata[:,:,z])
#     title("v = $a GHz")
#     e += 3
# end
# e = 3
# for z in [1 5 10 15]
#     subplot(4,3,e)
#     a = nu[z]/1e9
#     axis("off")
#     imshow(x[:,:,z])
#     title("v = $a GHz")
#     e += 3
# end
#
# subplots_adjust(bottom=0.1, right=0.8, top=0.9)
# cax = axes([0.85, 0.1, 0.025, 0.8])
# colorbar(cax=cax)


figure(10)
clf()
title("multiplot")
e = 1
for z in [1 5 10 15]
    subplot(3,4,e)
    a = nu[z]/1e9
    axis("off")
    imshow(sky[:,:,z])
    title("v = $a GHz")
    e += 1
end

e = 5
for z in [1 5 10 15]
    subplot(3,4,e)
    a = nu[z]/1e9
    axis("off")
    imshow(mydata[:,:,z])

    e += 1
end
e = 9
for z in [1 5 10 15]
    subplot(3,4,e)
    a = nu[z]/1e9
    axis("off")
    imshow(x[:,:,z])

    e += 1
end

subplots_adjust(bottom=0.1, right=0.8, top=0.9)
cax = axes([0.85, 0.1, 0.025, 0.8])
colorbar(cax=cax)
