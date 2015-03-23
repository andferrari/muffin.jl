using HDF5, JLD
using PyPlot

x = load("data.jld", "x")
nbitermax = load("data.jld", "nbitermax")
nfreq = load("data.jld", "nfreq")
errorrec = load("data.jld", "errorrec")
errorest = load("data.jld", "errorest")
errorraw = load("data.jld", "errorraw")
err = load("data.jld", "err")
tol1 = load("data.jld", "tol1")
tol2 = load("data.jld", "tol2")
tol3 = load("data.jld", "tol3")
tol4 = load("data.jld", "tol4")
tol5 = load("data.jld", "tol5")

figure(1)
for z = 1:nfreq
    subplot(5,2,z)
    imshow(x[:,:,z])
end

figure(2)
for niter = 1:1000
    clf()
    for z = 1:nfreq
        subplot(5,2,z)
        plot(err[1:niter,z])
    end
end

figure(3)
for z = 1:nfreq
    subplot(5,2,z)
    imshow(errorrec[:,:,z])
end

figure(4)
plot([1:nfreq],errorest,color="blue")
plot([1:nfreq],errorraw,color="red")
