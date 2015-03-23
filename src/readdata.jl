using HDF5, JLD
using PyPlot

x = load("../data/results/data.jld", "x")
nbitermax = load("../data/results/data.jld", "nbitermax")
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

figure(1)
for z = 1:nfreq
    subplot(5,2,z)
    imshow(x[:,:,z])
end

figure(2)
for z = 1:nfreq
    subplot(5,2,z)
    imshow(errorrec[:,:,z])
end

figure(3)
plot([1:nfreq],errorest,color="blue")
plot([1:nfreq],errorraw,color="red")

figure(4)
for niter = 1:nbitermax
    println(10*niter)
    clf()
    for z = 1:nfreq
        subplot(5,2,z)
        plot(err[1:10*niter,z])
    end
end
