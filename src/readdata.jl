using HDF5, JLD
using PyPlot
include("testobj.jl");

file = string("../data/results/data_2g.jld")
# m31 => data_m31_muesp0001_20db.jld
# 2g  => data_2gauss_muesp0001.jld

x = load(file, "x")
lastiter = load(file, "lastiter")
nfreq = load(file, "nfreq")
errorrec = load(file, "errorrec")
errorest = load(file, "errorest")
errorraw = load(file, "errorraw")
err = load(file, "err")
tol1 = load(file, "tol1")
tol2 = load(file, "tol2")
tol3 = load(file, "tol3")
tol4 = load(file, "tol4")
tol5 = load(file, "tol5")
mydata = load(file, "mydata")
sky = load(file, "sky")
#snr = load(file, "snr")
#nu = load(file, "nu")
nu0 = load(file, "nu0")
spectrex = load(file, "spectrex")
spectresky = load(file, "spectresky")
#noise = load("../file", "noise")

#############################
nw = 15
nu = zeros(Float64,nw)

for i = 1:nw
    nu[i] = 1.025e9 + (i-1)*50e6
end
nxy = size(x)[1]
##############################

#############################################################################
#############################################################################

include("alpha_rec.jl")

#############################################################################
#############################################################################
######################################   Objet reconstruit
figure(1)
clf()
title("Reconstructed object")
e = 0
for z in [1 5 10 15]
    e += 1
    subplot(2,2,e)
    a = nu[z]/1e9
    axis("off")
    imshow(x[:,:,z],vmin=0,vmax=0.9)
    title("v = $a[z] GHz")
end
subplots_adjust(bottom=0.1, right=0.8, top=0.95)
cax = axes([0.85, 0.1, 0.025, 0.8])
colorbar(cax=cax)
#############################################################################
#############################################################################
######################################  RMSE
figure(2)
clf()
errr = zeros(Float64,15)
for z = 1:nfreq
    errr[z] = sqrt(sum((x[:,:,z] - sky[:,:,z]).^2)/sum(sky[:,:,z].^2))
end
#plot(nu./1e9,err[350,:]',marker="o",linewidth=2)
plot(nu./1e9,errr[:],marker="o",linewidth=2)
xlabel(L"\nu \; (GHz)",fontsize=18)
ylabel(L"Relative \, rmse",fontsize=16)
axis([1, 1.8, 0.15, 0.22])
#############################################################################
#############################################################################
######################################  cste / alpha / beta
figure(3)
clf()
title("Cste Reconstruction")
subplot(121)
colorbar(imshow(-alpharec[:,:,1]))
subplot(122)
colorbar(imshow(-alpharecsky[:,:,1]))

figure(4)
clf()
title("Alpha Reconstruction")
subplot(121)
colorbar(imshow(-alpharec[:,:,2]))
subplot(122)
colorbar(imshow(-alpharecsky[:,:,2]))

figure(5)
clf()
title("Beta Reconstruction")
subplot(121)
colorbar(imshow(-alpharec[:,:,3]))
subplot(122)
colorbar(imshow(-alpharecsky[:,:,3]))
############################################# plot alpha
#############################################
figure(6)
clf()
m = 88
n = 168
p = linspace(m,n,n-m+1)
q = alpharecsky[128,m:n,2]'
r = alpharec[128,m:n,2]'
dirty = alpharecdata[128,m:n,2]'
plot(p,q,linewidth=2)
plot(p,r,linewidth=2)
plot(p,dirty,linewidth=2)
legend( ("Sky","Reconstructed sky", "Dirty image"), loc="upper left")
xlabel(L"pixel",fontsize=18)
ylabel(L"Spectral \, index \, \alpha",fontsize=18)
############################################# plot beta
#############################################
figure(7)
clf()
p = linspace(88,168,81)
q = alpharecsky[128,88:168,3]'
r = alpharec[128,88:168,3]'
dirty = alpharecdata[128,88:168,3]'
plot(p,q)
plot(p,r)
plot(p,dirty)
legend( ("Sky","Reconstructed sky", "Dirty image"), loc="lower center")
xlabel("pix")
ylabel(L"Spectral index $\beta$")
############################################# plot spectre
#############################################
figure(8)
clf()

subplot(1,2,1)
    p = nu./nu0
    i = 128
    j = 108

    q = 1*(squeeze(squeeze(x[i,j,:],1),1))
    d = (squeeze(squeeze(sky[i,j,:],1),1))

    plot(log10(p),log10(d),linewidth=2)
    plot(log10(p),log10(q),linewidth=2,linewidth=2)

    label1 = "Sky spectrum"
    label2 = "Reconstructed spectrum"

    legend( (label1, label2), loc="upper right")
    xlabel(L"$log(\nu/\nu0)$",fontsize=16)

subplot(1,2,2)
    p = nu./nu0
    i = 128
    j = 148

    q = 1.*(squeeze(squeeze(x[i,j,:],1),1))
    d = (squeeze(squeeze(sky[i,j,:],1),1))

    plot(log10(p),log10(d),linewidth=2)
    plot(log10(p),log10(q),linewidth=2)

    #plot(log10(p),log10(s))

    label1 = "Sky spectrum"
    label2 = "Reconstructed spectrum"


    legend( (label1, label2), loc="upper left")
    xlabel(L"$log(\nu/\nu0)$",fontsize=16)
#############################################
#############################################

#############################################################################
#############################################################################
######################################  sky / dirty / rec_obj
figure(9)
clf()
title("multiplot")
e = 1
for z in [1 5 10 15]
    subplot(3,4,e)
    axis("off")
    a = nu[z]/1e9
    imshow(sky[:,:,z],vmin=0,vmax=1.3)
    title(string(L"$\nu \, = \,$","$a ", L"\, GHz"),fontsize=14)
    e += 1
end

e = 5
for z in [1 5 10 15]
    subplot(3,4,e)
    axis("off")
    imshow(mydata[:,:,z]./maximum(mydata).*1.2,vmin=0,vmax=1.3)

    e += 1
end
e = 9
for z in [1 5 10 15]
    subplot(3,4,e)
    axis("off")
    imshow(x[:,:,z],vmin=0,vmax=1.3)

    e += 1
end

subplots_adjust(bottom=0.1, right=0.8, top=0.9)
cax = axes([0.85, 0.1, 0.025, 0.8])
colorbar(cax=cax)
#############################################################################
#############################################################################
######################################  Reconstruction error
figure(10)
clf()
e = 1
for z in [1 5 10 15]
    subplot(2,2,e)
    axis("off")
    a = nu[z]/1e9
    #imshow(sqrt(abs(sky[:,:,z]-x[:,:,z])),vmin = 0, vmax = round(maximum(sqrt(abs(sky[:,:,:]-x[:,:,:]))),1))
    imshow(((sky[:,:,z]-x[:,:,z])),vmin = 0, vmax = 0.1)
    title(string(L"$\nu \, = \,$","$a ", L"\, GHz"))
    e += 1
end

subplots_adjust(bottom=0.1, right=0.8, top=0.9)
cax = axes([0.85, 0.1, 0.025, 0.8])
colorbar(cax=cax)
#############################################################################
#############################################################################
######################################  sky / dirty / rec_obj / error
figure(11)
clf()
title("multiplot")
e = 1
for z in [1 5 10 15]
    subplot(4,4,e)
    axis("off")
    a = nu[z]/1e9
    imshow(sky[:,:,z],vmin=0,vmax=1.3)
    title(string(L"$\nu \, = \,$","$a ", L"\, GHz"),fontsize=14)
    e += 1
end
e = 5
for z in [1 5 10 15]
    subplot(4,4,e)
    axis("off")
    imshow(mydata[:,:,z]./maximum(mydata).*1.2,vmin=0,vmax=1.3)
    e += 1
end
e = 9
for z in [1 5 10 15]
    subplot(4,4,e)
    axis("off")
    imshow(x[:,:,z],vmin=0,vmax=1.3)
    e += 1
end

subplots_adjust(bottom=0.1, right=0.8, top=0.9)
cax = axes([0.85, 0.3, 0.025, 0.6])
colorbar(cax=cax)


e = 13
for z in [1 5 10 15]
    subplot(4,4,e)
    axis("off")
    imshow(((sky[:,:,z]-x[:,:,z])),vmin = 0, vmax = 0.08)
    e += 1
end

subplots_adjust(bottom=0.1, right=0.8, top=0.9)
cax = axes([0.85, 0.1, 0.025, 0.17])
colorbar(cax=cax)
#############################################################################
#############################################################################
######################################  sky / dirty / error
figure(12)
clf()
title("multiplot")
e = 1
for z in [1 5 10 15]
    subplot(3,4,e)
    axis("off")
    a = nu[z]/1e9
    imshow(sky[:,:,z],vmin=0,vmax=1.3)
    title(string(L"$\nu \, = \,$","$a ", L"\, GHz"),fontsize=14)
    e += 1
end
e = 5
for z in [1 5 10 15]
    subplot(3,4,e)
    axis("off")
    imshow(mydata[:,:,z]./maximum(mydata).*1.2,vmin=0,vmax=1.3)
    e += 1
end

subplots_adjust(bottom=0.1, right=0.8, top=0.9)
cax = axes([0.85, 0.41, 0.025, 0.46])
colorbar(cax=cax)

e = 9
for z in [1 5 10 15]
    subplot(3,4,e)
    axis("off")
    imshow(sqrt(abs(sky[:,:,z]-x[:,:,z])),vmin = 0, vmax = round(maximum(sqrt(abs(sky[:,:,:]-x[:,:,:]))),1))
    e += 1
end
#vmax = round(maximum(sky[:,:,:]-x[:,:,:]),1)
subplots_adjust(bottom=0.1, right=0.8, top=0.9)
cax = axes([0.85, 0.13, 0.025, 0.18])
colorbar(cax=cax)
#############################################################################
#############################################################################
