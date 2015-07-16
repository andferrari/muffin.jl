using HDF5, JLD
using PyPlot

file = "/Users/deguignet/Documents/Julia/imagerec.jld"

admmst = load(file,"admmst")
x = admmst.x
lastiter = load(file, "algost.lastiter")
nfreq = load(file, "algost.nfreq")

tol1 = load(file, "toolst.tol1")
tol2 = load(file, "toolst.tol2")
tol3 = load(file, "toolst.tol3")
tol4 = load(file, "toolst.tol4")
tol5 = load(file, "toolst.tol5")
mydata = load(file, "skyst.mydata")
sky = load(file, "skyst.sky")

nu = load(file, "psfst.nu")
nu0 = load(file, "psfst.nu0")


d = size(x)[1]
if isodd(size(x)[1]) == true
    d12 = round((size(x)[1])/2)
else d12 = d/2
end

#############################################################################
#############################################################################

include("alpha_rec.jl")

#############################################################################
#############################################################################
######################################   Objet reconstruit
# figure(1)
# clf()
# title("Reconstructed object")
# e = 0
# for z in [1 5 10 15]
#     e += 1
#     subplot(2,2,e)
#     a = nu[z]/1e9
#     axis("off")
#     imshow(x[:,:,z],vmin=0,vmax=0.9)
#     title("v = $a[z] GHz")
# end
# subplots_adjust(bottom=0.1, right=0.8, top=0.95)
# cax = axes([0.85, 0.1, 0.025, 0.8])
# colorbar(cax=cax)
#############################################################################
#############################################################################
######################################  RMSE
# figure(2)
# clf()
# errr = zeros(Float64,15)
# for z = 1:nfreq
#     errr[z] = sqrt(sum((x[:,:,z] - sky[:,:,z]).^2)/sum(sky[:,:,z].^2))
# end
# # plot(nu./1e9,err[350,:]',marker="o",linewidth=2)
# plot(nu./1e9,errr[:],marker="o",linewidth=2)
# xlabel(L"\nu \; (GHz)",fontsize=18)
# ylabel(L"Relative \, rmse",fontsize=16)
# axis([1, 1.8, 0.1, 0.2])
#############################################################################
#############################################################################
######################################  cste / alpha / beta
# figure(3)
# clf()
# title("Cste Reconstruction")
# subplot(121)
# colorbar(imshow(-alpharec[:,:,1]))
# subplot(122)
# colorbar(imshow(-alpharecsky[:,:,1]))
#
# figure(4)
# clf()
# title("Alpha Reconstruction")
# subplot(121)
# colorbar(imshow(-alpharec[:,:,2]))
# subplot(122)
# colorbar(imshow(-alpharecsky[:,:,2]))
#
# figure(5)
# clf()
# title("Beta Reconstruction")
# subplot(121)
# colorbar(imshow(-alpharec[:,:,3]))
# subplot(122)
# colorbar(imshow(-alpharecsky[:,:,3]))
############################################# plot alpha
#############################################
figure(6)
clf()
m = d12-40
n = d12+40
p = linspace(m,n,int64(n-m+1))
q = alpharecsky[d12,m:n,2]'
r = alpharec[d12,m:n,2]'
dirty = alpharecdata[d12,m:n,2]'
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
p = linspace(d12-40,d12+40,81)
q = alpharecsky[d12,d12-40:d12+40,3]'
r = alpharec[d12,d12-40:d12+40,3]'
dirty = alpharecdata[d12,d12-40:d12+40,3]'
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
    i = d12
    j = d12-20

    q = 1*(squeeze(squeeze(x[i,j,:],1),1))
    d = (squeeze(squeeze(sky[i,j,:],1),1))

    plot(log10(p),log10(d),linewidth=2)
    plot(log10(p),log10(q),linewidth=2,linewidth=2)

    label1 = "Sky spectrum"
    label2 = "Reconstructed spectrum"

    legend( (label1, label2), loc="upper right")
    # xlabel(L"$log(\nu/\nu0)$,fontsize=16)
    # xlabel(L"$log(\nu/\nu0)$",fontsize=16)

subplot(1,2,2)
    p = nu./nu0
    i = d12
    j = d12+20

    q = 1.*(squeeze(squeeze(x[i,j,:],1),1))
    d = (squeeze(squeeze(sky[i,j,:],1),1))

    plot(log10(p),log10(d),linewidth=2)
    plot(log10(p),log10(q),linewidth=2)

    # plot(log10(p),log10(s))

    label1 = "Sky spectrum"
    label2 = "Reconstructed spectrum"


    legend( (label1, label2), loc="upper left")
    # xlabel(L"$log(\nu/\nu0)$,fontsize=16)
    # xlabel(L"$log(\nu/\nu0)$",fontsize=16)
#############################################
#############################################

#############################################################################
#############################################################################
######################################  sky / dirty / rec_obj
# figure(9)
# clf()
# title("multiplot")
# e = 1
# for z in [1 5 10 15]
#     subplot(3,4,e)
#     axis("off")
#     a = nu[z]/1e9
#     imshow(sky[:,:,z],vmin=0,vmax=1.3)
#     title(string(L"$\nu \, = \,$","$a ", L"\, GHz"),fontsize=14)
#     e += 1
# end
#
# e = 5
# for z in [1 5 10 15]
#     subplot(3,4,e)
#     axis("off")
#     imshow(mydata[:,:,z]./maximum(mydata).*1.2,vmin=0,vmax=1.3)
#
#     e += 1
# end
# e = 9
# for z in [1 5 10 15]
#     subplot(3,4,e)
#     axis("off")
#     imshow(x[:,:,z],vmin=0,vmax=1.3)
#
#     e += 1
# end
#
# subplots_adjust(bottom=0.1, right=0.8, top=0.9)
# cax = axes([0.85, 0.1, 0.025, 0.8])
# colorbar(cax=cax)
#############################################################################
#############################################################################
######################################  Reconstruction error
# figure(10)
# clf()
# e = 1
# for z in [1 5 10 15]
#     subplot(2,2,e)
#     axis("off")
#     a = nu[z]/1e9
#     # imshow(sqrt(abs(sky[:,:,z]-x[:,:,z])),vmin = 0, vmax = round(maximum(sqrt(abs(sky[:,:,:]-x[:,:,:]))),1))
#     imshow(((sky[:,:,z]-x[:,:,z])),vmin = 0, vmax = 0.1)
#     title(string(L"$\nu \, = \,$","$a ", L"\, GHz"))
#     e += 1
# end
#
# subplots_adjust(bottom=0.1, right=0.8, top=0.9)
# cax = axes([0.85, 0.1, 0.025, 0.8])
# colorbar(cax=cax)
#############################################################################
#############################################################################
######################################  sky / dirty / rec_obj / error
# figure(11)
# clf()
# title("multiplot")
# e = 1
# for z in [1 5 10 15]
#     subplot(4,4,e)
#     axis("off")
#     a = nu[z]/1e9
#     imshow(sky[:,:,z],vmin=0,vmax=1.3)
#     title(string(L"$\nu \, = \,$","$a ", L"\, GHz"),fontsize=14)
#     e += 1
# end
# e = 5
# for z in [1 5 10 15]
#     subplot(4,4,e)
#     axis("off")
#     imshow(mydata[:,:,z]./maximum(mydata).*1.2,vmin=0,vmax=1.3)
#     e += 1
# end
# e = 9
# for z in [1 5 10 15]
#     subplot(4,4,e)
#     axis("off")
#     imshow(x[:,:,z],vmin=0,vmax=1.3)
#     e += 1
# end
#
# subplots_adjust(bottom=0.1, right=0.8, top=0.9)
# cax = axes([0.85, 0.3, 0.025, 0.6])
# colorbar(cax=cax)
#
#
# e = 13
# for z in [1 5 10 15]
#     subplot(4,4,e)
#     axis("off")
#     imshow(((sky[:,:,z]-x[:,:,z])),vmin = 0, vmax = 0.08)
#     e += 1
# end
#
# subplots_adjust(bottom=0.1, right=0.8, top=0.9)
# cax = axes([0.85, 0.1, 0.025, 0.17])
# colorbar(cax=cax)
# #############################################################################
# #############################################################################
# ######################################  sky / dirty / error
# figure(12)
# clf()
# title("multiplot")
# e = 1
# for z in [1 5 10 15]
#     subplot(3,4,e)
#     axis("off")
#     a = nu[z]/1e9
#     imshow(sky[:,:,z],vmin=0,vmax=1.3)
#     title(string(L"$\nu \, = \,$","$a ", L"\, GHz"),fontsize=14)
#     e += 1
# end
# e = 5
# for z in [1 5 10 15]
#     subplot(3,4,e)
#     axis("off")
#     imshow(mydata[:,:,z]./maximum(mydata).*1.2,vmin=0,vmax=1.3)
#     e += 1
# end
#
# subplots_adjust(bottom=0.1, right=0.8, top=0.9)
# cax = axes([0.85, 0.41, 0.025, 0.46])
# colorbar(cax=cax)
#
# e = 9
# for z in [1 5 10 15]
#     subplot(3,4,e)
#     axis("off")
#     imshow(sqrt(abs(sky[:,:,z]-x[:,:,z])),vmin = 0, vmax = round(maximum(sqrt(abs(sky[:,:,:]-x[:,:,:]))),1))
#     e += 1
# end
# #vmax = round(maximum(sky[:,:,:]-x[:,:,:]),1)
# subplots_adjust(bottom=0.1, right=0.8, top=0.9)
# cax = axes([0.85, 0.13, 0.025, 0.18])
# colorbar(cax=cax)
#############################################################################
#############################################################################
