using HDF5, JLD
using PyPlot
using FITSIO


###############################
########## Rec Image ##########
file = string("/Users/deguignet/Documents/Julia/imagerec_400ite_woshar.jld")
x = load(file, "admmst.x");
x = x[1024-512:1024+512-1,1024-512:1024+512-1,:];

file2 = string("/Users/deguignet/Documents/Julia/imagerec_700ite.jld")
x2 = load(file2, "admmst.x");
x2 = x2[1024-512:1024+512-1,1024-512:1024+512-1,:];

###############################
####### Original object #######
skyfile = string("/Users/deguignet/Documents/Julia/I_HALO_CL_RS_SKY.FITS")
f = FITS(skyfile)
sky = read(f[1]);
close(f)

###############################
############ Dirty ############
dirtyfile = string("/Users/deguignet/Documents/Julia/example_sim_dirty.fits")
f = FITS(dirtyfile)
mydata = read(f[1]);
close(f)

###############################
############# PSF #############
psffile = string("/Users/deguignet/Documents/Julia/example_sim_psf.fits")
f = FITS(psffile)
psf = read(f[1]);
close(f)

###############################
nxy  = size(x)[1]
lastiter = load(file, "algost.lastiter")
nfreq = load(file, "algost.nfreq")
tol1 = load(file, "toolst.tol1")
tol2 = load(file, "toolst.tol2")
tol3 = load(file, "toolst.tol3")
tol4 = load(file, "toolst.tol4")
tol5 = load(file, "toolst.tol5")
nu = load(file, "psfst.nu")
nu0 = load(file, "psfst.nu0")

sumsky2 = Array(Float64,nfreq)
for z in 1:nfreq
    sumsky2[z] = sum(sky[:,:,z].^2)
end
###############################

###############################
############ RMSE #############
# rmse = Array(Float64,nfreq)
# for z in 1:nfreq
#     rmse[z] = sqrt(sum( (x[:,:,z] - sky[:,:,z]).^2)/sumsky2[z])
# end
#
#
# figure(1)
# clf()
# plot(nu,rmse)
##########################
#
# err = Array(Float64,nxy,nxy,nfreq)
# for z in 1:nfreq
#     err[:,:,z] = (sky[:,:,z] - x[:,:,z]).^2./(sky[:,:,z].^2)
# end


##########################


##########################
#
# figure(2)
# clf()
# for z in 1:nfreq
#     subplot(6,6,z)
#     a = nu[z]/1e9
#     imshow(x[:,:,z],cmap=ColorMap("ocean"))
#     axis("off")
#     title("v = $a GHz")
#
# end

#############################################################################
#############################################################################
######################################  sky / dirty / rec_obj / error
# figure(3)
# clf()
# title("multiplot")
# e = 1
# for z in [8 16 24 32]
#     subplot(4,4,e)
#     axis("off")
#     a = nu[z]/1e9
#     imshow((abs(sky[:,:,z]).^0.25),cmap=ColorMap("spectral"))
#     title(string(L"$\nu \, = \,$","$a ", L"\, GHz"),fontsize=14)
#     e += 1
# end
# subplots_adjust(bottom=0.1, right=0.8, top=0.9)
# cax = axes([0.85, 0.3, 0.025, 0.6])
# colorbar(cax=cax)
# e = 5
# for z in [8 16 24 32]
#     subplot(4,4,e)
#     axis("off")
#     imshow((mydata[:,:,z]),cmap=ColorMap("spectral"))
#     e += 1
# end
# subplots_adjust(bottom=0.1, right=0.8, top=0.9)
# cax = axes([0.85, 0.3, 0.025, 0.6])
# colorbar(cax=cax)
#
# e = 9
# for z in [8 16 24 32]
#     subplot(4,4,e)
#     axis("off")
#     imshow((abs(x[:,:,z])/maximum(x[:,:,z])*maximum(sky[:,:,z])).^0.5,cmap=ColorMap("spectral"))
#     e += 1
# end
#
# subplots_adjust(bottom=0.1, right=0.8, top=0.9)
# cax = axes([0.85, 0.3, 0.025, 0.6])
# colorbar(cax=cax)
#
#
# e = 13
# for z in [8 16 24 32]
#     subplot(4,4,e)
#     axis("off")
#     imshow((abs(sky[:,:,z]-x[:,:,z])/(sum(sky[:,:,z]))).^0.25,cmap=ColorMap("spectral"))
#     e += 1
# end
#
# subplots_adjust(bottom=0.1, right=0.8, top=0.9)
# cax = axes([0.85, 0.1, 0.025, 0.17])
# colorbar(cax=cax)
#############################################################################
#############################################################################

# ################# PSF
# figure(4)
# clf()
# e=1
# for z in [8 16 24 32]
#
#     subplot(2,2,e)
#     axis("off")
#     a = nu[z]/1e9
#     imshow((abs(psf[:,:,z]).^0.25),cmap=ColorMap("spectral"))
#         title(string(L"$\nu \, = \,$","$a ", L"\, GHz"),fontsize=14)
#     e+=1
# end
#
# ################# Object
# figure(5)
# clf()
# e=1
# for z in [8 16 24 32]
#
#     subplot(2,2,e)
#     axis("off")
#     a = nu[z]/1e9
#     imshow((abs(sky[:,:,z]).^0.25),cmap=ColorMap("spectral"))
#         title(string(L"$\nu \, = \,$","$a ", L"\, GHz"),fontsize=14)
#     e+=1
# end


##################################
##################################
figure(6)
clf()
subplot(2,3,1)
    axis("off")
    colorbar(imshow((abs(x[:,:,1]).^0.25),cmap=ColorMap("spectral")))
    title("400 iterations w/o shaar")
    ylabel("puissance 1/4")
subplot(2,3,2)
    axis("off")
    colorbar(imshow((abs(x2[:,:,1]).^0.25),cmap=ColorMap("spectral")))
    title("700 iterations")
subplot(2,3,3)
    axis("off")
    colorbar(imshow((abs(sky[:,:,1]).^0.25),cmap=ColorMap("spectral")))
    title("original object")

    subplot(2,3,4)
        axis("off")
        colorbar(imshow((abs(x[:,:,1])),cmap=ColorMap("spectral")))
        title("400 iterations w/o shaar")
        ylabel("puissance 1")

    subplot(2,3,5)
        axis("off")
        colorbar(imshow((abs(x2[:,:,1])),cmap=ColorMap("spectral")))
        title("700 iterations")
    subplot(2,3,6)
        axis("off")
        colorbar(imshow((abs(sky[:,:,1])),cmap=ColorMap("spectral")))
        title("original object")


# subplot(1,3,1)
#     axis("off")
#     imshow(((x[:,:,1]).^1),cmap=ColorMap("spectral"))
# subplot(1,3,2)
#     axis("off")
#     imshow(((x2[:,:,1]).^1),cmap=ColorMap("spectral"))
# subplot(1,3,3)
#     axis("off")
#     imshow(((sky[:,:,1]).^1),cmap=ColorMap("spectral"))


##################################
##################################
# figure(7)
# clf()
# xt = linspace(1,lastiter,lastiter)
# plot(xt,log(tol2))
# plot(xt,log(tol3))
# plot(xt,log(tol4))
# plot(xt,log(tol5))

##################################
##################################
figure(6)
clf()
subplot(1,2,1)
    axis("off")
    (imshow((abs(x[:,:,1]).^0.25),cmap=ColorMap("spectral"),vmax=0.1))

subplot(1,2,2)
    axis("off")
    (imshow((abs(sky[:,:,1]).^0.25),cmap=ColorMap("spectral"),vmax=0.12))

    # subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    subplots_adjust(bottom=0.1, right=0.8, top=1.1)
    cax = axes([0.85, 0.3, 0.025, 0.6])
    colorbar(cax=cax)
