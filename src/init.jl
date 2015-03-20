# initialisation

using FITSIO
using PyPlot
@everywhere using Images

function writecdf{T<:FloatingPoint}(filename::ASCIIString,datacube::Array{T,3},paraview = false)
    if isfile(filename)
        error("file $filename already exists")
    end
    nx, ny, nfreq = size(datacube)
    nccreate(filename,"datacube","x",nx,"y",ny,"frequency",nfreq)
    ncwrite(datacube,filename,"datacube")
    ncclose()
    if paraview ==  true
        path = pwd()
        run(`open -a paraview --args --data=$path/$filename`)
    end
end


function cubefilter{T<:FloatingPoint}(imagecube::Array{T,3},psfcube::Array{T,3})
    nximg, nyimg, nfreq = size(imagecube)
    nxpsf, nypsf, nfreq = size(psfcube)
    rescube = Array(Float64,nximg,nyimg,nfreq)
    for k = 1:nfreq
        rescube[:,:,k] = imfilter_fft(imagecube[:,:,k],psfcube[:,:,k],"circular")
    end
    return rescube
end

function cubefilter{T<:FloatingPoint}(imagecube::Array{T,2},psfcube::Array{T,3})
    nximg, nyimg = size(imagecube)
    nxpsf, nypsf, nfreq = size(psfcube)
    rescube = Array(Float64,nximg,nyimg,nfreq)
    for k = 1:nfreq
        rescube[:,:,k] = imfilter_fft(imagecube,psfcube[:,:,k],"circular")
    end
    return rescube
end

# load psf fits file created by meqtrees
file = FITS("../data/meerkat_psf_33pix_100ch.fits")
#file = FITS("../data/meerkat_psf_17pix_100ch.fits")

data = read(file[1])
close(file)
psfcube = squeeze(data,3)
nxpsf, nypsf, nfreq = size(psfcube)

# load gray sky model fits file
file = FITS("../data/cluster.fits")
data = read(file[1])
cluster = squeeze(squeeze(data,4),3)

datacube = cubefilter(cluster,psfcube)
