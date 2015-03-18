# initialisation

using FITSIO
using Images

# module to write cdf file for paraview visualization
module Writecdf
export writecdf

using NetCDF
    function writecdf{T<:Float64}(filename::ASCIIString,datacube::Array{T,3},paraview = false)
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
end

# load psf fits file created by meqtrees
file = FITS("../data/meerkat_psf_32pix_100ch.fits")
data = read(file[1])
close(file)
psfcube = squeeze(data,3)
nxpsf, nypsf, nfreq = size(psfcube)

# load gray sky model fits file
file = FITS("../data/cluster.fits")
data = read(file[1])
cluster = squeeze(squeeze(data,4),3)
nxmod, nymod = size(cluster)

imgcube = Array(Float64,nxmod,nymod,nfreq)
for k = 1:nfreq
    imgcube[:,:,k] = imfilter_fft(cluster,psfcube[:,:,k])
end
