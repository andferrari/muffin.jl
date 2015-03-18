# initialisation

using FITSIO

# write cdf file for paraview visualization
module Writecdf
export writecdf

using NetCDF
    function writecdf{T<:Float64}(filename::String,psfcube,datacube::Array{T,3})
        nx, ny, nfreq = size(datacube)
        nccreate(filename,"psf","x",nx,"y",ny,"frequency",nfreq)
        ncwrite(psfcube,filename,"psf")
        ncclose()
    end
end

# load psf fits file created by meqtrees
file = FITS("../data/meerkat_psf_32pix_100ch.fits")
data = read(file[1])
close(file)
psfcube = squeeze(data,3)
nxpsf, nypsf, nfreq = size(psfcube)

file = FITS("../data/cluster.fits")
data = read(file[1])
cluster = squeeze(squeeze(data,4),3)
nxmod, nymod = size(cluster)

data = Array(Float64,nxmod,nymod,nfreq)
for k = 1:nfreq
    data[:,:,k] = conv2(cluster,psfcube[:,:,k])
end
