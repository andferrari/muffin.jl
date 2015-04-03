# initialisation

using FITSIO
@everywhere using Images
using Wavelets
using GHF

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

function cubeaverage{T<:FloatingPoint}(imagecube::Array{T,3},M::Int)
    nxpsf, nypsf, nfreq = size(imagecube)
    if nfreq < M
        error("Not enough channels to average!")
    end
    nfreqavg = itrunc(nfreq/M)
    rescube = Array(Float64, nxpsf, nypsf, nfreqavg)
    for k = 1:nfreqavg
        rescube[:,:,k] = sum(imagecube[:,:,(k-1)*M+1:k*M], 3)/M
    end
    return rescube
end

function cropcubexy{T<:FloatingPoint}(imagecube::Array{T,3},M::Int)
    nxpsf, nypsf, nfreq = size(imagecube)
    if ~(nxpsf == nypsf)
        error("Image must be square!")
    end
    if nxpsf < M
        error("Image too small to be croped!")
    end

    rescube = Array(Float64, M, M, nfreq)
    if (iseven(nxpsf) && iseven(M)) | (isodd(nxpsf) && isodd(M))
        nstart = (nxpsf-M)/2 +1
        nend = (nxpsf+M)/2
        rescube = imagecube[nstart:nend,nstart:nend,:]
    else
        error("size of image and crop must have same parity!")
    end
    return rescube
end

function sky2cube{T<:FloatingPoint}(sky::Array{T,2},nu::Array{T,1})
    nx, ny = size(sky)
    nbands = length(nu)
    corl = minimum([nx,ny])/5.0
    rho(x,y) = exp(-norm([x, y])/corl)
    tx = linspace(0,nx-1,nx)
    ty = linspace(0,ny-1,ny)
    field = genghf(tx, ty, rho)
    skycube = Array(Float64, nx, ny, nbands)

    sm, sM = extrema(sky)
    fm, fM = extrema(field)

    alpha = 0.25 + 0.1(sky-sm)/(sM-sm) + 0.1(field-fm)/(fM-fm)

    nu0 = (nu[end]-nu[1])/2
    for k =1:nbands
        skycube[:,:,k] = sky.* (nu[k]/nu0) .^(-alpha)
    end
    return skycube
end

# freq
nw = 15
nu = zeros(Float64,nw)
for i = 1:nw
    nu[i] = 1.025e9 + (i-1)*50e6
end




# load psf fits file created by meqtrees
 file = FITS("../data/meerkat_m30_25pix.psf.fits")
#file = FITS("../data/makems_meerkat_Jeremy.MS.psf.channel.75ch.fits")
data = float64(read(file[1]))
close(file)
psfcube = squeeze(data,3)
psfavg = cubeaverage(psfcube,5)
#psf = cropcubexy(psfavg,255)
psfcube = psfavg

# load gray sky model fits file
# file = FITS("../data/cluster.fits")
file = FITS("../data/M31.fits")
data = float64(read(file[1]))
close(file)
sky0 = squeeze(data,3)
sky = sky2cube(sky0,nu)
datacube = cubefilter(sky,psfcube)
