# initialisation

using FITSIO

# load psf fits

f = FITS("../data/makems_wb_meerkat_m30_4h60s.MS.psf.channel.32ch.fits")
close(f)
