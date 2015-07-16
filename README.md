
=========

`Muffin.jl` is a julia implementation of MUFFIN: Multi-Frequency Sparse Radio Interferometric imaging

see http://arxiv.org/abs/1504.06847

#Getting Started

####Installation

`Muffin.jl` uses the following Julia packages :

* FITSIO.jl : reading and writing FITS.
* Images.jl : Fast Fourier Transform implementation.
* Wavelets.jl : package for fast wavelet transforms. 
* HDF5.jl : for writing JLD ("Julia data") variables.
* GHF.jl : generation of Gaussian homogeneous spatial field.
* PyPlot.jl : provides a Julia interface to the Matplotlib plotting library.

The installation is as simple as : 

```julia   
 Pkg.clone("https://github.com/andferrari/GHF.jl.git")  
 Pkg.add("PyPlot")   
 Pkg.add("FITSIO")   
 Pkg.add("Images")    
 Pkg.add("Wavelets")   
 Pkg.add("HDF5")    
```
 
 
To install `Muffin.jl`, type from a Julia session the following command :

```julia
	Pkg.clone("https://github.com/andferrari/Muffin.jl.git")
```





####Usage

To load the `Muffin.jl` module, type from a Julia session :

```julia
	using Muffin
```

To use parallel computing, start Julia with **nprocs** local process and load the module :

```julia
	$ julia -p nprocs
	julia> @everywhere using Muffin
```

You just need to add the keyword **parallel** in the MUFFIN function, for an example :

	
```julia
	julia> psfst, skyst, algost, admmst, toolst = muffin(nitermax = 10,parallel="true");
```
	
#Functions

####Main Function 
`muffin(...)`

* **muffin** is called with parameters definition :

```julia
	  psfst, skyst, algost, admmst, toolst = muffin(folder, dataobj, datapsf, nitermax, rhop, rhot, rhov, rhos, μt, μv, mueps)
```

* It returns 5 structures :
	* psfst : datas related to the PSF
	* skyst : datas related to the SKY
	* algost : contains algorithm parameters
	* admmst : datas related to the admm method
	* toolst : contains error calculation arrays
	
* Parameters are :
	* folder : path to the folder containing FITS files. Default : installation folder, containing FITS files for the demo
	* dataobj : FITS object.
	* datapsf : FITS psf.
	* nitermax : maximum number of  ADMM iterations. Default : `500`
	* rhop : ADMM parameter for positivity constraint. Default : `1`
	* rhot : ADMM parameter for spatial constraint. Default : `5`
	* rhov : ADMM parameter for spectral constraint. Default : `2`
	* rhos : ADMM parameter for spectral constraint. Default : `1`
	* μt : Spatial regularization parameter. Default : `5e-1`
	* μv : Spectral regularization parameter. Default : `1`
	* mueps : Ridge/Tikhonov regularization parameter : `1e-3`
	* ws : warm start
	* parallel : to start parallel Muffin
 
 
####Save data
`savedata(...)`

* **savedata** is called by :

```julia
	savedata(savepath, psfst, skyst, algost, admmst, toolst)
```

where `savepath` is the path/JLD file name where datas will be saved. 


#Muffin Demo

`Muffin.jl` contains a demo file `muffindemo.jl` in the default installation folder.
 
To run the demo, type :

	using Muffin
	demo = string(Pkg.dir("Muffin"),"/src/muffindemo.jl")
	include(demo);
	
If you want to visualize your results, just type after the simulation :

	plot = string(Pkg.dir("Muffin"),"/src/readdata.jl")
	include(plot);


The demo will run with the following parameters :

	myfolder = string(Pkg.dir("Muffin"))
	mydataobj =  "data/M31.fits"
	mydatapsf =  "data/meerkat_m30_25pix.psf.fits"
	mynitermax = 1
	myrhop = 1
	myrhot = 5
	myrhov = 2
	myrhos = 1
	myμt = 5e-1
	myμv = 1
	mymueps = 1e-3
	mysavepath = string(myfolder,myfolder[1],"data/demo_results.jld")
	
and the results will be saved in the `demo_results.jld` file.










