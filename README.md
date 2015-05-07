Muffin.jl
=========

`Muffin.jl` is a julia implementation of MUFFIN: Multi-Frequency Sparse Radio Interferometric imaging

see http://arxiv.org/abs/1504.06847

##Getting Started

####Installation

``Muffin.jl` uses the following Julia packages :

* FITSIO.jl : reading and writing FITS.
* Images.jl : Fast Fourier Transform implementation.
* Wavelets.jl : package for fast wavelet transforms. 
* HDF5.jl : for writing JLD ("Julia data") variables.
* GHF.jl : generation of Gaussian homogeneous spatial field.
* PyPlot.jl : provides a Julia interface to the Matplotlib plotting library

The installation is as simple as : 

.. code:: julia

 Pkg.clone("https://github.com/andferrari/GHF.jl.git")  
 Pkg.add("PyPlot")   
 Pkg.add("FITSIO")   
 Pkg.add("Images")    
 Pkg.add("Wavelets")   
 Pkg.add("HDF5")   
 
 
To install ``Muffin.jl``, type from a Julia session the following command :

	Pkg.clone("https://github.com/andferrari/Muffin.jl.git")





####Usage

To load the ``Muffin.jl`` module, type from a Julia session :

.. code:: julia

	using Muffin

To use parallel computing, start Julia with **nprocs** local process and load the module :

.. code:: julia

	$ julia -p nprocs  
	julia> using Muffin
	
	
	
##Functions

####Main Function 
``muffin(...)``

* **muffin** is called with parameters definition :

.. code:: julia


	  psfst, skyst, algost, admmst, toolst = muffin(nitermax, rhop, rhot, rhov, rhos, μt, μv, mueps)

* It returns 5 structures :
	* psfst : datas related to the PSF
	* skyst : datas related to the SKY
	* algost : contains algorithm parameters
	* admmst : datas related to the admm method
	* toolst : contains error calculation arrays
 











