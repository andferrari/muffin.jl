@everywhere function myidwt(M)
 idwt(M[1],M[2])
end

# proximity operators
function prox_visi(y_est::Array)
    return (2.*visi + rho_y.*y_est)./(rho_y + 2)
end

function prox_u(u::Array,μ)
    return (max(1-μ./abs(u),0).*u)
end


function estime_x(x)

 yy = tau_y + rho_y*yc
 # tt = tau_T + rho_T*T
 # pp = tau_P + rho_P*P
 # ss = tau_S + rho_S*S

 #reg = pp + ss

 for z = 1:nw
  x1[:,:,z] = nfft_adjoint(Plan[z],yy[:,z])/N
  # for i = 1:nspat
  #  x2[:,:,z,i] = idwt(tt[:,:,z,i], wavelet(SpatialWaveletsBasisList[i]))
  # end
 end

println("debut_xtt","  ")

  xtt = x1  #+ reg  + sum(x2,4)

println("fin_xtt")



K=1
# K = μ_ϵ +rho_P + rho_S + 9*rho_T
 for z = 1:nw
  xxest = matinv[:,:,z]*nfft(Plan[z],xtt[:,:,z])
  xest = nfft_adjoint(Plan[z],xxest)
  xx = (xtt[:,:,z] - rho_y/K*xest)/K
  x[:,:,z] = real(xx)
  Fx[:,z] = nfft(Plan[z],xx)
 end
  return x,Fx
end

# function estime_S(S)
#   for z=1:nw
#    S[:,:,z] = (idct(tau_V[:,:,z]/rho_V + V[:,:,z]) + x[:,:,z] - tau_S[:,:,z]/rho_S)/(rho_V*9+rho_S)
#    sH[:,:,z] = dct(S[:,:,z])
#   end
#   return S,sH
# end
function estime_S(S)
  for z=1:nw
   S[:,:,z] = (idct(tau_V[:,:,z] + rho_V*V[:,:,z]) + rho_S*x[:,:,z] - tau_S[:,:,z])/(rho_V*nspat+rho_S)
   sH[:,:,z] = dct(S[:,:,z])
  end
  return S,sH
end
#update of v  S




@everywhere function MyPlan(Mtmp)
 NFFTPlan(Mtmp[1],(Mtmp[2],Mtmp[2]))
end

function nudft3d(tab_u::Array,tab_v::Array)
# 3D nudft [number of bases (spatiale frequencies) per wavelength * number of pixels * number of wavelength]
# U and V must be of the form : [number of bases * number of wavelength]
# defined as in equation 3 of PAINTER
if N<1||N<1
  error("nudft3d:Image size must be positive")
end
  nfu  	= length(tab_u)
  nfv  	= length(tab_v)
  if(nfu==nfv&&nfu>0&&nfv>0)
    # test if U and V have same size and have same number of wavelength, otherwise PB
    size_tab_u = size(tab_u)
    size_tab_v = size(tab_v)
    if(size_tab_u[2]==size_tab_v[2])
      (Nb,nw)= size_tab_u
	  # # # SERIAL
	  #       F3D = complex(zeros(Nb,N*N,nw))
	  # tic()
	  #       for n      = 1:nw
	  #         nudftmat   = non_uniform_dft(tab_u[:,n],tab_v[:,n])
	  #         F3D[:,:,n] = nudftmat
	  #       end
	  # toc()
	  # tic()
	  UVMAT  = {(tab_u[:,n],tab_v[:,n],N,Nb) for n =1:nw }
	  F3D    = pmap(non_uniform_dft_par,UVMAT)
	  # toc()
    return F3D
    else
    error("nudft3d: U ad V vectors have not of same number of wavelength")
    end
  else
    error("nudft3d: U ad V vectors are not of same size")
  end
end

function non_uniform_dft(u_in::Array,v_in::Array)
# Non Uniform DFT Matrix:
# for a N*N pixels image
# u, v vector of normalized spatial frequencies (between -0.5 and 0.5)
u 		= vec(u_in)
v 		= vec(v_in)

kx 		= [-N/2:-1+N/2]
k1  	= repmat(kx,1,N)
k2  	= k1'
k1v 	= vec(k1)'/N
k2v 	= vec(k2)'/N
k1m 	= repmat(k1v,Nb,1)
k2m 	= repmat(k2v,Nb,1)

um 		= repmat(u,1,N*N)
vm 		= repmat(v,1,N*N)

F    	= exp(-2im*pi*(k1m.*um+k2m.*vm))./N
end

@everywhere function non_uniform_dft_par(UVMAT)
#UVMAT: u_in::Array,v_in::Array
# Non Uniform DFT Matrix:
# for a N*N pixels image
# u, v vector of normalized spatial frequencies (between -0.5 and 0.5)
u 		= vec(UVMAT[1])
v 		= vec(UVMAT[2])
N 		= UVMAT[3]
Nb 		= UVMAT[4]

kx 		= [-N/2:-1+N/2]
k1  	= repmat(kx,1,N)
k2  	= k1'
k1v 	= vec(k1)'/N
k2v 	= vec(k2)'/N
k1m 	= repmat(k1v,Nb,1)
k2m 	= repmat(k2v,Nb,1)

um 		= repmat(u,1,N*N)
vm 		= repmat(v,1,N*N)

return exp(-2im*pi*(k1m.*um+k2m.*vm))./N
end
