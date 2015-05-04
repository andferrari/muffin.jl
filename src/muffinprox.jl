#################################
###### proximity operators ######
function prox_u(u::SharedArray,μ)
    return (max(1-μ./abs(u),0).*u)
end

function prox_u(u::Array,μ)
    return (max(1-μ./abs(u),0).*u)
end
##################################

function forconjgrad(xw::SharedArray{Float64,3},bw::Array,mypsfw::Array,mypsfadjw::Array,mu::Float64,wlt::SharedArray{Float64,3},nfreq::Int64)
    @sync @parallel for z in 1:nfreq
                    xw[:,:,z] = conjgrad(xw[:,:,z], bw[:,:,z] + wlt[:,:,z], mypsfw[:,:,z], mypsfadjw[:,:,z], mu, tol=1e-4, itermax = 1e3)
                    end
    return xw
end

##########################################
###### Conjugate Gradient Algorithm ######
function conjgrad(xw::Array,bw::Array,mypsfw::Array,mypsfadjw::Array,mu::Float64;tol = 1e-6,itermax = 1e3)
    tic()
    r = bw - (imfilter_fft(imfilter_fft(xw, mypsfw,"circular"), mypsfadjw,"circular") + mu*xw)
    rm = r
    p = r
    iter = 0
    loop = true
    a = toq()
    # println("init conjgrad", "   ", a)

    while loop
        iter += 1
        tic()
        Qp = imfilter_fft(imfilter_fft(p, mypsfw,"circular"), mypsfadjw,"circular") + mu*p
        a = toq()
        # println("iteration int conjgrad","  ",a)

        alpha = vecnorm(r)^2/sum(Qp.*p)
        xw = xw + alpha*p
        r = r - alpha*Qp
        betaa = (vecnorm(r)/vecnorm(rm))^2
        rm = r
        p = r + betaa*p
        crit = vecnorm(r)
        if iter > itermax
            error("itermax reached")
        end
        if crit < tol
            loop = false
            # println("iteration conjgrad","  ",iter)
        end
    end
    return xw
end
##########################################

##########################################
# function estime_s(s,tmp)
#     for i in 1:nxy, j in 1:nxy
#      admmst.spectralwlt[i,j,:]= idct(tmp[:,i,j])
#     end
#     s = (admmst.spectralwlt + rhos*admmst.x - admmst.taus)/(rhov*nspec + rhos)
#     return s
# end
#
# function estime_sh(s)
#     vecs = permutedims(s,[3,1,2])
#     for i in 1:nxy, j in 1:nxy
#      admmst.sh[i,j,:]   = dct(vecs[:,i,j] )
#     end
#     return admmst.sh
# end
##########################################


function estime_s(s,tmp,nxy,nspec,admmst)
    for i in 1:nxy, j in 1:nxy
     admmst.spectralwlt[i,j,:]= idct(tmp[:,i,j])
    end
    s = (admmst.spectralwlt + admmst.rhos*admmst.x - admmst.taus)/(admmst.rhov*nspec + admmst.rhos)
    return s
end

function estime_sh(s,nxy,admmst)
    vecs = permutedims(s,[3,1,2])
    for i in 1:nxy, j in 1:nxy
     admmst.sh[i,j,:]   = dct(vecs[:,i,j] )
    end
    return admmst.sh
end
