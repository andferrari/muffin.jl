#################################
###### proximity operators ######
function prox_u(u::SharedArray,μ::Float64)
    return (max(1-μ./abs(u),0).*u)
end

function prox_u(u::Array,μ::Float64)
    return (max(1-μ./abs(u),0).*u)
end
##################################


##########################################
###### Conjugate Gradient Algorithm ######
function conjgrad(xw::Array{Float64,2},bw::Array{Float64,2},mypsfw::Array{Float64,2},
                  mypsfadjw::Array{Float64,2},mu::Float64; tol = 1e-6,itermax = 1e3)

    r = bw - (imfilter_fft(imfilter_fft(xw, mypsfw,"circular"), mypsfadjw,"circular") + mu*xw)
    rm = r
    p = r
    iter = 0
    loop = true


    while loop
        iter += 1

        Qp = imfilter_fft(imfilter_fft(p, mypsfw,"circular"), mypsfadjw,"circular") + mu*p

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
        end
    end
    return xw
end
##########################################


function estime_s(s::Array{Float64,3},tmp::Array{Float64,3},nxy::Int64,nspec::Int64,spectralwlt::Array{Float64,3},
                 x::SharedArray{Float64,3},taus::Array{Float64,3},rhov::Float64,rhos::Float64)
    for i in 1:nxy, j in 1:nxy
     spectralwlt[i,j,:]= idct(tmp[:,i,j])
    end
    s = (spectralwlt + rhos*x - taus)/(rhov*nspec + rhos)
    return s
end

function estime_sh(s::Array{Float64,3},sh::Array{Float64,3},nxy::Int64)
    vecs = permutedims(s,[3,1,2])
        for i in 1:nxy, j in 1:nxy
            sh[i,j,:] = dct(vecs[:,i,j])
        end
    return sh
end

##########################################
##########################################
function estime_ssh(s::Array{Float64,3},sh::Array{Float64,3},tmp::Array{Float64,3},nxy::Int64,nspec::Int64,spectralwlt::Array{Float64,3},
                 x::SharedArray{Float64,3},taus::Array{Float64,3},rhov::Float64,rhos::Float64)

    for i in 1:nxy, j in 1:nxy
        spectralwlt[i,j,:]= idct(tmp[:,i,j])
    end
    s = (spectralwlt + rhos*x - taus)/(rhov*nspec + rhos)

    vecs = permutedims(s,[3,1,2])

    for i in 1:nxy, j in 1:nxy
        sh[i,j,:] = dct(vecs[:,i,j])
    end

    return s,sh
end
