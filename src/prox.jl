# proximity operators
function prox_visi(y_est::Array)
    return (2.*visi + rho_y.*y_est)./(rho_y + 2)
end

function prox_u(u::Array,μ)
    return (max(1-μ./abs(u),0).*u)
end


function
b = fty + taup + rhop*p

x[:,:,:] = 1
xm = x                                                        # x0
r = b - (cubefilter(cubefilter(x, mypsf), mypsfadj))  # r0
rm = r
p = r                                                # p0

for z = 1:nfreq
    while vecnorm(x[:,:,z] - xm[:,:,z], 2) > 1E-8

        alpha = (r[:,:,z],r[:,:,z])/dot((cubefilter(cubefilter(p[:,:,z], mypsf), mypsfadj)),p[:,:,z])
        x = x + alpha*p
        r = r - alpha*(cubefilter(cubefilter(p, mypsf), mypsfadj))
        betaa = dot(r[:,:,z],r[:,:,z])/dot(rm[:,:,z],rm[:,:,z])
        rm = r
        p = r + betaa*p
        println(vecnorm(x[:,:,z] - xm[:,:,z])
        xm = x
    end
end
