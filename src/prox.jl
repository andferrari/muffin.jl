# proximity operators
function prox_visi(y_est::Array)
    return (2.*visi + rho_y.*y_est)./(rho_y + 2)
end

function prox_u(u::Array,μ)
    return (max(1-μ./abs(u),0).*u)
end


@everywhere function gradD(x::Array,fty::Array,taup::Array,mypsf::Array,mypsfadj::Array,p::Array,z::Int64)



                    b = fty + taup + rhop*p

                    xm = x
                    r = b - (imfilter_fft(imfilter_fft(x, mypsf), mypsfadj) + (muesp + rhop)*eye(nfty)*x)  # r0
                    rm = r
                    p = r                                                 # p0

                    iter = 0
                    loop = true
                    figure(1)
                    subplot(5,2,z)
                    crit=Float64[]
                    while loop

                        iter += 1
                        alpha = sum(r.*r)/sum((imfilter_fft(imfilter_fft(p, mypsf), mypsfadj) + (muesp + rhop)*eye(nfty)*p).*p)
                        x = x + alpha*p
                        r = r - alpha*(imfilter_fft(imfilter_fft(p, mypsf), mypsfadj) + (muesp + rhop)*eye(nfty)*p)
                        betaa = sum(r.*r)/sum(rm.*rm)
                        rm = r
                        p = r + betaa*p
                        push!(crit,(norm(x - xm)))
                        plot(crit)


                        if iter != 1
                            if crit[iter-1] < crit[iter]
                            loop = false
                            end
                        end
                        xm = x
                    end

                    return x
            end
