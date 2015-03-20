# proximity operators
function prox_visi(y_est::Array)
    return (2.*visi + rho_y.*y_est)./(rho_y + 2)
end

function prox_u(u::Array,μ)
    return (max(1-μ./abs(u),0).*u)
end


@everywhere function conjgrad(xw::Array,bw::Array,mypsfw::Array,mypsfadjw::Array,mu::Float64;tol = 1e-6,itermax = 1e3)

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
