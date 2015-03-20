# proximity operators
function prox_visi(y_est::Array)
    return (2.*visi + rho_y.*y_est)./(rho_y + 2)
end

function prox_u(u::Array,μ)
    return (max(1-μ./abs(u),0).*u)
end


function conjgrad(xw::Array,bw::Array,mypsfw::Array,mypsfadjw::Array,mu::Float64;tol = 1e-6,itermax = 1e3)

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


#
# function conjgrad( M ;tol = 1e-6,itermax = 1e3)
#
#     r = M[z][2] - (imfilter_fft(imfilter_fft(M[z][1], M[z][3],"circular"), M[z][4],"circular") + M[z][5]*M[z][1])
#     rm = r
#     p = r
#     iter = 0
#     loop = true
#
#     while loop
#         iter += 1
#         Qp = imfilter_fft(imfilter_fft(p, M[z][3],"circular"), M[z][4],"circular") + M[z][5]*p
#
#         alpha = vecnorm(r)^2/sum(Qp.*p)
#         M[z][1] = M[z][1] + alpha*p
#         r = r - alpha*Qp
#         betaa = (vecnorm(r)/vecnorm(rm))^2
#         rm = r
#         p = r + betaa*p
#         crit = vecnorm(r)
#         if iter > itermax
#             error("itermax reached")
#         end
#         if crit < tol
#             loop = false
#         end
#     end
#     return M[z][1]
# end
