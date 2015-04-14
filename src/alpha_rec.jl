####################################################################
####################### alpha reconstruction #######################
####################################################################

alpharec = zeros(Float64,nxy,nxy,3)
nu0 = (nu[end]+nu[1])/2
# X = log(nu/nu0)
# b = ones(15)
count = 0
countnonzero = 0
# for i = 1:nxy, j = 1:nxy
#
#     if minimum(squeeze(squeeze(x[i,j,:],1),1)) <= 0
#     c = 1
#
#     elseif minimum(squeeze(squeeze(x[i,j,:],1),1)) > 0
#     a = squeeze(squeeze(log(x[i,j,:]),1),1)
#     y = [a b]
#         if det(y'*y) != 0
#             alpharec[i,j,1] = (inv(y'*y)*y'*X)[1]
#             alpharec[i,j,2] = (inv(y'*y)*y'*X)[2]
#             countnonzero += 1
#         else
#             count += 1
#         end
#     end
# end

N = 15
ν = nu./nu0
A = [ones(N) log(ν) log(ν).^2]
for i = 1:nxy, j = 1:nxy
    b = squeeze(squeeze(log(x[i,j,:]),1),1)
    c = A\b
    alpharec[i,j,1] = c[1]
    alpharec[i,j,2] = c[2]
    alpharec[i,j,3] = c[3]
end


figure(6)
clf()
title("Cste Reconstruction")
colorbar(imshow(-alpharec[:,:,1]))

figure(7)
clf()
title("Alpha Reconstruction")
colorbar(imshow(-alpharec[:,:,2]))

figure(8)
clf()
title("Beta Reconstruction")
colorbar(imshow(-alpharec[:,:,3]))

# figure(7)
# clf()
# title("Real alpha")
# colorbar(imshow(alpha[:,:]))
