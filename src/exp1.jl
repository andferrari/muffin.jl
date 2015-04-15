using PyPlot

alpharecsky = zeros(Float64,nxy,nxy,3)
nu0 = (nu[end]+nu[1])/2


N = 15
ν = nu./nu0
A = [ones(N) log10(ν) log10(ν).^2]
for i = 1:nxy, j = 1:nxy

    if minimum(squeeze(squeeze(sky[i,j,:],1),1)) <= 0
    alpharecsky[i,j,:] = 0

    elseif minimum(squeeze(squeeze(sky[i,j,:],1),1)) > 0
    b = squeeze(squeeze(log10(sky[i,j,:]),1),1)
    c = A\b
    alpharecsky[i,j,1] = c[1]
    alpharecsky[i,j,2] = c[2]
    alpharecsky[i,j,3] = c[3]
    end
end




figure(6)
clf()
title("Cste Reconstruction")
colorbar(imshow(-alpharecsky[:,:,1]))

figure(7)
clf()
title("Alpha Reconstruction")
colorbar(imshow(-alpharecsky[:,:,2]))

figure(8)
clf()
title("Beta Reconstruction")
colorbar(imshow(-alpharecsky[:,:,3]))

figure(9)
clf()
title("Coupe alpha")
p = linspace(88,168,81)
q = alpharecsky[128,88:168,2]'
plot(p,q)


w = createobjtheo(alpharecsky[:,:,2],alpharecsky[:,:,3])
v = createobjtheo(alpharec[:,:,2],alpharec[:,:,3])
figure(10)
clf()
p = nu./nu0
m = 148
q = (squeeze(squeeze(x[128,m,:],1),1))

d = (squeeze(squeeze(sky[128,m,:],1),1))

s = (nu./nu0).^(alpharecsky[128,m,2]+alpharecsky[128,m,3].*log10(nu./nu0))

#maax = r[1]
#q = q/q[1]*maax
#s = s/s[1]*maax
#d = d/d[1]*maax


plot(log10(p),log10(d))
plot(log10(p),log10(q))
plot(log10(p),log10(s))


label1 = "Sky spectrum"
label2 = "Reconstructed spectrum"
label3 = "Fitted Sky spectrum"

legend( (label1, label2, label3), loc="upper right")
xlabel(r"$log(\nu./\nu0)$")


figure(9)
#clf()
title("Coupe alpha")
p = linspace(88,168,81)
q = alpharecsky[128,88:168,3]'
plot(p,q)
