####################################################################
####################### alpha reconstruction #######################
####################################################################
using PyPlot

alpharec = zeros(Float64,nxy,nxy,3)
alpharecsky = zeros(Float64,nxy,nxy,3)
nu0 = (nu[end]+nu[1])/2

N = 15
ν = nu./nu0
A = [ones(N) log10(ν) log10(ν).^2]
for i = 1:nxy, j = 1:nxy

    if minimum(squeeze(squeeze(x[i,j,:],1),1)) <= 0
    alpharec[i,j,:] = 0

    elseif minimum(squeeze(squeeze(x[i,j,:],1),1)) > 0
    b = squeeze(squeeze(log10(x[i,j,:]),1),1)
    c = A\b
    alpharec[i,j,1] = c[1]
    alpharec[i,j,2] = c[2]
    alpharec[i,j,3] = c[3]
    end
end

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
subplot(121)
colorbar(imshow(-alpharec[:,:,1]))
subplot(122)
colorbar(imshow(-alpharecsky[:,:,1]))

figure(7)
clf()
title("Alpha Reconstruction")
subplot(121)
colorbar(imshow(-alpharec[:,:,2]))
subplot(122)
colorbar(imshow(-alpharecsky[:,:,2]))

figure(8)
clf()
title("Beta Reconstruction")
subplot(121)
colorbar(imshow(-alpharec[:,:,3]))
subplot(122)
colorbar(imshow(-alpharecsky[:,:,3]))

############################################# plot alpha
#############################################
figure(9)
clf()
p = linspace(88,168,81)
q = alpharecsky[128,88:168,2]'
r = alpharec[128,88:168,2]'
plot(p,q)
plot(p,r)
legend( (L"$\alpha$ sky", L"$\alpha$ reconstructed"), loc="upper center")
xlabel("pix")
############################################# plot beta
#############################################
figure(10)
clf()
p = linspace(88,168,81)
q = alpharecsky[128,88:168,3]'
r = alpharec[128,88:168,3]'
plot(p,q)
plot(p,r)
legend( (L"$\beta$ sky", L"$\beta$ reconstructed"), loc="upper right")
xlabel("pix")
#############################################
#############################################






############################################# plot spectre
#############################################
#w = createobjtheo(alpharecsky[:,:,2],alpharecsky[:,:,3])
#v = createobjtheo(alpharec[:,:,2],alpharec[:,:,3])

figure(11)
clf()
subplot(1,2,1)
p = nu./nu0
i = 128
j = 108

q = 1*(squeeze(squeeze(x[i,j,:],1),1))
d = (squeeze(squeeze(sky[i,j,:],1),1))
#s = (nu./nu0).^(alpharecsky[128,m,2]+alpharecsky[128,m,3].*log10(nu./nu0))

plot(log10(p),log10(d))
plot(log10(p),log10(q))
#plot(log10(p),log10(s))

label1 = "Sky spectrum"
label2 = "Reconstructed spectrum"
#label3 = "Fitted Sky spectrum"

legend( (label1, label2), loc="upper right")
xlabel(L"$log(\nu/\nu0)$")

subplot(1,2,2)
p = nu./nu0
i = 128
j = 148

q = 1.*(squeeze(squeeze(x[i,j,:],1),1))
d = (squeeze(squeeze(sky[i,j,:],1),1))
#s = (nu./nu0).^(alpharecsky[128,m,2]+alpharecsky[128,m,3].*log10(nu./nu0))

plot(log10(p),log10(d))
plot(log10(p),log10(q))

#plot(log10(p),log10(s))

label1 = "Sky spectrum"
label2 = "Reconstructed spectrum"
#label3 = "Fitted Sky spectrum"

legend( (label1, label2), loc="upper left")
xlabel(L"$log(\nu/\nu0)$")
#############################################
#############################################
