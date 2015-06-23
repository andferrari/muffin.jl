####################################################################
####################### alpha reconstruction #######################
####################################################################
using PyPlot

alpharec = zeros(Float64,nxy,nxy,3)
alpharecsky = zeros(Float64,nxy,nxy,3)
alpharecdata = zeros(Float64,nxy,nxy,3)
nu0 = (nu[end]+nu[1])/2

N = 11
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

for i = 1:nxy, j = 1:nxy

    if minimum(squeeze(squeeze(mydata[i,j,:],1),1)) <= 0
    alpharecdata[i,j,:] = 0

elseif minimum(squeeze(squeeze(mydata[i,j,:],1),1)) > 0
    b = squeeze(squeeze(log10(mydata[i,j,:]),1),1)
    c = A\b
    alpharecdata[i,j,1] = c[1]
    alpharecdata[i,j,2] = c[2]
    alpharecdata[i,j,3] = c[3]
    end
end
