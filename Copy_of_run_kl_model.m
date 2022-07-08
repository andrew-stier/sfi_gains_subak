clear;
addpath('/home/Dropbox/Dropbox/sfi_gains');
load("steady_state_kremer_lansing_spin.mat")
rng(2023);
N=100;
nrstates=4;
pestradius=2;
harvestradius=1;
temp=0.05;
nblock=4;
a=0.5;
b=20*a;
c=.01;
sigma = .8*b;
T=250;
shockrate = 1;
tF=5;
counter=50;
[spins,harvests] = temperature_Kremer_Lansing_Model(N,nrstates,pestradius,harvestradius,temp,nblock,T,a,b,tF,sigma, shockrate,counter,spin);
% figure();imagesc(spins{T});colorbar()

ht = [];
for t = 1:T
    ht = [ht mean(mean(harvests{t}(~isnan(harvests{t}))))];
end
figure()
plot(ht)

g = [];
for t=1:T
    htemp = harvests{t};
    htemp = htemp(~isnan(htemp));
    g = [g ginicoeff(htemp)];
end
figure()
plot(g)

d = [];
for t=1:T
    d = [d sum(sum(isnan(spins{t})))/N^2];
end
figure()
plot(d)

sp = spins{T};
sp(isnan(sp))=0;
imagesc(sp);colorbar()
