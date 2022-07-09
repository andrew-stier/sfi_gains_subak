clear;
addpath('/home/Dropbox/Dropbox/sfi_gains/sfi_gains_subak/');

N=100;
nrstates=4;
pestradius=2;
harvestradius=1;
temp=0.05;
nblock=4;
a=0.5;
b=9.6;
T=400;
counter=50;
[spin,harvest] = Kremer_Lansing_Model(N,nrstates,pestradius,harvestradius,temp,nblock,T,a,b,counter);
save( 'steady_state_kremer_lansing_spin.mat', 'spin' )
figure(1)
imagesc(spin)

