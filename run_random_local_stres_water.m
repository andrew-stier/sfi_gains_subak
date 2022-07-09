clear;
addpath('/home/Dropbox/Dropbox/sfi_gains/sfi_gains_subak/');
rng(2022)
N=100;
nrstates=4;
pestradius=2;
harvestradius=1;
temp=0.05;
nblock=4;
localwaterstress = .5;
a=0.5;
b=9.6;
T=200;
counter=50;
[spinl,harvest] = random_local_water_stress_Kremer_Lansing_Model(N,nrstates,pestradius,harvestradius,localwaterstress,temp,nblock,T,a,b,counter);
%save( 'steady_state_kremer_lansing_spin.mat', 'spin' )
figure()
imagesc(spinl)

load("steady_state_kremer_lansing_spin.mat")
[ MIl,Lstatl,xil ] = NormalizedCorreletionSpinLattice(spinl,5);
[ MI,Lstat,xi ] = NormalizedCorreletionSpinLattice(spin,5);