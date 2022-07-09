clear;
addpath('/home/Dropbox/Dropbox/sfi_gains/sfi_gains_subak/');
load("steady_state_kremer_lansing_spin.mat")
rng(2023);
N=100;
nrstates=4;
pestradius=2;
harvestradius=1;
temp=0.05;
nblock=4;
a=0.5;
b=9.6;
T=100;
shockrate = 1;
counter=50;
hts = [];
gts = [];
dts = [];
xis = zeros(10,100);
for tF = 1:10
    for sg = 1:5
        display(sprintf('tf: %d, sg: %d',tF,sg));
        sigma = sg/5*b;
        [spins,harvests] = temperature_Kremer_Lansing_Model(N,nrstates,pestradius,harvestradius,temp,nblock,T,a,b,tF,sigma, shockrate,counter,spin);
        ht = [];
        for t = 1:T
            ht = [ht mean(mean(harvests{t}(~isnan(harvests{t}))))];
        end
        g = [];
        for t=1:T
            htemp = harvests{t};
            htemp = htemp(~isnan(htemp));
            g = [g ginicoeff(htemp)];
        end
        d = [];
        for t=1:T
            d = [d sum(sum(isnan(spins{t})))/N^2];
        end
        hts = [hts;ht];
        gts = [gts;g];
        dts = [dts;d];        
        sp = spins{T};
        sp(isnan(sp)) = -999;
        [ MI,Lstat,xi ] = NormalizedCorreletionSpinLattice(sp,5);
        xis(tF,sg) = xi;
    end
end