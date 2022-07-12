clear;
addpath('../sfi_gains_subak/');
load("steady_state_kremer_lansing_spin.mat")
rng(seed);
N=100;
nrstates=4;
pestradius=2;
harvestradius=1;
temp=0.05;
nblock=4;
a=0.5;
b=9.6;
T=800;
shockrate = 1;
counter=50;
sg_buckets = 20;
tf_max = 10;
xis = {};
hts = [];
gts = [];
dts = [];
spins_all = {};
harvests_all = {};
pests_all = {};
shocks_all = {};
params = {};
cnt = 1;
for tF = 1:tf_max
    for sg = 1:sg_buckets
        params{cnt} = [tF sg];
        cnt = cnt +1;
    end
end
parfor idx = 1:(cnt-1)
    prms = params{idx};
    tF = prms(1);
    sg = prms(2);
    display(sprintf('tf: %d, sg: %d',tF,sg));
    sigma = sg/sg_buckets*b;
    [spins,harvests,shocks,pests] = temperature_Kremer_Lansing_Model(N,nrstates,pestradius,harvestradius,temp,nblock,T,a,b,tF,sigma, shockrate,counter,spin);
    spins_all{idx} = spins;
    harvests_all{idx} = harvests;
    shocks_all{idx} = shocks;
    pests_all{idx} = pests;     
    sp = spins{T};
    sp(isnan(sp)) = -999;
    [ MI,Lstat,xi ] = NormalizedCorreletionSpinLattice(sp,5);
    xis{idx} = xi;
end
save(sprintf('../volatility_shock_run_kl_model_seed_%d.mat',seed))