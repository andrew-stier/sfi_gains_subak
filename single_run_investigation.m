%clear;
addpath('/home/Dropbox/Dropbox/sfi_gains/sfi_gains_subak/');
load("steady_state_kremer_lansing_spin.mat")
rng(2022);
N=100;
nrstates=4;
pestradius=2;
harvestradius=1;
temp=0.05;
nblock=4;
a=0.5;
b=9.6;
T=400;
shockrate = 1;
counter=50;
hts = [];
gts = [];
dts = [];
sg_buckets = 20;
xis = zeros(10,sg_buckets);
tF = 3;
sg = 11.5;
display(sprintf('tf: %d, sg: %d',tF,sg));
sigma = sg/sg_buckets*b;
[spins,harvests,shocks] = temperature_Kremer_Lansing_Model(N,nrstates,pestradius,harvestradius,temp,nblock,T,a,b,tF,sigma, shockrate,counter,spin);
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

neighborradius = 1;
align_bins = 10;
alignment = 1:align_bins;
alignment = alignment/align_bins;
failure_pdfs = [];
for t=2:T
    failed = zeros(N,N);
    f_aligned = zeros(N,N);
        for i=1:N
            for j=1:N
                ilimit1=max(1,i-neighborradius);
                ilimit2=min(N,i+neighborradius);
                SpinNeigh=[]; % vector of state values in neigborhood to compute pest load
                for qq=ilimit1:ilimit2
                    width=pestradius-abs(qq-i);
                    jlimit1=max(1,j-width);
                    jlimit2=min(N,j+width);
                    SpinNeigh=[SpinNeigh spins{t-1}(qq,jlimit1:jlimit2)];
                end
                f_aligned(i,j) = sum(SpinNeigh == spins{t-1}(i,j))/length(SpinNeigh);
                failed(i,j) = isnan(spins{t}(i,j));
            end
        
        end
        p_fail = [];
        for b=1:align_bins
            if b>1
                p_fail = [p_fail sum(sum((f_aligned<=alignment(b)) & (f_aligned>alignment(b-1)) & failed))];
            else
                p_fail = [p_fail sum(sum((f_aligned<=alignment(b)) & failed))];
            end
        end
        p_fail = p_fail/sum(p_fail);
        failure_pdfs = [failure_pdfs; p_fail];
end
figure();plot(alignment,nansum(failure_pdfs)/nansum(nansum(failure_pdfs)));set(gca, 'YScale', 'log');xlabel('fraction of neighbors aligned');ylabel('probability of farm failure')

n_failed = [];
for t=1:T
    n_failed = [n_failed sum(sum(isnan(spins{t})))];
end

die_off = 5;

sp = harvests{8};sp(sp<0)=0;sp(1,1)=0;sp(1,2)=2.5;figure();imagesc(sp);colorbar()
% [sizes, ts] = PatchSize(spins{1});
% [sizes, ts] = PatchSize(spins{die_off-1});
% [sizes, ts] = PatchSize(spins{die_off});
% 
% max_size = log10(max(sizes));
% bins = 10;
% x = 1:bins;
% x = x/bins*max_size;
% probs = [];
% for i = 1:bins
%     if i>1
%         probs = [probs sum((log10(sizes)<=x(i)) & (log10(sizes)>x(i-1)))];
%     else
%         probs = [probs sum(log10(sizes)<=x(i))];
%     end
% end
% probs = probs/sum(probs);
% 
% figure();
% plot(x,probs);
% ylim([0 .37]);
% xlabel('log10(Patch Size)');
% ylabel('Pr(patch size)');