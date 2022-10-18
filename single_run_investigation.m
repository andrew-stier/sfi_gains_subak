%clear;
%addpath('.../sfi_gains_subak/');
load("steady_state_kremer_lansing_spin.mat");
rng(2022);
N=100;
nrstates=4;
pestradius=2;
harvestradius=1;
temp=0.0;
nblock=4;
a=0.5;
b=9.6;
T=800;
shockrate = 1;
counter=50;
hts = [];
gts = [];
dts = [];
sg_buckets = 20;
xis = zeros(10,sg_buckets);
tF = 3;
sg = 16.5; %11.5
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

neighborradius = 2;
align_bins = 10;
max_distance = 20;
distance = 0:max_distance;
psize_bins = 15;
psize = 1:psize_bins;
psize = psize/psize_bins;
alignment = 1:align_bins;
alignment = alignment/align_bins;
alignment = [0 alignment];
failure_pdfs = [];
failure_pdfs_distance = [];
failure_pdfs_psize = [];
% code joris
n_farms_want_to_switch_time_t = zeros(1, T);
n_times_want_switches_per_farm = ones(N,N);
% end code joris
failure_size_distance = {};
for t=2:T
    failed = zeros(N,N);
    failed_prev = zeros(N,N);
    f_aligned = zeros(N,N);
    patch_size = zeros(N,N);
    patch_boundary = zeros(N,N);
    [size,ts,clusters]=PatchSize(spins{t-1});
    clusters = reshape(clusters,N,N);
    for i=1:N
        for j=1:N
            spself = spins{t-1}(i,j);  % added joris
            harvself = harvests{t-1}(i,j); % added joris
            ilimit1=max(1,i-neighborradius);
            ilimit2=min(N,i+neighborradius);
            SpinNeigh=[]; % vector of state values in neigborhood to compute pest load
            for qq=ilimit1:ilimit2
                width=neighborradius-abs(qq-i);
                jlimit1=max(1,j-width);
                jlimit2=min(N,j+width);
                SpinNeigh=[SpinNeigh spins{t-1}(qq,jlimit1:jlimit2)];
                % code joris
                for rr=jlimit1:jlimit2
                    spneigh = spins{t-1}(qq,rr);
                    harvneigh = harvests{t-1}(qq,rr);
                    % find if harvest neighbours is higher with
                    % different cropping pattern and increment counters
                    if (spneigh ~= spself) && (harvneigh > harvself)
                        n_times_want_switches_per_farm(i,j) = n_times_want_switches_per_farm(i,j) + 1;
                        n_farms_want_to_switch_time_t(t) = n_farms_want_to_switch_time_t(t) + 1;
                        break;
                    end
                end
                % end code joris
            end
            f_aligned(i,j) = sum(SpinNeigh == spins{t-1}(i,j))/length(SpinNeigh);
            failed(i,j) = isnan(spins{t}(i,j));
            failed_prev(i,j) = isnan(spins{t-1}(i,j));
            cid = clusters(i,j);
            if cid>0
                patch_size(i,j) = size(cid);
                patch_boundary(i,j) = length(find(f_aligned(clusters == cid)<1));
            else
                patch_size(i,j) = nan;
                patch_boundary = nan;
            end
        end        
    end
    distance_to_periphary = zeros(N,N);
    for i=1:N
        for j=1:N
            r_search = 0;
            found_periphary = 0;
            if f_aligned(i,j)<1
                distance_to_periphary(i,j) = 0;
            else
                while found_periphary==0
                    r_search = r_search+1;
                    ilimit1=max(1,i-r_search);
                    ilimit2=min(N,i+r_search);
                    alignNeigh=[]; % vector of state values in neigborhood to compute pest load
                    for qq=ilimit1:ilimit2
                        width=r_search-abs(qq-i);
                        jlimit1=max(1,j-width);
                        jlimit2=min(N,j+width);
                        alignNeigh=[alignNeigh f_aligned(qq,jlimit1:jlimit2)];
                    end

                    if length(find(alignNeigh<1))>0
                        distance_to_periphary(i,j) = r_search;
                        found_periphary=1;
                    end
                end
            end    
        end
    end

    p_fail = [];
    for b=1:length(alignment)
        if b>1
            p_fail = [p_fail sum(sum((f_aligned<=alignment(b)) & (f_aligned>alignment(b-1)) & failed & ~failed_prev))];
        else
            p_fail = [p_fail sum(sum((f_aligned<=alignment(b)) & failed & ~failed_prev))];
        end
    end 
    failure_pdfs = [failure_pdfs; p_fail];

    p_fail_distance = [];
    for b=1:length(distance)
        if b>1
            p_fail_distance = [p_fail_distance nansum(nansum((distance_to_periphary<=distance(b)) & (distance_to_periphary>distance(b-1)) & failed & ~failed_prev))];
        else
            p_fail_distance = [p_fail_distance nansum(nansum((distance_to_periphary<=distance(b)) & failed & ~failed_prev))];
        end
    end
    failure_pdfs_distance= [failure_pdfs_distance; p_fail_distance];
    
    p_boundary_size = patch_boundary./patch_size;
    p_fail_psize= [];
    for b=1:length(psize)
        if b>1
            p_fail_psize= [p_fail_psize nansum(nansum(((p_boundary_size<=psize(b)) & (p_boundary_size>psize(b-1)) & failed & ~failed_prev)))];
        else
            p_fail_psize= [p_fail_psize nansum(nansum(((p_boundary_size<=psize(b)) & failed & ~failed_prev)))];
        end
    end
    failure_pdfs_psize = [failure_pdfs_psize; p_fail_psize];    
end

figure();plot(alignment,nansum(failure_pdfs)/nansum(nansum(failure_pdfs)));set(gca, 'YScale', 'log');xlabel('fraction of neighbors aligned');ylabel('probability of farm failure')
figure();plot(alignment,nansum(failure_pdfs)/nansum(nansum(failure_pdfs)));set(gca, 'YScale', 'log');xlabel('fraction of neighbors aligned');ylabel('probability of farm failure')
% code joris
figure();plot(n_farms_want_to_switch_time_t / (N*N));xlabel('Timestep');ylabel('Fraction of farmers that want to switch')
figure();imagesc(n_times_want_switches_per_farm / T)
% end code joris
figure();plot(alignment,nansum(failure_pdfs)/nansum(nansum(failure_pdfs)));set(gca, 'YScale', 'log');xlabel('fraction of neighbors aligned');ylabel('probability of farm failure')
figure();plot(distance,nansum(failure_pdfs_distance)/nansum(nansum(failure_pdfs_distance)));set(gca, 'YScale', 'log');xlabel('distance to periphary');ylabel('probability of farm failure')
figure();plot(psize,nansum(failure_pdfs_psize)/nansum(nansum(failure_pdfs_psize)));set(gca, 'YScale', 'log');xlabel('patch boundary / patch size');ylabel('probability of farm failure')


figure();plot(dts);

[size,ts,clusters]=PatchSize(spins{T});

n_failed = [];
for t=1:T
    n_failed = [n_failed sum(sum(isnan(spins{t})))];
end

die_off = 5;

sp = spins{T};sp(sp<0)=0;sp(1,1)=0;sp(1,2)=2.5;figure();imagesc(sp);colorbar()

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