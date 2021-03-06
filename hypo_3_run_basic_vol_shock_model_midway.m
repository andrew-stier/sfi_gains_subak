%addpath('../sfi_gains_subak/');
load("steady_state_kremer_lansing_spin.mat")
seed=99;
rng(seed);
xis = {};
params = {};
hts_all = {};
gts_all = {};
dts_all = {};
avg_same_crop_before_failure = [];
var_same_crop_before_failure = [];
avg_time_in_pattern_no_fail = [];
var_avg_time_in_pattern_no_fail = [];
avg_time_in_pattern_fail = [];
var_avg_time_in_pattern_fail = [];
failure_pdfs_psize_all = {};
failure_pdfs_distance_all = {};
failure_pdfs_all = {};
failed_t_all = {};
water_vector_all = {};    
pest_vector_all = {};
die_vector_all = {};
average_water_t_all = {};
average_pest_t_all = {};
increase_average_water_all = {};
increase_average_pest_all = {};
sg_buckets = 5; %shock size, 1-5, was 1-10
tf_max = 4; % amount of H<0 before fail, 1-4, was 1-10
cnt = 1;
for tF = 1:tf_max
    for sg = 1:sg_buckets
        params{cnt} = [tF sg];
        cnt = cnt +1;
    end
end
for idx = 1:(cnt-1)
    %%%%%%%%%%%%%%%%%%%%%%%%
    %setup variables
    N=100;
    spin = spin(1:N, 1:N);
    nrstates=4;
    pestradius=2;
    harvestradius=1;
    temp=0.05;
    nblock=4;
    a=0.5;
    b=9.6;
    T=500;
    shockrate = 1;
    counter=50;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    prms = params{idx};
    tF = prms(1);
    sg = prms(2);
    display(sprintf('tf: %d, sg: %d',tF,sg));
    sigma = sg/sg_buckets*b;
    [spins,harvests,shocks,pests,waters] = temperature_Kremer_Lansing_Model(N,nrstates,pestradius,harvestradius,temp,nblock,T,a,b,tF,sigma, shockrate,counter,spin);    
    shocks_all{idx} = shocks;
    sp = spins{T};
    sp(isnan(sp)) = -999;
    [ MI,Lstat,xi ] = NormalizedCorreletionSpinLattice(sp,5);
    xis{idx} = xi;
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
    dts_all{idx} = d;
    gts_all{idx} = g;
    hts_all{idx} = ht;
    
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
    num_steps_same_crop_before_fail = NaN(N,N);
    n_switches = ones(N,N);
    t_to_fail = T*ones(N,N);
    
    for i=1:N
        for j=1:N
            for t=2:T
                if isnan(spins{t}(i,j))
                    t_failure = t;
                    t_to_fail(i,j) = t_failure;
                    last_crop = spins{t_failure-1}(i,j);
                    for t_prev=t_failure-1:-1:1
                        if spins{t_prev}(i,j) ~= last_crop
                            num_steps_same_crop_before_fail(i,j) = t_failure - t_prev;
                            break;
                        elseif t_prev == 1
                            num_steps_same_crop_before_fail(i,j) = t_failure - t_prev;
                            break;
                        end
                    end
                    break;
                end
            end
            for t=2:T
                last_crop = spins{t_to_fail(i,j)-1}(i,j);
                for t_prev=t_to_fail(i,j)-1:-1:1
                    prev_crop = spins{t_prev}(i,j);
                    if last_crop ~= prev_crop
                        n_switches(i,j) = n_switches(i,j) + 1;
                    end
                    last_crop = prev_crop;
                end
                break;
            end
        end
    end
    
    %figure();pcolor(num_steps_same_crop_before_fail);title('Number of steps with same cropping pattern before failure');colorbar()
    %figure();pcolor(n_switches);title('Number of switches per farmer');colorbar()
    %figure();pcolor(t_to_fail);title('Time to failure');colorbar()
    
    
    %figure();plot(alignment,nansum(failure_pdfs)/nansum(nansum(failure_pdfs)));set(gca, 'YScale', 'log');xlabel('fraction of neighbors aligned');ylabel('probability of farm failure')
    
    avg_switches_before_fail = NaN(N,N);
    for i=1:N
        for j=1:N
            avg_switches_before_fail(i,j) = t_to_fail(i,j) / n_switches(i,j);
        end
    end
    
    figure();imagesc(avg_switches_before_fail);title('Average time in cropping pattern');colorbar()
    
    av_switches_fail = [];
    av_switches_no_fail = [];
    for i=1:N
        for j=1:N
            if t_to_fail(i,j) < T
                av_switches_fail = [av_switches_fail avg_switches_before_fail(i,j)];
            else
                av_switches_no_fail = [av_switches_no_fail avg_switches_before_fail(i,j)];
            end
        end
    end
    
    
    avg_same_crop_before_failure(prms(1), prms(2)) = nanmean(num_steps_same_crop_before_fail(:));
    var_same_crop_before_failure(prms(1), prms(2)) = nanvar(num_steps_same_crop_before_fail(:));
    
    avg_time_in_pattern_no_fail(prms(1), prms(2)) = nanmean(av_switches_no_fail);
    var_avg_time_in_pattern_no_fail(prms(1), prms(2)) = nanvar(av_switches_no_fail);
    
    avg_time_in_pattern_fail(prms(1), prms(2)) = nanmean(av_switches_fail);
    var_avg_time_in_pattern_fail(prms(1), prms(2)) = nanvar(av_switches_fail);



        for t=2:T
        failed = zeros(N,N);
        failed_t=zeros(1,T);
        averagepest_t=zeros(1,T);
        increase_average_pest=zeros(1,T);
        averagewater_t=zeros(1,T);
        increase_average_water=zeros(1,T);
        failed_prev = zeros(N,N);
        f_aligned = zeros(N,N);
        patch_size = zeros(N,N);
        patch_boundary = zeros(N,N);
        [size,~,clusters]=PatchSize(spins{t-1});
        clusters = reshape(clusters,N,N);
        for i=1:N
            for j=1:N
                ilimit1=max(1,i-neighborradius);
                ilimit2=min(N,i+neighborradius);
                SpinNeigh=[]; % vector of state values in neigborhood to compute pest load
                for qq=ilimit1:ilimit2
                    width=neighborradius-abs(qq-i);
                    jlimit1=max(1,j-width);
                    jlimit2=min(N,j+width);
                    SpinNeigh=[SpinNeigh spins{t-1}(qq,jlimit1:jlimit2)];
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
        for bb=1:length(alignment)
            if bb>1
                p_fail = [p_fail sum(sum((f_aligned<=alignment(bb)) & (f_aligned>alignment(bb-1)) & failed & ~failed_prev))];
            else
                p_fail = [p_fail sum(sum((f_aligned<=alignment(bb)) & failed & ~failed_prev))];
            end
        end 
        failure_pdfs = [failure_pdfs; p_fail];
    
        p_fail_distance = [];
        for bb=1:length(distance)
            if bb>1
                p_fail_distance = [p_fail_distance nansum(nansum((distance_to_periphary<=distance(bb)) & (distance_to_periphary>distance(bb-1)) & failed & ~failed_prev))];
            else
                p_fail_distance = [p_fail_distance nansum(nansum((distance_to_periphary<=distance(bb)) & failed & ~failed_prev))];
            end
        end
        failure_pdfs_distance= [failure_pdfs_distance; p_fail_distance];
        
        p_boundary_size = patch_boundary./patch_size;
        p_fail_psize= [];
        for bb=1:length(psize)
            if bb>1
                p_fail_psize= [p_fail_psize nansum(nansum(((p_boundary_size<=psize(bb)) & (p_boundary_size>psize(bb-1)) & failed & ~failed_prev)))];
            else
                p_fail_psize= [p_fail_psize nansum(nansum(((p_boundary_size<=psize(bb)) & failed & ~failed_prev)))];
            end
        end
        failure_pdfs_psize = [failure_pdfs_psize; p_fail_psize];  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        failed_t(1,t)=(sum(sum(failed))-sum(failed_t(1,:)));
        averagepest_t(1,t-1)=mean(pests{t}(~isnan(spins{t})),"omitnan");
        if t>3
            increase_average_pest(1,t-1)=(averagepest_t(1,t-1)-averagepest_t(1,t-2));
        end
        averagewater_t(1,t-1)=mean(waters{t}(~isnan(spins{t})),"omitnan");
        if t>3
            increase_average_water(1,t-1)=(averagewater_t(1,t-1)-averagewater_t(1,t-2));
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    failure_pdfs_psize_all{idx} = failure_pdfs_psize;
    failure_pdfs_distance_all{idx} = failure_pdfs_distance;
    failure_pdfs_all{idx} = failure_pdfs;
    failed_t_all{idx} = failed_t;
    average_water_t_all{idx} = averagewater_t;
    average_pest_t_all{idx} = averagepest_t;
    increase_average_water_all{idx} = increase_average_water;
    increase_average_pest_all{idx} = increase_average_pest;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pest_per_ind=zeros(N,N);    
    water_per_ind=zeros(N,N);  
    for i=1:N    
        for j=1:N
            pest_t=zeros(1,T-1);
            water_t=zeros(1,T-1);
            for t=2:T
                if ~isnan(spins{t}(i,j))
                    pest_t(1,t-1)=pests{t}(i,j);
                    water_t(1,t-1)=waters{t}(i,j);
                else
                    pest_t(1,t-1)=nan;
                    water_t(1,t-1)=nan;
                end
            end
            pest_per_ind(i,j)=mean(pest_t(1,:),"omitnan");
            water_per_ind(i,j)=mean(water_t(1,:),"omitnan");
        end
    end  
    pest_vector=reshape(pest_per_ind.',1,[]);
    die_vector=reshape(spins{T}.',1,[]);
    water_vector=reshape(water_per_ind.',1,[]);    
    water_vector_all{idx} = water_vector;    
    pest_vector_all{idx} = pest_vector;
    die_vector_all{idx} = die_vector;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

clear g bb failure_pdfs failure_pdfs_distance failure_pdfs_psize p_fail_psize p_boundary_size psize failed failed_prev p_fail_distance;
clear distance_to_periphary distance p_fail_distance p_fail alignment found_periphary alignNeigh SpinNeigh r_search qq;
clear ilimit1 ilimit2 jlimit1 jlimit2 width f_aligned patch_size patch_boundary size clusters htemp g d ht;
clear harvests spins shocks pests waters;
clear sigma sp MI Lstat xi pest_vector die_vector pest_t pest_per_ind water_vector water_per_ind water_t;
clear averagepest_t increase_average_pest averagewater_t increase_average_water counter cid tF spin sg harvest cnt;

save(sprintf('volatility_shock_run_kl_model_seed_%d.mat',seed))