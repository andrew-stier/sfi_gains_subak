load("volatility_shock_run_kl_model_seed_99.mat")
figure();
imagesc(avg_time_in_pattern_no_fail);
c=colorbar();
xlabel('sg');
yticks((1:10)*20-10);
yticklabels({'tF=1','tF=2','tF=3','tF=4'});
ylabel('tf');c.Label.String = 'average time same cropping pattern no failure';


figure();
imagesc(avg_time_in_pattern_fail);
c=colorbar();
xlabel('sg');
yticks((1:10)*20-10);
yticklabels({'tF=1','tF=2','tF=3','tF=4'});
ylabel('tf');c.Label.String = 'average time same cropping pattern failure';

figure();
imagesc(avg_same_crop_before_failure);
c=colorbar();
xlabel('timestep');
yticks((1:10)*20-10);
yticks((1:10)*20-10);
yticklabels({'tF=1','tF=2','tF=3','tF=4'});
xlabel('sg');
ylabel('Time to failure');c.Label.String = 'average time same cropping pattern before failure';

figure();
imagesc(var_same_crop_before_failure);
c=colorbar();
xticklabels((1:4)*5/20);
xlabel('sg');
yticks((1:10)*20-10);
yticklabels({'tF=1','tF=2','tF=3','tF=4'});
ylabel('tf');c.Label.String = 'variance time same cropping pattern before failure';

%dt_diff = diff(dts,1,2);
%[sel, c] = max( dt_diff >.4, [], 2 );
%shift = zeros(10,20);
%cnt = 1;
%for i =1:10
%    for j = 1:20
%        shift(i,j)=c(cnt);
%        cnt = cnt + 1;
%    end
%end
%shift(shift==1) = nan;
%imagesc(shift);colorbar();
%colormap( [0 0 0; parula(255)] );
%shiftT = shift;