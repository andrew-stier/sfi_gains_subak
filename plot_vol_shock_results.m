load("Climate_shock_run_kl_model_3.mat")
figure();
imagesc(dts(1:end,1:400));
c=colorbar();
xlabel('timestep');
yticks((1:10)*20-10);
yticklabels({'tF=1','tF=2','tF=3','tF=4','tF=5','tF=6','tF=7','tF=8','tF=9','tF=10'});
ylabel('time to farm failure (tF)');c.Label.String = 'fraction of farms failed';


figure();
imagesc(gts);
c=colorbar();
xlabel('timestep');
yticks((1:10)*20-10);
yticklabels({'tF=1','tF=2','tF=3','tF=4','tF=5','tF=6','tF=7','tF=8','tF=9','tF=10'});
ylabel('time to farm failure (tF)');c.Label.String = 'gini coefficient';

figure();
imagesc(hts);
c=colorbar();
xlabel('timestep');
yticks((1:10)*20-10);
yticklabels({'tF=1','tF=2','tF=3','tF=4','tF=5','tF=6','tF=7','tF=8','tF=9','tF=10'});
ylabel('time to farm failure (tF)');
c.Label.String = 'harvest amount';

figure();
imagesc(xis);
c=colorbar();
xlabel('volatility (% of b)');
xticklabels((1:4)*5/20);
ylabel('time to farm failure');
c.Label.String = 'xi';

dt_diff = diff(dts,1,2);
[sel, c] = max( dt_diff >.4, [], 2 );
shift = zeros(10,20);
cnt = 1;
for i =1:10
    for j = 1:20
        shift(i,j)=c(cnt);
        cnt = cnt + 1;
    end
end
shift(shift==1) = nan;
imagesc(shift);colorbar();
colormap( [0 0 0; parula(255)] );
shiftT = shift;