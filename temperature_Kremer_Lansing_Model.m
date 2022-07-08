function [spins,harvests] = temperature_Kremer_Lansing_Model(N, nrstates, pestradius, harvestradius, temp, nblock, T, a, b, tF, sigma, shock, counter, varargin)
    % This program simulates the evolution of cropping pattern (started from
    % random) and stop at time step T
    % N: dimension of the lattice
    % nrstates: number of cropping patterns
    % pestradius: the spatial extent at which the pests can affect harvests
    % harvestradius: farmers are comparing their harvest to other farmers
    % within this radius, a harvestradius of 1 includes 4 neighbors on the
    % lattice
    % temp: probability to chose state randomly
    % nblock: maximum noise block to be included
    % T: number of timesteps
    % counter: show the time step if counter>1
    % a: pest stress
    % b: water stress
    % example: N=100;
    % nrstates=4;
    % pestradius=2;
    % harvestradius=1;
    % temp=0.05;
    % nblock=4;
    % a=0.5;
    % b=9.6;
    % T=4*N;
    % counter=50;
    % [spin,harvest] = Kremer_Lansing_Model(N,nrstates,pestradius,harvestradius,temp,nblock,T,a,b,counter);
    % figure(1)
    % imagesc(spin)
    %INITIALIZE
    h0=5; % maximal achievable harvest (payoff)
    spins = {};
    harvests = {};
    p = zeros(N,N); % pest load ()
    w = zeros(N,N); % waterstress
    h = zeros(N,N);% harvest
    if ~isempty(varargin)
        s = varargin{1};
    else
        s = randi(nrstates,N); % random states assigned
    end
    s2 = s; % updated states
    t=0;
    % TIME EVOLUTION
    while t<=T        
        bt = b;
        if rand()<shock
            bt = normrnd(b,sigma);            
        end
        if counter>1
            if mod(t,counter)==0
                display(t)
            end
        end
        f=[]; % update fraction of nodes in particular states
        for iz=1:nrstates
            f=[f length(find(s==iz))];
        end
        f=f/N^2;
        for i=1:N
            for j=1:N
                ilimit1=max(1,i-pestradius);
                ilimit2=min(N,i+pestradius);
                SpinNeigh=[]; % vector of state values in neigborhood to compute pest load
                for qq=ilimit1:ilimit2
                    width=pestradius-abs(qq-i);
                    jlimit1=max(1,j-width);
                    jlimit2=min(N,j+width);
                    SpinNeigh=[SpinNeigh s(qq,jlimit1:jlimit2)];
                end
            % update water and pests and harvest
                if isnan(s(i,j))
                    w(i,j) = nan;
                else
                    w(i,j) = f(s(i,j));
                end
                p(i,j) = 1/ ( 0.1+ (length(find(SpinNeigh==s(i,j)))-1) /(length(SpinNeigh)-1) ) ;
                h(i,j) = h0-a*p(i,j)-bt*w(i,j);
                if h(i,j) <0
                    h(i,j) = 0;
                end
                if t>tF
                    failed = 0;
                    for ts = (t-tF):t
                        if harvests{ts}(i,j)<=0
                            failed = failed +1;
                        end
                    end
                    if failed==tF
                        h(i,j) = 0;
                        s(i,j) = nan;
                    end
                end
            end
        end
        % go throuh nodes randomly (not necessary!)
        xs=randperm(N);
        ys=randperm(N);
        for i=xs
            for j=ys
                % check nearest neighbors’ harvests of past timestep
                lowerlimiti=max(1,i-harvestradius);
                upperilimiti=min(N,i+harvestradius);
                HarvestNeigh=[];
                Neigh=[];
                for qq=lowerlimiti:upperilimiti
                    width=harvestradius-abs(qq-i);
                    lowerlimitj=max(1,j-width);
                    upperlimitj=min(N,j+width);
                    HarvestNeigh=[HarvestNeigh h(qq,lowerlimitj:upperlimitj)];
                    Neigh=[Neigh; qq*ones(upperlimitj-lowerlimitj+1,1), (lowerlimitj:upperlimitj)'];
                end
                iii=find(HarvestNeigh==max(HarvestNeigh)); % find neigbor with maximal harvest
                if length(iii)>1, iii=iii(randsample(length(iii),1)); end
                % update state variable (copy most successful neighbor if he is better than you)
                if h(Neigh(iii,1), Neigh(iii,2)) > h(i,j)
                    s2(i,j)= s(Neigh(iii,1), Neigh(iii,2));
                else 
                    s2(i,j)=s(i,j);
                end
                if isnan(s(i,j))
                    s2(i,j) = nan;
                end
            end
        end
        if temp>0
            Nu=0;
            while Nu<temp*N^2
                mu=randsample(N,2,'true');
                LL=randsample(nblock,1);
                if mu(1)+LL-1<=N && mu(2)+LL-1<=N
                    s2(mu(1):mu(1)+LL-1,mu(2):mu(2)+LL-1)=ones(LL,LL)*randi(nrstates,1);
                    Nu=Nu+LL^2;
                end
            end
        end
        s=s2; % update state
        spins{t+1} = s;
        harvests{t+1} = h;
        h2 = h;
        h2(h<0)=0;
        %temp = ginicoeff(reshape(h2,N*N,1))/2;
        t=t+1;
    end
        spin=s;
        harvest=h;
        cooperation=f;
end