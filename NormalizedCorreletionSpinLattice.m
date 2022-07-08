function [ MI,Lstat,xi ] = NormalizedCorreletionSpinLattice(spin,Dcut)
    % Function NormalizedCorreletionSpinLattice gives correlation between spins separated by distance d, for d=0,1,2,..,Dcu
    % Correlation is quantified by mutual information (MI).
    % spin: an NxM matrix which specifies the cropping pattern of each site.
    % Sites are connected as a two-dimensional lattice, i.e. each site is
    % connected to four nearest neighbors.
    % Cropping patterns are denoted by nonzero integers, s=1,2,3,...,ns
    % sites with no data are represented by negative integers, e.g. -9999
    % MI: Correlation function
    % xi: Correlation distance
    % Lstat: size of the statistical sample used to calculate the
    % correlation at each distance
    % example: spin=randi(4,10,10); % 10x10 lattice with 4 cropping patterns
    % Dcut=5;
    % [ MI,Lstat,xi ] = NormalizedCorreletionSpinLattice(spin,Dcut);
    sz=size(spin); %dimension of spin
    N=sz(1);
    M=sz(2);
    ns=max(max(spin)); % number of different cropping patterns
    Lstat=zeros(1,Dcut+1);
    MI=zeros(1,Dcut+1);
    Q=zeros(ns,ns,Dcut+1);
    for x1=1:N
        for y1=1:M
            for x2=max(x1-Dcut,1):min(x1+Dcut,N)
                for y2=max(y1-Dcut,1):min(y1+Dcut,M)
                    D=abs(x1-x2)+abs(y1-y2);
                    if D>=0 && D<=Dcut
                        if spin(x1,y1)>0 && spin(x2,y2)>0
                            Q(spin(x1,y1),spin(x2,y2),D+1)=Q(spin(x1,y1),spin(x2,y2),D+1)+1;
                        end
                    end
                end
            end
        end
    end
    SQ=size(Q);
    for dd=1:SQ(3)
        Lstat(dd)=sum(sum(Q(:,:,dd)));
        P=Q(:,:,dd)/Lstat(dd);
        MI(dd)=MutualInformation(P,ns);
    end
    MI=MI/MI(1);
    R=0:Dcut;
    xi=sqrt((R.^2*MI')/sum(MI));
end

function [ MI ] = MutualInformation( P,ns )
    MI=0;
    for rr=1:ns
        for ss=1:ns
            Px=sum(P(rr,:));
            Py=sum(P(ss,:));
            if Px~=0 && Py~=0 && P(rr,ss)~=0
                MI=MI+P(rr,ss)*log2( P(rr,ss)/(Px*Py) );
            end
        end
    end
end
