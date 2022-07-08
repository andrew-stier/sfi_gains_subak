function [S,T] = PatchSize(spin)
    % Function PatchSize gives size (S) and cropping pattern (T) of each patch.
    % spin: an NxM matrix which specifies the cropping pattern of each site.
    % Sites are connected as a two-dimensional lattice, i.e. each site is
    % connected to four nearest neighbors.
    % Cropping patterns are denoted by nonzero integers,
    % sites with no data are represented by negative integers
    % example: spin=randi(4,5,5); % 5x5 lattice with 4 cropping patterns
    % [S,T]=PatchSize(spin);
    sz=size(spin); % dimension of spin (NxM)
    N=sz(1);
    M=sz(2);
    k=4; % each site is connected to four nearest neighbors
    List=zeros(N*M,k); % Adjacency list
    for jj=0:N-1
        for kk=0:M-1
            v=(kk)*N+jj+1;
            link=1;
            if jj-1>-1
                vN=(kk)*N+(jj-1)+1;
                List(v,link)=vN;
                link=link+1;
            end
            if jj+1<N
                vS=(kk)*N+(jj+1)+1;
                List(v,link)=vS;
                link=link+1;
            end
            if kk+1<M
                vE=(kk+1)*N+(jj)+1;
                List(v,link)=vE;
                link=link+1;
            end
            if kk-1>-1
                vW=(kk-1)*N+(jj)+1;
                List(v,link)=vW;
                link=link+1;
            end
        end
    end
    degree=zeros(N*M,1); % degree of the lattice
    for ii=1:N*M
        degree(ii)=length(find(List(ii,:)~=0));
    end
    cluster=zeros(1,N*M); % cluster which each site belongs to
    Nc=1; % index of cluster
    T=zeros(N,1);
    % assignes cluster index to each site
    for nn=1:N*M
        Origin=nn;
        if cluster(Origin)==0
            crop=spin(mod(Origin-1,N)+1,ceil(Origin/N));
            if crop>0
                cluster(Origin)=Nc;
                T(Nc)=crop;
                Queue=zeros(1,N*M);
                Queue(1)=Origin;
                q=1;
                uu=1;
                while Queue(uu)>0
                    current=Queue(uu);
                    Connected=List(current,1:degree(current));
                    for j=1:length(Connected)
                        if spin(mod(Connected(j)-1,N)+1,ceil(Connected(j)/N))==crop && cluster(Connected(j))==0
                            cluster(Connected(j))=Nc;
                            Queue(q+1)=Connected(j);
                            q=q+1;
                        end
                    end
                    uu=uu+1;
                end
                Nc=Nc+1;
            end
        end
    end
    T=T(1:max(cluster));
    S=zeros(max(cluster),1);
    for cc=1:max(cluster)
        S(cc)=length(find(cluster==cc));
    end
end