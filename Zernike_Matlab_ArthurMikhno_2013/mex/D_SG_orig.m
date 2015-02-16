function G = D_SG_orig(num_facets,i,N,C,D,Vol,F)

G=zeros(N+1,N+1);
S=zeros(num_facets,N+1,N+1,N+1);

for j=0:(N-i)
    for k=0:(N-i-j)
        
        S = D_SG_orig_part(num_facets,i,j,k,C,D,S,Vol,F);
        % Geometric moments after summing through all facets
        G(j+1,k+1)=sum(S(:,i+1,j+1,k+1));
    end
end