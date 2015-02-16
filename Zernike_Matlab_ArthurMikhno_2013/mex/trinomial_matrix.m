
function tri_matrix=trinomial_matrix(N)
% Computes the trinomial of
% the three input arguments

tri_matrix=zeros(N+1,N+1,N+1);
for i=0:N
    for j=0:(N-i)
        for k=0:(N-i-j)
            tri_matrix(i+1,j+1,k+1)=trinomial(i,j,k);
        end
    end
end