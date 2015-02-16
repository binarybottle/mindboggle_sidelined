function c=mon_comb(tri_matrix,vertex,N)
% Computes the value for a specified monomial
% combination for a given vertex

% Vertex coordinates
x=vertex(1);
y=vertex(2);
z=vertex(3);
% Computation
c=zeros(N+1,N+1,N+1);
for i=0:N
    for j=0:(N-i)
        for k=0:(N-i-j)
            mon=power(x,i)*power(y,j)*power(z,k);
            c(i+1,j+1,k+1)=tri_matrix(i+1,j+1,k+1)*mon;
        end
    end
end
