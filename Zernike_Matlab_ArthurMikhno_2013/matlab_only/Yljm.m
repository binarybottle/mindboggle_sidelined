function y=Yljm(l,j,m)

% Computes a constant used in Zernike moments

aux_1=power(-1,j)*(sqrt(2*l+1)/power(2,l));
aux_2=trinomial(m,j,l-m-2*j)*nchoosek(2*(l-j),l-j);
aux_3=sqrt(trinomial(m,m,l-m));
y=(aux_1*aux_2)/aux_3;

