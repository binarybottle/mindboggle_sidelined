function q=Qklnu(k,l,nu)

% Computes Q, neccesary constant
% for the moments computation

aux_1=power(-1,k+nu)/power(4,k);
aux_2=sqrt((2*l+4*k+3)/3);
aux_3=trinomial(nu,k-nu,l+nu+1)*nchoosek(2*(l+nu+1+k),l+nu+1+k);
aux_4=nchoosek(2*(l+nu+1),l+nu+1);
q=(aux_1*aux_2*aux_3)/aux_4;