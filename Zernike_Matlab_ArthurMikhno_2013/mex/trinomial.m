
function t=trinomial(i,j,k)
% Computes the trinomial of
% the three input arguments

aux_1=factorial(i+j+k);
aux_2=factorial(i)*factorial(j)*factorial(k);
t= aux_1/aux_2;

