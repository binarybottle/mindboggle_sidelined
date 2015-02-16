
function Descriptors=feature_extraction(Z,N)

% The features are the magnitude of the Zernike moments.
% Collect the moments into (2*l+1)-dimensional vectors 
% and define the rotationally invariant 3D Zernike descriptors 
% Fnl as norms of vectors Znl.

F=zeros(N+1,N+1)+NaN;
% Calcuate the magnitude
for n=0:N
    for l=0:n
        if mod(n-l,2)==0
            aux_1=Z(n+1,l+1,1:l+1);
            aux_1=aux_1(:)';
            if l>0
                % Taking the values for m<0
                aux_2=conj(aux_1(2:(l+1)));
                for m=1:l
                    aux_2(m)=aux_2(m)*power(-1,m);
                end
                % Sorting the vector
                aux_2=rot90(aux_2,2);
                aux_1=cat(2,aux_2,aux_1);
            end
            F(n+1,l+1)=norm(aux_1,2);
        end
    end
end
% Generate vector of Zernike descriptors
% Column vector that store the magnitude
% of the Zernike moments for n,l.
idx=find(F>=0);
Descriptors=F(idx);