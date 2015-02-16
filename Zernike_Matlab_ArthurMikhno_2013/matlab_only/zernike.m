%% ZERNIKE 3D FROM GEOMETRIC MOMENTS

function Z=zernike(G,N)

% Compute the 3D Zernike moments
% from the geometric moments
V=zeros(N+1,N+1,N+1);
W=zeros(N+1,N+1,N+1);
X=zeros(N+1,N+1,N+1);
Y=zeros(N+1,N+1,N+1);
Z=zeros(N+1,N+1,N+1);

%% FOR a,b,c
% Computing V
i=sqrt(-1);
for a=0:floor(N/2)
    for b=0:(N-2*a)
        for c=0:(N-2*a-b)
            % Calculation for a given a,b,c
            tmp=0;
            for alpha=0:(a+c)
                tmp=tmp+power(i,alpha)*nchoosek(a+c,alpha)*G(2*a+c-alpha+1,alpha+1,b+1);
            end
            V(a+1,b+1,c+1)=tmp;
        end
    end
end

% Computing W
i=sqrt(-1);
for a=0:floor(N/2)
    for b=0:(N-2*a)
        for c=0:(N-2*a-b)
            % Calculation for a given a,b,c
            tmp=0;
            for alpha=0:a
                tmp=tmp+power(-1,alpha)*power(2,a-alpha)*nchoosek(a,alpha)*V(a-alpha+1,b+1,c+2*alpha+1);
            end
            W(a+1,b+1,c+1)=tmp;
        end
    end
end

% Computing X
i=sqrt(-1);
for a=0:floor(N/2)
    for b=0:(N-2*a)
        for c=0:(N-2*a-b)
            % Calculation for a given a,b,c
            tmp=0;
            for alpha=0:a
                tmp=tmp+nchoosek(a,alpha)*W(a-alpha+1,b+2*alpha+1,c+1);
            end
            X(a+1,b+1,c+1)=tmp;
        end
    end
end


%% FOR n,l,m,nu
% Computing Y
for l=0:N
    for nu=0:floor((N-l)/2)
        for m=0:l
            % Calculation for a given l,nu,m
            tmp=0;
            for j=0:floor((l-m)/2)
                tmp=tmp+Yljm(l,j,m)*X(nu+j+1,l-m-2*j+1,m+1);
            end
            Y(l+1,nu+1,m+1)=tmp;
        end
    end
end

%Computing Zernike moments
for n=0:N
    for l=0:n
        if mod(n-l,2)==0
            for m=0:l
                % actual computation
                tmp=0;
                k=(n-l)/2;
                for nu=0:k
                    tmp=tmp+Qklnu(k,l,nu)*conj(Y(l+1,nu+1,m+1));
                end
                Z(n+1,l+1,m+1)=(3/(4*pi))*tmp;
                
            end
        end
    end
end


%correction for improper sign problem in sum([n l m]) even terms
for n=0:N
    for l=0:n
        for m=0:l
            % actual computation
            if mod(sum([n l m]),2) == 0
                Z(n+1,l+1,m+1) = real(Z(n+1,l+1,m+1)) - imag(Z(n+1,l+1,m+1))*i;
            end
            
            if mod(sum([n l m]),2) == 1
                Z(n+1,l+1,m+1) = -real(Z(n+1,l+1,m+1)) + imag(Z(n+1,l+1,m+1))*i;
            end
            
        end
    end
end


































