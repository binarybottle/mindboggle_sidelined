function D=Dabc_orig(C,N)
% Pre-computes the values for D
% and store them in matrix form
% to be used for Geometric Moments

D = zeros(N+1,N+1,N+1);
for a=0:N
    for b=0:N
        for c=0:N
            temp=0;
            for i2=0:a
                for j2=0:b
                    for k2=0:c
                        temp=temp+C(2,i2+1,j2+1,k2+1)*C(3,a-i2+1,b-j2+1,c-k2+1);
                    end
                end
            end          
            D(a+1,b+1,c+1)=temp;
        end
    end
end