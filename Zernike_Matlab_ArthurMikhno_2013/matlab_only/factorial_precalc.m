function F = factorial_precalc(N)

for i = 1:N
    
    for j=0:(N-i)
        for k=0:(N-i-j)
            
            F(i + 1) = factorial(i);
            F(j + 1) = factorial(j);
            F(k + 1) = factorial(k);
            F(i+j+k+2 + 1) = factorial(i+j+k+2);
            
        end
    end
end
