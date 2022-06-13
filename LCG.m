function X = LCG(x0,a,c,M,n)

if gcd(M,c) ~= 1
    warning('M & c are not relative prime. Cycle length is not maximum.')
end

for p = factor(16)
    if mod(a,p) ~= 1
        warning('mod(a,p) ~= 1. Cycle length is not maximum.')
        break
    end
end

   
X(1) = x0;
    for i = 1:n
    X(i+1) = mod( a*X(i) + c, M);
    end
X = X(2:end);


end
