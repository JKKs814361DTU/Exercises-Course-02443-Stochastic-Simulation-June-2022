function Ui = U_ecd(A,S)
n = length(S);
S = [S,S(1)];
Ui=0;
for i = 1:n
    Ui = Ui+A(S(i),S(i+1));
end