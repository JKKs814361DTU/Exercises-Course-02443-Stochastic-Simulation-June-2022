%% Exercise 2 Discrete random variables
%In the excercise you can use a build in procedure for generating
%random numbers. Compare the results obtained in simulations with
%expected results. Use histograms (and tests).
clc;clear;close all;

%% 1.
%Choose a value for the probability parameter p in the geometric
%distribution and simulate 10,000 outcomes. You can experiment
%with a small, moderate and large value if you like.

N = 10000; %number of outcomes
U = rand(N,1); %uniform random numbers
Table =["p","T_KS","p_KS"];
for p = [1e-3,0.45,0.999]; %probability


% Generate geometric distribution
f_geom = @(n) (1-p).^(n-1).*p; %True pdf
F_geom = @(n) 1-(1-p).^n; %True cdf
X_geom = floor(log(U)./log(1-p))+1; %Generated geometric
%%% plot histogram
x = min(X_geom):max(X_geom);
figure;
hold on
h = histogram(X_geom,'BinMethod','integers');
plot(x,f_geom(x)*N,'*r')
legend('generated','true')
title("Geometric distribution with p = "+string(p))
xlabel("n")
ylabel("Occurance")
hold off
filename = 'ex2_p_'+string(p)+'histogram_geom.png';
saveas(gcf,filename)


%chi2 test
n_expected = f_geom(x)*N;
n_observed = h.BinCounts;
n_observed =histcounts(X_geom,length(x));


%Two-sample Kolmogorov-Smirnov test
X_true_geom = [];
for i = 1:length(n_expected)
    X_true_geom = [X_true_geom; i*ones(round(n_expected(i)),1)];
end
[h,p_KS,T_KS] = kstest2(X_geom,X_true_geom)

Table = [Table;[p,T_KS,p_KS]];
end
writematrix(Table, 'ex2.xlsx');

%% 2.
clc;clear;close all;
%Simulate the 6 point distribution with
p = [7/48,5/48,1/8,1/16,1/4,5/16];
k = length(p);
%Table for results
Table(1,:) =["Method","Generated rnd","Execution time","T chi^2","p-val chi^2",...
    "T KS","p-val KS"];
Table(2:7,1) = ["Direct crude";"Rejection";"Alias";"Direct crude";"Rejection";"Alias"];
for w = 1:2
NN = [1e2,1e6]
N = NN(w); %number of outcomes


% Outcomes for KS test
n_expected = p*N;
X_true = [];
for i = 1:length(n_expected)
    X_true = [X_true; i*ones(round(n_expected(i)),1)];
end


%% a. Direct (crude) method

U = rand(N,1); %uniform random numbers
F = [];
X_crude = [];
for i = 1:length(p)
    F(i) = sum(p(1:i));
end
F = [0,F];
tic;
for i = 1:N
    for xj = 1:(k)
        if U(i)> F(xj)
            if U(i) <= F(xj+1)
                X_crude(i) = xj;        
                break
        X_crude(i) = k;
            end
        end
    end
end
Table(2+(w-1)*3,3) = toc; 
[h1,Table(2+(w-1)*3,4), Table(2+(w-1)*3,5),Table(2+(w-1)*3,6),Table(2+(w-1)*3,7)] ...
    = plot_test(X_crude,X_true,n_expected,k,"Crude")

Table(2+(w-1)*3,2) = N;
%% rejection method
tic;
X_rejection = [];
c = max(p);
N_gen_rnd = 0;

for i = 1:N
    U2 = 2; %ensure that one iteration runs
    I = 1;
    while U2 > p(I)/c
        U1 = rand();
        U2 = rand();
        N_gen_rnd = N_gen_rnd+2;
        I = floor(U1*k)+1;
    end
    X_rejection(i) = I;  
end
Table(3+(w-1)*3,3) = toc;

Table(3+(w-1)*3,2) = N_gen_rnd;
[h2,Table(3+(w-1)*3,4), Table(3+(w-1)*3,5),Table(3+(w-1)*3,6),Table(3+(w-1)*3,7)]...
    = plot_test(X_rejection,X_true,n_expected,k,"Rejection")
%% Alias method
%Generate look up tables
tic;
eps = 1e-8;
L = 1:k;
F = k*p;
G = find(F>=1); S = find(F<=1);

while ~isempty(S)
    i = G(1); j=S(1); 
    L(j) = i; F(i) = F(i)-(1-F(j)); 
    if F(i)<1-eps
        G(1) =[];
        S = [S i];
    end
    S(1) = [];

end

%Generation

X_Alias = [];

for i = 1:N
    U1 = rand();
    U2 = rand();
    I = floor(U1*k)+1;
    if U2 <= F(I)
        X_Alias(i) = I;
    else
        X_Alias(i) =L(I);
    end
end
Table(4+(w-1)*3,3) = toc;
Table(4+(w-1)*3,2) = 2*N;
[h3,Table(4+(w-1)*3,4), Table(4+(w-1)*3,5),Table(4+(w-1)*3,6),Table(4+(w-1)*3,7)]...
    = plot_test(X_rejection,X_true,n_expected,k,"Alias")
end
writematrix(Table, 'ex2_2.xlsx');
%% Ploting and tests
function [h, T_chi,pval_chi,T_KS,p_KS] = plot_test(X,X_true,n_expected,k,method)
    %%% plot histogram
    x = 1:6;
    figure;
    hold on
    h = histogram(X,'BinMethod','integers');
    plot(x,n_expected,'*r')
    legend('generated','true',Location='northwest')
    title(method+" method")
    xlabel("n")
    ylabel("Occurance")
    hold off
    filename = 'ex2_N_'+string(length(X_true))+'histogram_'+method+'.png';
    saveas(gcf,filename)
    
    %chi2 test
    
    n_observed = h.BinCounts;
    n_classes = k;
    [pval_chi, T_chi] = chi2(n_observed,n_expected,n_classes);
    
    %Two-sample Kolmogorov-Smirnov test
    
    [H,p_KS,T_KS] = kstest2(X,X_true);
end 