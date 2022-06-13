clc; clear;close all;

%% 1. 
% Estimate the integral R01 e^xdx by simulation (the crude Monte Carlo
%estimator). Use eg. an estimator based on 100 samples and present
%the result as the point estimator and a confidence interval.
N = 100;

U = rand(1,N);

X = exp(U);

E_X = mean(X);
Var_X = var(X);
[E_L_X,E_U_X] = CI(X,0.05);

histogram(X)

%% 2. 
% Estimate the integral R01 exdx using antithetic variables, with
% comparable computer ressources
N = 100;
U = rand(1,N);
Y = (exp(U)+exp(1-U))/2;

E_Y = mean(Y);
Var_Y = var(Y);
[E_L_Y,E_U_Y] = CI(Y,0.05);

histogram(Y)

%% 3. 
% Estimate the integral R01 exdx using a control variable, with
%comparable computer ressources.

N = 100;

U = rand(1,N);

X = exp(U);
C = -cov(X,U);
c = C(1,2);
Z = X+c*(U-1/2);

E_Z = mean(Z);
Var_Z = var(Z);
[E_L_Z,E_U_Z] = CI(Z,0.05);

histogram(Z)

%% 4. 
% Estimate the integral R01 exdx using stratified sampling, with
% comparable computer ressources.
N = 100;
U = rand(sqrt(N),sqrt(N));

W = zeros(sqrt(N),1);
for i = 1:sqrt(N)
    for j = 1:sqrt(N)
    W(i) = W(i) + exp((U(i,j)+(j-1))/10);
    end
    W(i) = W(i)/10;
end

E_W = mean(W);
Var_W = var(W);
[E_L_W,E_U_W] = CI(W,0.05);

histogram(W)

Table(1,:) =["Method","E(X)","Var(X)","CI_{lower}","CI_{upper}"];
Table(1,:) = ["Analytical",exp(1)-1,"NA","NA","NA"];
Table(3,:) = ["Curde Monte Carlo",E_X,Var_X,E_L_X,E_U_X];
Table(4,:) = ["Antithetic Variables",E_Y,Var_Y,E_L_Y,E_U_Y];
Table(5,:) = ["Control Variable",E_Z,Var_Z,E_L_Z,E_U_Z];
Table(6,:) = ["Stratified sampling",E_W,Var_W,E_L_W,E_U_W];


errorbar(categorical(Table(3:6,1)),[E_X,E_Y,E_Z,E_W], ...
    [E_L_X,E_L_Y,E_L_Z,E_L_W]-[E_X,E_Y,E_Z,E_W], ...
    [E_U_X,E_U_Y,E_U_Z,E_U_W]-[E_X,E_Y,E_Z,E_W],'o')
ylabel('E(X)')
yline(exp(1)-1,'-','exp(1)-1');
saveas(gcf,'ex5_1.png')
writematrix(Table, 'ex5_1.xlsx');
%% 7. 
% For a standard normal random variable Z ∼ N(0, 1) using the crude
%Monte Carlo estimator estimate the probability X > a. Then try
%importance sampling with a normal density with mean a and
% variance σ2. For the expirements start using σ2 = 1, use different
% values of a (e.g. 2 and 4), and different sample sizes. If time
% permits experiment with other values for σ2. Finally discuss the
% efficiency of the methods
clc; clear;close all;

N = 1e2;
Table(1,:) =["Method","a","sigma","P(X>a)"];

for a = [2,4]
    sigma_g = 1;
    P_true = cdf('Normal',a,0,1,'upper');
    Table(end+1,2:4) = [a,sigma_g,P_true];
    X = random('Normal',0,1,1,N);
    
    P_crude = mean(X>a)
    Table(end+1,2:4) = [a,sigma_g,P_crude];
    %importance sampling
    
    Y = random('Normal',a,sigma_g,1,N);
    h = Y>a;
    f = pdf('Normal',Y,0,1);
    g = pdf('Normal',Y,a,sigma_g);
    
    P_importance = mean(h.*f./g)
    
    
    Table(end+1,2:4) = [a,sigma_g,P_importance];

end
Table(2:7,1) = ["Analytical";"Crude Monte Carlo"; "Importance sampling";...
    "Analytical";"Crude Monte Carlo"; "Importance sampling"]
writematrix(Table, 'ex5_7.xlsx');
%% 8. 
% Use importance sampling with g(x) = λ exp (−λ ∗ x) to calculate
% the integral R01 exdx . Try to find the optimal value of λ by
% calculating the variance of h(x)f(x)/g(x) and verify by simulation.
% Note that importance sampling with the exponential distribution
% will not reduce the variance.

%importance sampling

clc;clear;
%lambda = fminsearch(@V,1e3)
lambda = 1.34;
N = 100;
mu = 1/abs(lambda);
Y = random('Exponential',mu,1,N);

h = 0<=Y & Y<=1;
f = exp(Y);

g = pdf('Exponential',Y,mu);

P_importance = mean(h.*f./g)
[E_L,E_U] =CI(h.*f./g,0.05)
P_true = exp(1)-1;
v = var(h.*f./g)
function v = V(lambda)
N = 100;
mu = 1/abs(lambda);
Y = random('Exponential',mu,1,N);

h = 0<=Y & Y<=1;
f = exp(Y);

g = pdf('Exponential',Y,mu);

P_importance = mean(h.*f./g);
[E_L,E_U] =CI(h.*f./g,0.05);
P_true = exp(1)-1;
v = var(h.*f./g);
end


% Table(7,:) = ["Importance sampling",P_importance,var(h.*f./g),E_L,E_U];
% 
% errorbar(categorical(Table(3:7,1)),[E_X,E_Y,E_Z,E_W,P_importance], ...
%     [E_L_X,E_L_Y,E_L_Z,E_L_W,E_L]-[E_X,E_Y,E_Z,E_W,P_importance], ...
%     [E_U_X,E_U_Y,E_U_Z,E_U_W,E_U]-[E_X,E_Y,E_Z,E_W,P_importance],'o')
% ylabel('E(X)')
% yline(exp(1)-1,'-','exp(1)-1');
% saveas(gcf,'ex5_8.png')
% writematrix(Table, 'ex5_8.xlsx');