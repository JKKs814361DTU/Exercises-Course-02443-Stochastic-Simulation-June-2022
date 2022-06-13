% Exercise 1
clc;close all; clear;
rng(0)
%% 1 Linear Congruential Generator of Random Numbers

%% a) Random numbers generation

N = 10000;

X = zeros(1, N);
U = zeros(1, N);
X(1) = 3;

% Parameters
a = 7^4*2; % 0 < a < M
c = 1; % 0 <= c < M
M = 3^31-1; 
Mstr = '65536'
a = 129;
c = 27061;
M = 65536;
a = 5;
c = 12;
M = 13;
Mstr = '13'
filename = 'Exercise_1_1_M  ='+string(M)+'a='+string(a)+', c= ' +string(c);
for i = 2:N
    X(i) = mod(a*X(i-1)+c, M);
    U(i) = X(i)/M;
end

U_1 = [3/M,U(1:(end-1))];

%% b)  Evaluate the quality of the generator

classes = 10;
Count = histcounts(U,10);
Expected = N/classes;

%histogram
%subplot(1,2,1)
figure;
h = histogram(U,classes);
ylabel('Occurances')
xlabel('Generated classes of X')
title('LCG. Parameters: M  =' +string(Mstr)+', a='+string(a)+', c= ' ...
    +string(c));
saveas(gcf,'histogram_'+filename,'png')
n_observed = h.BinCounts;
%%
%subplot(1,2,2)
figure;
% Scatter plot
%subplot(1,2,2)
plot(U_1,U,'*')
t = 'LCG.Random numbers U_i against U_{i+1}, X_{i+1} = ('+string(a)+...
    ' X_i + '+string(c)+')(mod '+string(Mstr)+')';
title(t);
ylabel('U_{i+1}')
xlabel('U_i')
hold off

saveas(gcf,'scatter_'+filename,'png')
% Chi square Test

Count = histcounts(U,10);
Expected = N/classes;
Tchi = 0;
estparameters = 3;
dof = classes-1-estparameters;

for i=1:classes
    
    T = ((Count(i)-Expected)^2)/Expected;
    Tchi = Tchi + T;
    
end

pvalueChi2 = chi2cdf(Tchi,dof,'upper');%it should be either 1-chi2cdf 
                                        % or with attribute 'upper'

% KS Test
% The Kolmogorov-Smirnov test is used to decide if a sample 
% comes from a population with a specific distribution.

xx = linspace(0,1,N);

xs = sort(U);
Fe = zeros(1,N);
for i=1:N
    Fe(i) = sum(xs(i)>= U)/N;
end

pd = makedist('Uniform');
Fx = cdf(pd,xx);

figure;
plot(xs,Fe,'-*',xx,Fx,'--');
title('Emperical cdf F_e(x) vs. true uniform F_{uniform}(x) for N ='...
    +string(N)+'random numbers');
legend('F_e(x)','F_{uniform}(x)','Location','northwest');
saveas(gcf,'KScdf_'+filename,'png')
D = max((Fe-Fx));


% With all parameters known, the adjusted test statistic is

Dadj = (sqrt(N)+0.12+(0.11/sqrt(N)))*D;

% and the p-value is

if Dadj < 1.138
    pvalueKS = 1;
elseif Dadj < 1.224
    pvalueKS = 0.15;
elseif Dadj < 1.358
    pvalueKS = 0.10;
elseif Dadj < 1.480
    pvalueKS = 0.05;
elseif Dadj < 1.628
    pvalueKS = 0.01;
else
    pvalueKS = 0;
end

p_KS = 1-kolmcdf(D);
% Run test I - ABOVE/BELOW

% Wald-Wolfowitz Test (also called Wald-Wolfowitz run test) is 
% a non-parametric 
% hypothesis test used to test the randomness of a two-valued data sequence. 
% It tests to see if the sequence are mutually independent.


Ume = median(U);
Below = U<Ume;
Above = U>Ume;
n1 = sum(Above(:)==1);
n2 = sum(Below(:)==1);


props = regionprops(Above, 'Area')
ra = [props.Area];
Ra = length(ra)
props = regionprops(Below, 'Area')
rb = [props.Area];
Rb = length(rb);


WWt = Ra + Rb; %Wald-Wolfowitz test statistic

mu = 2*n1*n2/(n1+n2);
sig = sqrt(2*(n1*n2*(2*n1*n2-n1-n2))/((n1+n2)^(2)*(n1+n2-1)));
Ru = mu+1.96*sig; % At 5% level of significance
Rl = mu-1.96*sig; % At 5% level of significance

% If Rl<WWt<Ru H0 is accepted and the sequence is random
RunTestI = Rl<WWt & WWt<Ru;
p_RT1 = normcdf(WWt,mu,sig,'upper');
% Run test II - UP/DOWN from Knuth
r = [ra,rb];

for i = 1:5
    R(i,1) = sum(r==i);
end
R(6,1) = sum(r>=6);
A = [4529.4 9044.9 13568 18091 22615 27892;...
    9044.9 18097 27139 36187 45234 55789;...
    13568 27139 40721 54281 67852 83685;... 
    18091 36187 54281 72414 90470 111580;...
    22615 45234 67852 90470 113262 139476;...
    27892 55789 83685 111580 139476 172860];
B = [1/6 5/24 11/120 19/720 29/5040 1/840]';

V = (1/(N-6))*(R-N*B)'*A*(R-N*B); % To compare with Chi2(6) distribution
p_RT2 = chi2cdf(V,6,'upper');

% Run test III 

Rt3 = []; % count runs

for i=1:(length(U)-1)
    if U(i)<U(i+1)
        Rt3 = [Rt3 1]; %increase
    else
        Rt3 = [Rt3 0]; %decrease
    end
end
Rt3 = logical(Rt3);
props = regionprops(Rt3, 'Area')
rl = [props.Area];

props = regionprops(~Rt3, 'Area')
rl = [rl,[props.Area]];
X =length(rl);

Z=(X-(2*N-1)/(3))/(sqrt((16*N-29)/(90)));
p_RT3 = normcdf(Z,0,1,'upper');

%Correlation
c = zeros(N,1);
for h = 1:N
    for i = 1:(N-h)
        c(h)=c(h)+(1)/(N-h)*U(i)*U(i+h);
    end
end
pd = makedist("Normal",0.25,7/144/N);
[h,p_correl,T_correl] = chi2gof(c,"CDF",pd)

% Assuming the numbers are random, the correlation between pairs at any
% distance should be independent of the distance (k) and ≈0.25
table =...
    [filename,Tchi,pvalueChi2,Dadj,p_KS,WWt,p_RT1,V,p_RT2,Z,p_RT3,p_correl];
% [numbers, strings, raw] = xlsread('ex1.xlsx');
% lastRow = size(raw, 1);
% nextRow = lastRow + 1;
% cellReference = sprintf('A%d', lastRow);
% xlswrite('ex1.xlsx', [raw;numbers;table], 'Sheet1', cellReference);
writematrix(table, 'ex1.xlsx', 'WriteMode','append');
%% c) Experimenting with different values of 'a', 'c', 'M'

% From "Simulation, Sheldon M. Ross, Elsevier, 2013"
%
% Since each of the numbers xn assumes one
% of the values 0, 1, . . . ,m − 1, it follows that after some finite number (of at most
% m) of generated values a value must repeat itself; and once this happens the whole
% sequence will begin to repeat. Thus, we want to choose the constants a and m so
% that, for any initial seed x0, the number of variables that can be generated before
% this repetition occurs is large.
% In general the constants a and m should be chosen to satisfy three criteria:
% 1. For any initial seed, the resultant sequence has the “appearance” of being a
% sequence of independent uniform (0, 1) random variables.
% 2. For any initial seed, the number of variables that can be generated before
% repetition begins is large.
% 3. The values can be computed efficiently on a digital computer.
% A guideline that appears to be of help in satisfying the above three conditions
% is that m should be chosen to be a large prime number that can be fitted to the
% computer word size

% An example of good parameters is:
% a = 7^5; 
% c = 1; 
% M = 2^31-1; 

% An example of bad parameters is:
% a = 5;
% c = 12;
% M = 13;



%% 2. Apply a system available generator and perform the various
%statistical tests you did under Part 1 point (b) for this
%generator too.
clc;clear;close all;rng(0)
N=10000;

U = rand(1,N);
U_1 = [0,U(1:(end-1))];
filename = string('Exercise_1_1_Matlab');

classes = 10;
Count = histcounts(U,10);
Expected = N/classes;

%histogram
%subplot(1,2,1)
figure;
h = histogram(U,classes);
ylabel('Occurances')
xlabel('Generated classes of U')
title('Matlab inbuilt generator');
saveas(gcf,'histogram_'+filename,'png')
n_observed = h.BinCounts;
%%
%subplot(1,2,2)
figure;
% Scatter plot
%subplot(1,2,2)
plot(U_1,U,'*')
t = 'Matlab generated U_i against U_{i+1}';
title(t);
ylabel('U_{i+1}')
xlabel('U_i')
hold off

saveas(gcf,'scatter_'+filename,'png')
% Chi square Test

Count = histcounts(U,10);
Expected = N/classes;
Tchi = 0;
estparameters = 0;
dof = classes-1-estparameters;

for i=1:classes
    
    T = ((Count(i)-Expected)^2)/Expected;
    Tchi = Tchi + T;
    
end

pvalueChi2 = chi2cdf(Tchi,dof,'upper');%it should be either 1-chi2cdf or 
% with attribute 'upper'

% KS Test
% The Kolmogorov-Smirnov test is used to decide if a sample 
% comes from a population with a specific distribution.

xx = linspace(0,1,N);

xs = sort(U);
Fe = zeros(1,N);
for i=1:N
    Fe(i) = sum(xs(i)>= U)/N;
end

pd = makedist('Uniform');
Fx = cdf(pd,xx);

figure;
plot(xs,Fe,'-*',xs,Fx,'--');
title('Emperical cdf F_e(x) vs. true uniform F_{uniform}(x) for N ='...
    +string(N)+'random numbers');
legend('F_e(x)','F_{uniform}(x)','Location','northwest');
saveas(gcf,'KScdf_'+filename,'png')
D = max((Fe-Fx));


% With all parameters known, the adjusted test statistic is

Dadj = (sqrt(N)+0.12+(0.11/sqrt(N)))*D;

% and the p-value is

if Dadj < 1.138
    pvalueKS = 1;
elseif Dadj < 1.224
    pvalueKS = 0.15;
elseif Dadj < 1.358
    pvalueKS = 0.10;
elseif Dadj < 1.480
    pvalueKS = 0.05;
elseif Dadj < 1.628
    pvalueKS = 0.01;
else
    pvalueKS = 0;
end

p_KS = 1-kolmcdf(Dadj);
% Run test I - ABOVE/BELOW

% Wald-Wolfowitz Test (also called Wald-Wolfowitz run test) is a non-parametric 
% hypothesis test used to test the randomness of a two-valued data sequence. 
% It tests to see if the sequence are mutually independent.


Ume = median(U);
Below = U<Ume;
Above = U>Ume;
n1 = sum(Above(:)==1);
n2 = sum(Below(:)==1);


props = regionprops(Above, 'Area')
ra = [props.Area];
Ra = length(ra)
props = regionprops(Below, 'Area')
rb = [props.Area];
Rb = length(rb);


WWt = Ra + Rb; %Wald-Wolfowitz test statistic

mu = 2*n1*n2/(n1+n2);
sig = sqrt(2*(n1*n2*(2*n1*n2-n1-n2))/((n1+n2)^(2)*(n1+n2-1)));
Ru = mu+1.96*sig; % At 5% level of significance
Rl = mu-1.96*sig; % At 5% level of significance

% If Rl<WWt<Ru H0 is accepted and the sequence is random
RunTestI = Rl<WWt & WWt<Ru;
p_RT1 = normcdf(WWt,mu,sig,'upper');
% Run test II - UP/DOWN from Knuth
r = [ra,rb];

for i = 1:5
    R(i,1) = sum(r==i);
end
R(6,1) = sum(r>=6);
A = [4529.4 9044.9 13568 18091 22615 27892;...
    9044.9 18097 27139 36187 45234 55789;...
    13568 27139 40721 54281 67852 83685;... 
    18091 36187 54281 72414 90470 111580;...
    22615 45234 67852 90470 113262 139476;...
    27892 55789 83685 111580 139476 172860];
B = [1/6 5/24 11/120 19/720 29/5040 1/840]';

V = (1/(N-6))*(R-N*B)'*A*(R-N*B); % To compare with Chi2(6) distribution
p_RT2 = chi2cdf(V,6,'upper');

% Run test III 

Rt3 = []; % count runs

for i=1:(length(U)-1)
    if U(i)<U(i+1)
        Rt3 = [Rt3 1]; %increase
    else
        Rt3 = [Rt3 0]; %decrease
    end
end
Rt3 = logical(Rt3);
props = regionprops(Rt3, 'Area')
rl = [props.Area];

props = regionprops(~Rt3, 'Area')
rl = [rl,[props.Area]];
X =length(rl);

Z=(X-(2*N-1)/(3))/(sqrt((16*N-29)/(90)));
p_RT3 = normcdf(Z,0,1,'upper');


% Assuming the numbers are random, the correlation between pairs at any
% distance should be independent of the distance (k) and ≈0.25


%Correlation
c = zeros(N,1);
for h = 1:N
    for i = 1:(N-h)
        c(h)=c(h)+(1)/(N-h)*U(i)*U(i+h);
    end
end
pd = makedist("Normal",0.25,(7/144/N));
[h,p_correl] = chi2gof(c,"CDF",pd)

table = [filename,Tchi,pvalueChi2,Dadj,p_KS,WWt,p_RT1,V,p_RT2,Z,p_RT3,p_correl];
writematrix(table, 'ex1.xlsx', 'WriteMode','append');

