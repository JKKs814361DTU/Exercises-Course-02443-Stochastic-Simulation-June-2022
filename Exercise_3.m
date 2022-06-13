clc;clear; close all
%% Exercise 3
%% 1.
Table(1,:) =["Distribution","T chi^2","p-val chi^2",...
    "T KS","p-val KS"];
Table(2:8,1) = ["Exponential lambda = 2";"Normal Z1 (Box-Mueller)";...
    "Normal Z2 (Box-Mueller)";"Pareto k =2.05";"Pareto k = 2.5";...
    "Pareto k = 3"; "Pareto k = 4"];
N = 1e4; %how many rnd no
U = rand(N,1); %array of rands
%% a. 

lambda = 2; %arbitrary lambda

X_exp = -log(U)/lambda; %Generated varibles
f_exp = @(x) lambda*exp(-lambda*x); %analytical pdf
F_exp = @(x) 1-exp(-lambda*x); %analytical 
pd = makedist('Exponential','mu',1/lambda);
% plot histogram
x = linspace(min(X_exp),max(X_exp),100);
hold on
histogram(X_exp,'Normalization','pdf')
plot(x,f_exp(x),'-r')
%plot(x, pdf(pd, x),'-b');
ylabel('pdf(x)')
xlabel('x')
title("Simulated vs analytical exponential pdf; \lambda = 2")
legend('Simulated','Analytical','MATLAB')
hold off
saveas(gcf,'ex3_1_pdf_exponential.png')
% plot cdf
[f,x] = ecdf(X_exp);
figure;
hold on
plot(x,f,'-*r',x, F_exp(x),'-b')
ylabel('cdf(x)')
xlabel('x')
title("Simulated vs analytical exponential cdf; \lambda = 2")
legend('Simulated','Analytical','MATLAB')
hold off
saveas(gcf,'ex3_1_cdf_exponential.png')
% tests
[~,Table(2,3),stats] = chi2gof(X_exp,'CDF',pd)
Table(2,2) = stats.chi2stat;
[~,Table(2,5),Table(2,4),cv] = kstest(X_exp,'CDF',pd)

%% b.

% Box Muller

U1 = rand(1,N);
U2 = rand(1,N);

Z1 = zeros (1,N);
Z2 = zeros (1,N);

for i=1:N
    Z1(i) = sqrt(-2*log(U1(i)))*cos(2*pi*U2(i));
    Z2(i) = sqrt(-2*log(U1(i)))*sin(2*pi*U2(i));
end

pd=makedist("Normal");
%Z1
% plot histogram
figure;
x = linspace(min(Z1),max(Z1),100);
hold on
histogram(Z1,'Normalization','pdf')
plot(x, pdf(pd, x),'-b');
ylabel('pdf(x)')
xlabel('x')
title("Z1 vs analytical normal pdf;")
legend('Z1','Analytical')
hold off
saveas(gcf,'ex3_1_pdf_Z1.png')
% plot cdf
[f,x] = ecdf(Z1);
figure;
hold on
plot(x,f,'-*r',x, cdf(pd,x),'-b')
ylabel('cdf(x)')
xlabel('x')
title("Z1 vs analytical normal cdf;")
legend('Z1','Analytical')
hold off
saveas(gcf,'ex3_1_cdf_Z1.png')
% tests
[~,Table(3,3),stats] = chi2gof(Z1,'CDF',pd)
Table(3,2) = stats.chi2stat;
[~,Table(3,5),Table(3,4),cv] = kstest(Z1,'CDF',pd)
%Z2
% plot histogram
x = linspace(min(Z2),max(Z2),100);
figure;
hold on
histogram(Z2,'Normalization','pdf')
plot(x, pdf(pd, x),'-b');
ylabel('pdf(x)')
xlabel('x')
title("Z2 vs analytical normal pdf;")
legend('Z2','Analytical')
hold off
saveas(gcf,'ex3_1_pdf_Z2.png')
% plot cdf
[f,x] = ecdf(Z2);
figure;
hold on
plot(x,f,'-*r',x, cdf(pd,x),'-b')
ylabel('cdf(x)')
xlabel('x')
title("Z2 vs analytical normal cdf;")
legend('Z2','Analytical')
hold off
saveas(gcf,'ex3_1_cdf_Z2.png')
% tests
[h,Table(4,3),stats] = chi2gof(Z2,'CDF',pd)
Table(4,2) = stats.chi2stat;
[H,Table(4,5),Table(4,4),cv] = kstest(Z2,'CDF',pd)
close all;
%% c.

beta = 1;
K = [2.05, 2.5, 3,4];
for i =1:4
    k = K(i);

    %Simlated
    X_Pareto = beta*(U.^(-1/k));

    %Analytical
    x = linspace(beta,max(X_Pareto));
    f_Pareto = @(x) k*beta^k./(x.^(k+1));
    F_Pareto = @(x) 1-(beta./x).^k;
    % histogram
    figure;
    hold on
    histogram(X_Pareto,'Normalization','pdf')
    plot(x, f_Pareto(x),'-b');
    ylabel('pdf(x)')
    xlabel('x')
    title("Simulated vs analytical Pareto pdf; k ="+string(k))
    legend('Simulated','Analytical')
    hold off
    saveas(gcf,'ex3_1_pdf_Pareto_k_'+string(k)+'.png')

    % plot cdf
    [f_cdf,x_cdf] = ecdf(X_Pareto);
    figure;
    hold on
    plot(x_cdf,f_cdf,'-*r',x_cdf, F_Pareto(x_cdf),'-b')
    ylabel('cdf(x)')
    xlabel('x')
    title("Simulated vs analytical Pareto cdf;k ="+string(k))
    legend('Simulated','Analytical','Location','best')
    hold off
    saveas(gcf,'ex3_1_cdf_Pareto_k_'+string(k)+'.png')
    
    %KS test
    T_KS = max(f_cdf-F_Pareto(x_cdf));
    p_KS = 1-kolmcdf(T_KS);
    Table(4+i,5) = p_KS;
    Table(4+i,4) = T_KS;

end
writematrix(Table, 'ex3_1.xlsx');

%% 2
clc;clear;close all;
N = 1e4; %how many rnd no
U = rand(N,1); %array of rands
beta = 1;
E_true = @(k) beta.*k./(k-1);
Var_true =@(k) beta^2*k./(k-1).^2/(k-2);
K = linspace(2+1e-9,100);
for i = 1:100
    k =K(i);
    %Simlated
    X_Pareto = beta*(U.^(-1/k));
    E(i) = mean(X_Pareto);
    Var(i) = var(X_Pareto);
end
figure;
hold on
plot(K,abs(E-E_true(K)),'-b')
ylabel('|E_{simulated}-E_{true}|')
grid on
xlabel('k')
%title("Simulated vs analytical Pareto cdf;k ="+string(k))
%legend('Simulated','Analytical','Location','best')
hold off
saveas(gcf,'ex3_2_E.png')
figure;
hold on
plot(K,abs(Var-Var_true(K)),'-b')
ylabel('|Var_{simulated}-Var_{true}|')
grid on
xlabel('k')
%title("Simulated vs analytical Pareto cdf;k ="+string(k))
%legend('Simulated','Analytical','Location','best')
hold off
saveas(gcf,'ex3_2_var.png')
%% 3
clc;clear;close all;
% For the normal distribution generate 100 95% confidence
% intervals for the mean and variance, each based on 10
% observations

N = 10; 

for j=1:2:100
    U1 = rand(1,N);
    U2 = rand(1,N);
    
    Z1 = zeros (1,N);
    Z2 = zeros (1,N);
    
    for i=1:N
        Z1(i) = sqrt(-2*log(U1(i)))*cos(2*pi*U2(i));
        Z2(i) = sqrt(-2*log(U1(i)))*sin(2*pi*U2(i));
    end
    mu(j) = mean(Z1);
    sigma(j) = std(Z1);

    %Calculate confidence intervals
    delta(j) = sigma(j)/sqrt(N) * tinv(1-0.05/2, N-1);

    mu(j+1) = mean(Z2);
    sigma(j+1) = std(Z2);

    %Calculate confidence intervals
    delta(j+1) = sigma(j)/sqrt(N) * tinv(1-0.05/2, N-1);
    
end
figure
hold on
title("100 95% CIs of mean based on 10 obsevarions per CIs")
plot(1:100,mu-delta,'-*r',1:100,mu,'--xb',1:100,mu+delta,'-*r')
legend('CI','mean')
hold off
saveas(gcf,'ex3_3.png')

%CI of variance
[mean(delta)-std(delta)/sqrt(10)*tinv(1-0.05/2, 100-1),...
    mean(delta)+std(delta)/sqrt(10)*tinv(1-0.05/2, 100-1)]