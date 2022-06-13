function [theta_L,theta_U] = CI(Theta,alpha)

theta = mean(Theta);
n = length(Theta);

theta_L = theta+std(Theta)/sqrt(n)*tinv(alpha/2,n-1);

theta_U = theta+std(Theta)/sqrt(n)*tinv(1-alpha/2,n-1);


end