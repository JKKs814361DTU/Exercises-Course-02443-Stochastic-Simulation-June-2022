function [pval, T] = chi2(n_observed,n_expected,n_classes)
    T = sum((n_observed - n_expected).^2./n_expected);
    df = n_classes - 1;%-3;
    pval = chi2cdf(T,df,'upper');
end