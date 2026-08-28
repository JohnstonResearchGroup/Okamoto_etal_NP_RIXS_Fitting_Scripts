function out = chi2(P)
load Filename;
load gam;
load sig;
load sig2;
load z;
load sigma_noise;

A = load(Filename);
w = A(:,1);
get_intensities;

I = find(w >= -0.1 & w <= 0.5);

out = sum((A(I,2)./max(A(I,2))+normrnd(0,sigma_noise,size(w(I)))-y_total_LCP(I)).^2) ...
    + sum((A(I,3)./max(A(I,3))+normrnd(0,sigma_noise,size(w(I)))-y_total_RCP(I)).^2);