function out = chi2(P)
load Filename;
load gam;
load sig;
load sig2;
load z;

A = load(Filename);
w = A(:,1);
get_intensities;

I = find(w >= -0.1 & w <= 0.5);

boundary = 0; 
% for nn = 1:length(P)
%   if(P(nn) < 0) boundary = boundary + 10000; end
% end
out = sum((A(I,2)./max(A(I,2))-y_total_LCP(I)).^2) ...
    + sum((A(I,3)./max(A(I,3))-y_total_RCP(I)).^2) ...
    + boundary;

