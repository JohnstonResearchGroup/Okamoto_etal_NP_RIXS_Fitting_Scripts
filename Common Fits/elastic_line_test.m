clc;
close all;
clear;

sig = 0.009;
res = 0.024;

disp(sqrt(sig*sig+res*res))

Filename = '../Data/B69Kdata_0803.txt';
A = load(Filename);
w = A(:,1);
LCP = A(:,2);
RCP = A(:,3);

I = find(w > -0.1 & w < 0.025);

yel = lorentzian(w,0,sig);
I_el = zeros(size(w)); 

for nw = 1:1:length(w)
   gauss = gaussian(w,w(nw),res/2.355);
   I_el(nw) = trapz(w,-yel.*gauss);
end

figure(1); hold on;
set(gcf,'color','white')
plot(w(I),LCP(I)./max(LCP),'ok','MarkerFaceColor','k','MarkerSize',4)
plot(w(I),I_el(I)./max(I_el),'-r','LineStyle','-','LineWidth',1)

