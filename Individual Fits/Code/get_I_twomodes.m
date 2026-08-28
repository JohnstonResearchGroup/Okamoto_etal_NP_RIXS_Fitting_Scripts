Amp_ph = P(1); 
wph1 = P(2);
M1 = P(3);
wph2 = P(4);
M2 = P(5);
Amp_el = P(6); 

load gam;
load sig;
load z;
% 
yel = Amp_el*lorentzian(x,0,sig);
yph = zeros(size(x));
yback = zeros(size(x));
I = find(x>0);
yback(I) = P(7).*x(I)./((x(I)-P(8)).^2 + P(9)^2);

g1 = M1/wph1;
g2 = M2/wph2;

for n1 = 0:1:8
  for n2 = 0:1:8
    if(n1+n2 ~= 0)
        delta = lorentzian(x,n1*wph1+n2*wph2,sig);
    else
        delta = zeros(size(x));
    end
    cross_section = 0;
    
    for m1 = 0:1:50
      for m2 = 0:1:50
         D1 = B(max(n1,m1),min(n1,m1),g1)*B(m1,0,g1);
         D2 = B(max(n2,m2),min(n2,m2),g2)*B(m2,0,g2);
         cross_section = cross_section + D1*D2/(z-(m1-g1*g1)*wph1-(m2-g2*g2)*wph2);
      end
    end
    yph = yph + Amp_ph*real(cross_section*(conj(cross_section)))*delta;   
  end
end


y1 = yel + yph + yback;

Itotal = zeros(size(x));
Iel = zeros(size(x));
Iph = zeros(size(x));
Iback = zeros(size(x));

for nw = 1:1:length(x)
   gauss = gaussian(x,x(nw),0.0154488/2.355);
   Itotal(nw) = trapz(x,-y1.*gauss);
      Iel(nw) = trapz(x,-yel.*gauss);
      Iph(nw) = trapz(x,-yph.*gauss);      
    Iback(nw) = trapz(x,-yback.*gauss);
end

y = Itotal*trapz(y1)/trapz(Itotal);
if(trapz(x,Iel)~=0) yel = Iel*trapz(yel)/trapz(Iel); end
if(trapz(x,Iback)~=0) yback = Iback*trapz(yback)/trapz(Iback); end
if(trapz(x,Iph)~=0) yph = Iph*trapz(yph)/trapz(Iph); end

norm = max(y);
y = y/norm;
yph = yph/norm;
yback = yback/norm;
yel = yel/norm;
