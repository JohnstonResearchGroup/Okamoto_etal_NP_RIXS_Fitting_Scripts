
Amp_el = P(1); 
Amp_back = P(2);
beta = P(3);
Amp_ph_LCP = P(4);
Amp_ph_RCP = P(5);
wph1 = P(6);
wph2 = P(7);
M1_LCP = P(8);
M1_RCP = P(9);
M2 = P(10);

g1_LCP = M1_LCP/wph1;
g1_RCP = M1_RCP/wph1;
g2 = M2/wph2;

yel = Amp_el*lorentzian(w,0,sig);
yback = Amp_back*(1+tanh(w*beta));


yph_LCP = zeros(size(w));
yph_RCP = zeros(size(w));

for n1 = 0:1:8
  for n2 = 0:1:8
    if(n1+n2 ~= 0)
        delta = lorentzian(w,n1*wph1+n2*wph2,sig2);
     else
        delta = zeros(size(w));
    end
    cross_section_LCP = 0;
    cross_section_RCP = 0;
    
    for m1 = 0:1:20
      for m2 = 0:1:20
         D1_LCP = B(max(n1,m1),min(n1,m1),g1_LCP)*B(m1,0,g1_LCP);
         D1_RCP = B(max(n1,m1),min(n1,m1),g1_RCP)*B(m1,0,g1_RCP);
         D2 = B(max(n2,m2),min(n2,m2),g2)*B(m2,0,g2);
         cross_section_LCP = cross_section_LCP ...
                           + D1_LCP*D2/(z-(m1-g1_LCP*g1_LCP)*wph1-(m2-g2*g2)*wph2);
         cross_section_RCP = cross_section_RCP ...
                           + D1_RCP*D2/(z-(m1-g1_RCP*g1_RCP)*wph1-(m2-g2*g2)*wph2);
      end
    end
    yph_LCP = yph_LCP + Amp_ph_LCP*real(cross_section_LCP*(conj(cross_section_LCP)))*delta;  
    yph_RCP = yph_RCP + Amp_ph_RCP*real(cross_section_RCP*(conj(cross_section_RCP)))*delta;  
  end
end

y_total_LCP = yel + yback + yph_LCP;
y_total_RCP = yel + yback + yph_RCP;

%Normalize the spectra
norm2_LCP = max(y_total_LCP);
y_total_LCP = y_total_LCP/norm2_LCP;

norm2_RCP = max(y_total_RCP);
y_total_RCP = y_total_RCP./max(norm2_RCP);

% resolution broaden the spectra
I_LCP = zeros(size(w));
I_RCP = zeros(size(w)); 
I_back = zeros(size(w));
I_el = zeros(size(w));
I_ph_LCP = zeros(size(w));
I_ph_RCP = zeros(size(w));

for nw = 1:1:length(w)
   gauss = gaussian(w,w(nw),0.024/2.355);
   I_LCP(nw) = trapz(w,-y_total_LCP.*gauss);
   I_RCP(nw) = trapz(w,-y_total_RCP.*gauss);
   I_ph_LCP(nw) = trapz(w,-yph_LCP.*gauss);
   I_ph_RCP(nw) = trapz(w,-yph_RCP.*gauss);   
   I_back(nw) = trapz(w,-yback.*gauss);
   I_el(nw) = trapz(w,-yel.*gauss);
end

norm_LCP = max(I_LCP);
y_total_LCP = I_LCP/norm_LCP;
y_back_LCP = I_back/(norm_LCP*norm2_LCP);
y_el_LCP = I_el/(norm_LCP*norm2_LCP);
y_ph_LCP = I_ph_LCP/(norm_LCP*norm2_LCP);

norm_RCP = max(I_RCP);
y_total_RCP = I_RCP/norm_RCP;
y_back_RCP = I_back/(norm_RCP*norm2_RCP);
y_el_RCP = I_el/(norm_RCP*norm2_RCP);
y_ph_RCP = I_ph_RCP/(norm_RCP*norm2_RCP);