h = 2*1e-3;
n = numel(P_con);
hessian = zeros(n);
f0 = chi2(P_con);
for i=1:n
    ei = zeros(1,n); ei(i)=1;
    fp = chi2(P + h*ei);
    fm = chi2(P - h*ei);

    hessian(i,i) = (fp - 2*f0 + fm)/(h*h);

    for j=i+1:n
        ej = zeros(1,n); ej(j)=1;
        fpp = chi2(P + h*ei + h*ej);
        fmp = chi2(P - h*ei + h*ej);
        fpm = chi2(P + h*ei - h*ej);
        fmm = chi2(P - h*ei - h*ej);
        hessian(i,j) = (fpp - fmp - fpm + fmm)/(4*h*h);
        hessian(j,i) = hessian(i,j);
        disp([i,j]);
    end
end