function out = B(m,n,alpha)
% pre = exp(-alpha*alpha/2)*sqrt(factorial(n)*factorial(m));
% tmp = 0;
% for l = 0:1:n
%     num = ((-1)^(m+l))*(alpha^(2*l+m-n));
%     num = num/factorial(n-l);
%     num = num/factorial(l);
%     num = num/factorial(m-n+l);
%     tmp = tmp + num;
% end
% B1 = pre*tmp;
% B2 = ((-1)^m)*exp(-alpha*alpha/2)*(alpha^(m-n))*sqrt(factorial(n)/factorial(m))* ...
%      polyval(LaguerreGen(n, m-n), alpha*alpha);
out = ((-1)^m)*exp(-alpha*alpha/2)*(alpha^(m-n))*sqrt(factorial(n)/factorial(m))* ...
      polyval(LaguerreGen(n, m-n), alpha*alpha);