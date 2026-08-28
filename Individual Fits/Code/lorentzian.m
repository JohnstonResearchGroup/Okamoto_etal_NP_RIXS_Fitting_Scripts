function out = lorentzian(x,x0,sig)
out = (sig/(2*pi))./((x-x0).^2 + (sig/2)^2);