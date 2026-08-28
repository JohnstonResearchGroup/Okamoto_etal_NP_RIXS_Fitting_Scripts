function out = gaussian(x,x0,sig)
out = exp(-((x-x0)/(sqrt(2)*sig)).^2)/(sig*sqrt(2*pi));