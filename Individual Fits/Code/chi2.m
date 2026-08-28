function out = chi2(P)
load Filename;
load pol;
A = load(Filename);
x = A(:,1);
get_I_twomodes;

boundary = 0; 
if(P(1) < 0) boundary = boundary + 10000; end
if(P(6) < 0) boundary = boundary + 10000; end
if(P(7) < 0) boundary = boundary + 10000; end
if(P(8) < 0) boundary = boundary + 10000; end
if(P(9) < 0) boundary = boundary + 10000; end
out = sum((A(:,pol)-y).^2) + boundary;
