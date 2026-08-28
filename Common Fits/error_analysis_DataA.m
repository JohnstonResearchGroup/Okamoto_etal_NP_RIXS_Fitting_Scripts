clc;
close all;
clear;

%Set the parameters related to the RIXS cross-section
gam = 0.200;                    %Inverse corehole lifetime in eV.
sig = 0.009;                    %Lifetime parameter, elastic line
sig2= 0.016;                    %Lifetime parameter, phonons
Delta = 0;                      %Detuning energy.
z = Delta + 1i*gam;             %

Filename = '../Data/A69Kdata.txt';
A = load(Filename);
w = A(:,1);
I = find(w >= -0.3 & w <= 0.5);
LCP = A(:,2)./max(A(I,2));    %normalize!
RCP = A(:,3)./max(A(I,3));

%compute the noise
I = find(w >= -0.35 & w <= -0.2);
sigma_noise = std(LCP(I),RCP(I));

%Save the parameters for use in the fitting routines.
save Filename Filename
save gam gam;
save sig sig;
save z z;
save sigma_noise sigma_noise;

load parameters_fit_dataA.mat;

lb = [ 0,  0,   0,  0,  0,    0,    0, 0, 0, 0];   % lower bounds
ub = [10, 10, 100, 10, 10, 0.15, 0.15, 1, 1, 1]; % upper bounds

% Convert options for fmincon
opt_fmincon = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',5000, ...
    'MaxIterations',1000,'OptimalityTolerance',1e-6,'StepTolerance',1e-6);

N_boot = 30;
P_boot = [];
fcon_boot = [];
for n = 1:N_boot
    % Run fmincon starting from P (result of fminsearch)
    [P_con, fval_con, exitflag, output, lambda, grad, hessian2] = ...
        fmincon(@chi2_wnoise, P+0.05*P.*rand(size(P)), [], [], [], [], lb, ub, [], opt_fmincon);
    P_boot = [P_boot;P_con];
    fcon_boot = [fcon_boot;fval_con];
end
save error_analysis_dataA.mat P_boot fcon_boot