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
I = find(w <= 0.02);
LCP = A(:,2)./max(A(I,2));    %normalize!
RCP = A(:,3)./max(A(I,3));

%Save the parameters for use in the fitting routines.
save Filename Filename
save gam gam;
save sig sig;
save z z;

%P = [6.4348    6.3240   23.5216    1.9143    1.5776    0.0845    0.0259  0.1788    0.1767    0.0532]
P0 = [4.9659 5.23 36.1594 1.4307 1.1891 0.0850 0.0371 0.1690 0.1725 0.0400];

% Use fmincon to refine the fit and obtain parameter error estimates
% Define bounds (example: enforce positivity where appropriate)
lb = [ 0,  0,   0,  0,  0,    0,    0, 0, 0, 0];   % lower bounds
ub = [10, 10, 100, 10, 10, 0.15, 0.15, 1, 1, 1]; % upper bounds

% Convert options for fmincon
opt_fmincon = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',4000, ...
    'MaxIterations',100,'OptimalityTolerance',1e-12,'StepTolerance',1e-12);

% Run fmincon starting from P (result of fminsearch)
[P_con, fval_con, exitflag, output, lambda, grad, hessian2] = ...
    fmincon(@chi2, P0, [], [], [], [], lb, ub, [], opt_fmincon);

% Compute numerical Hessian of chi2 at P_con
% Replace P with refined estimate
P = P_con;
compute_hessian

% Estimate parameter covariance and standard errors
% For least-squares-like chi2, approximate covariance = 2 * inv(H)
% (factor depends on definition of chi2; adjust if needed)
stderrP = 2*sqrt(diag(inv(hessian)));

% Replace P with refined estimate
P = P_con;

% Save refined parameters and errors
save parameters_fit_dataA.mat P stderrP hessian fval_con exitflag output lambda
chi2(P)

get_intensities;

%Plot the results
figure('Renderer', 'painters', 'Position', [10 10 600 700]); 
set(gcf,'color','white')
subplot('position',[0.125,0.55,0.83,0.43]); hold on; 
plot(w,LCP,'ok','MarkerFaceColor','k','MarkerSize',4)
plot(w,y_total_LCP,'-r','LineWidth',1)
area(w,y_el_LCP,'EdgeColor','none','FaceColor','g','ShowBaseLine','off','FaceAlpha',0.25)
area(w,y_ph_LCP,'EdgeColor','none','FaceColor','m','ShowBaseLine','off','FaceAlpha',0.25)
area(w,y_back_LCP,'EdgeColor','none','FaceColor','b','ShowBaseLine','off','FaceAlpha',0.25)
box on;
axis([-0.1,0.5,0,1.2])
set(gca,'XTick',(-0.1:0.1:0.5),'YTick',(0:0.2:1),'FontSize',20,'FontName','Arial','LineWidth',1)
ylabel('Intensity (arb. units)','FontSize',20,'FontName','Arial')
xticklabels({' ', ' ', ' ', ' ', ' ', ' ', ' '})
text(-0.085,1.1,'(a) LCP','FontSize',20,'FontName','Arial')
legend('Exp.','Fit','Elastic Peak','Phonons','Background')
legend('Location','northeast','FontName','Arial','TextColor','k','FontSize',20,'Box','off');

subplot('position',[0.125,0.1,0.83,0.43]); hold on; 
plot(w,RCP,'ok','MarkerFaceColor','k','MarkerSize',4)
plot(w,y_total_RCP,'-r','LineWidth',1)
box on;
axis([-0.1,0.5,0,1.2])
xlabel('Energy Loss (eV)','FontSize',20,'FontName','Arial')
set(gca,'XTick',(0:0.1:0.5),'YTick',(0:0.2:1),'FontSize',20,'FontName','Arial','LineWidth',1)
xlabel('Energy Loss (eV)','FontSize',20,'FontName','Arial')
ylabel('Intensity (arb. units)','FontSize',20,'FontName','Arial')
text(-0.085,1.1,'(b) RCP','FontSize',20,'FontName','Arial')
area(w,y_el_LCP,'EdgeColor','none','FaceColor','g','ShowBaseLine','off','FaceAlpha',0.25)
area(w,y_ph_LCP,'EdgeColor','none','FaceColor','m','ShowBaseLine','off','FaceAlpha',0.25)
area(w,y_back_LCP,'EdgeColor','none','FaceColor','b','ShowBaseLine','off','FaceAlpha',0.25)

%Save the results.
saveas(gcf,'fit_dataA.png','png')
saveas(gcf,'fit_dataA.pdf','pdf')

dlmwrite('fitted_curves_dataA_LCP.txt',[w,y_total_LCP,y_el_LCP,y_ph_LCP,y_back_LCP],' ')
dlmwrite('fitted_curves_dataA_RCP.txt',[w,y_total_RCP,y_el_RCP,y_ph_RCP,y_back_RCP],' ')
