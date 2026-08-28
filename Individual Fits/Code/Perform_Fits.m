%=========================================================================================
% This MATLAB script will attempt fits of the quasi-elastic RIXS intensity of the data. 
% Script by S. Johnston, 28 March 2025.
%=========================================================================================

% Start with setting all environment variables
clear;                                      %Clear Memory
clc;                                        %Clear the command line
close all;                                  %Close all prior figures

%Set the default interpreter to Latex for plotting purposes
set(0,'defaulttextInterpreter','latex')     

pol = 2; %LCP
pol = 3; %RCP


%Set the input file. The file input is assumed to be two columns. The first column 
%is the energy loss. The second column is the RIXS intensity.
Filename = '../Data/B69Kdata_0803.txt';
A = load(Filename);
A(:,1) = A(:,1);
A(:,pol) = A(:,pol)./max(A(:,pol));

%Set the options for the fitting routine. 
opt = optimset('Display','iter','MaxFunEvals',400,'MaxIter',5000,'TolX',1e-12,'TolFun',1e-12);
% 
% %Set the parameters related to the RIXS cross-section
gam = 0.250;                    %Inverse corehole lifetime in eV.
sig = 0.016;                   %Resolution parameter
Delta = 0;                      %Detuning energy.
z = Delta + 1i*gam;             %

%Determine the background
I = find(A(:,1)<-0.4 & A(:,1) > -0.7);
Background = mean(A(I,pol));  

%Save the parameters for use in the fitting routines.
save Filename Filename
save gam gam;
save sig sig;
save z z;
save Background Background;
save pol pol;

% % Results for LCP (Column 2)
% P0 = [0.7569    0.0883    0.1516    0.0717    0.0894    1.7870    2.8936    0.0003    0.2966];

% Results for RCP (Column 3)
P0 = [0.8767    0.0887    0.1313    0.0727    0.0840    1.7376    2.5175    0.0000    0.1933];
P = P0;
[P,fval] = fminsearch(@chi2,P0,opt);
chi2(P)
x = A(:,1);
get_I_twomodes;

%Plot the results
figure('Renderer', 'painters', 'Position', [10 10 900 600]); hold on

plot(A(:,1),A(:,2)-Background*ones(size(A(:,1))),'ok','MarkerFaceColor','k','MarkerSize',3)
plot(x,y,'-r','linewidth',1)
area(x,yel,'EdgeColor','none','FaceColor','g','ShowBaseLine','off','FaceAlpha',0.25)
area(x,yph,'EdgeColor','none','FaceColor','m','ShowBaseLine','off','FaceAlpha',0.25)
area(x,yback,'EdgeColor','none','FaceColor','b','ShowBaseLine','off','FaceAlpha',0.25)

% Make figure and axes background white

set(gcf,'Color','white');
set(gca,'Color','white');
box on;

% Set axis font to black Times New Roman
set(gca, 'FontName', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 12, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
axis([-0.1,0.5,0,1.1])

% Also set labels and title (if present) to Times New Roman and black
hx = get(gca, {'Title','XLabel','YLabel'});
for k = 1:numel(hx)
    if ~isempty(hx{k})
        set(hx{k}, 'FontName', 'Times New Roman', 'Color', 'k','Interpreter','latex','FontSize',20);
    end
end

set(gca,'XTick',(-0.1:0.1:0.5),'YTick',(0:0.2:1),'FontSize',20)
xlabel('Energy Loss (eV)','FontSize',25)
ylabel('RIXS Intensity (arb. units)','FontSize',25)
legend('Exp.','Fit','Elastic Peak','Phonons','Background')
legend('Location','northeast','FontName','Times New Roman','TextColor','k','FontSize',20,'Box','off');

if(pol == 3)
    saveas(gcf,'./fits_RCP.png','png')
elseif(pol==2)
    saveas(gcf,'./fits_LCP.png','png')
end