clear;
clc;
close all;

set(groot, 'DefaultLineLineWidth', 1)

A = load('../Data/B69Kdata_0803.txt');

figure('Renderer', 'painters', 'Position', [10 10 900 600]);
hold on;
plot(A(:,1),A(:,2),'-r', 'DisplayName', 'LCP')
plot(A(:,1),A(:,3),'-b', 'DisplayName', 'RCP')

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
legend('Location','northeast','FontName','Times New Roman','TextColor','k','FontSize',20,'Box','off');
saveas(gcf,'../Figures/exp_data.png','png')
saveas(gcf,'../Figures/exp_data.eps','eps')