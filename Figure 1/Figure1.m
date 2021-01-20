%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Run this script to plot Figure 1 from the paper.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear Matlab
clc
clear
close all
%
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

fileDir = ['./'];
data = load([fileDir 'Tc_vs_lambda.dat']);

%conversion of T for t = 75meV -> [Kelvin]/conv = [t]
conv = 870.3393754624;

tcvslam1 = data(:,2); %Tc values for U = 0
tcvslam2 = data(:,3); %Tc values for U =7
lamlist = data(:,1);  %strengths of electron-phonon couplings (the given
                      % Fermi surface average mass enhancements)
                      
figure(1);hold on;
set(gcf,'color','white','Position',get(0,'Screensize'));
f1 = plot(lamlist,tcvslam1*conv,'Marker','s','MarkerSize',20,'Color','r','LineWidth',2.5,'DisplayName','U = 0');
set(f1, 'markerfacecolor', get(f1, 'color'));
f2 = plot(lamlist,tcvslam2*conv,'Marker','o','MarkerSize',20,'Color','b','LineWidth',2.5,'DisplayName','U = 7t');
set(f2, 'markerfacecolor', get(f2, 'color'));
xlabel('$\lambda_{m}$','FontSize',35);
ylabel('$T_{c}$ [K]','Fontsize',35);
set(gca,'XTick',[0:0.01:0.25],'YTick',[0:5:120],'FontSize',35,'FontName','Latin Modern Roman');
yticklabels({'0',' ',' ',' ','20',' ',' ',' ','40',' ',' ',' ','60',' ',' ',' ','80',' ',' ',' ','100',' ',' ',' ','120'});
xticklabels({'0',' ',' ',' ',' ','0.05',' ',' ',' ',' ','0.1',' ',' ',' ',' ','0.15',' ',' ',' ',' ','0.20',' ',' ',' ',' ','0.25'});
legend('FontSize',35,'FontName','Latin Modern Roman','Location','northwest')
lgd = legend;
lgd.NumColumns = 1;
box on;
axis([0.0 0.25 0.0 120]);
legend('boxoff')
hold off;