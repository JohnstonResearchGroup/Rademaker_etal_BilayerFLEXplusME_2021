%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Run this script to plot Figures 2a, 2b, and 2c from the paper. 
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

figs = ['a' 'b' 'c'];
Ulist = [7 0 7];
tzlist = [2.36 4.21 2.36];
lamlist = [0.0 0.19 0.19];

mytemp = 0.01125; %single temperature value chosen for bandstrcuture plots

%conversion of T for t = 75meV -> [Kelvin]/conv = [t]
conv = 870.3393754624;
%Size of momentum grid
Nk = [32 32];
for ii = 1:2
    K{ii} = (0:1:2*Nk(ii)-1)*pi/Nk(ii);
end
[KX,KY] = meshgrid(K{2},K{1});      %row: y, column: x

ff = 20;
for iiu=1:numel(figs)
    ff = ff +1;
    lambda_m = lamlist(iiu);
    U = Ulist(iiu);
    tz = tzlist(iiu);
    
    fileDir = ['./tz_' num2str(tz) '/U_' num2str(U) '/lambda_' num2str(lambda_m,'%g') '/'];
    
    T = mytemp;
    
    %load the calculated self energies and energy dispersion
    load([fileDir 'selfE_Nk=' num2str(Nk(1)) '_' num2str(Nk(2))  '_T=' num2str(T) '_U=' num2str(U) '_freqFFT.mat'],...
        'Nk','numwi','mu','WN','WNU','S0','P0','ek','-mat')
    
    P0pp = P0(:,:,1) + P0(:,:,2);
    P0mm = P0(:,:,1) - P0(:,:,2);
    S0pp = S0(:,:,1) + S0(:,:,2);
    S0mm = S0(:,:,1) - S0(:,:,2);
    S0ppC = conj(S0(:,:,1)) + conj(S0(:,:,2));
    S0mmC = conj(S0(:,:,1)) - conj(S0(:,:,2));
    Zpp = ones(2*Nk(1),2*Nk(2)) - ((S0pp - S0ppC)/(2*1i*pi*T));
    Zmm = ones(2*Nk(1),2*Nk(2)) - ((S0mm - S0mmC)/(2*1i*pi*T));
    %interacting band structure for band 1
    eplus = real(((ek(:,:,1)+ek(:,:,2))-(mu*ones(size(S0(:,:,1))))+(S0(:,:,1)+S0(:,:,2))));
    atemp = eplus(1,1:1:Nk(1)+1);
    btemp = eplus(1:1:Nk(2)+1,Nk(1)+1);
    btemp = transpose(btemp);
    ctemp = diag(eplus(Nk(2)+1:-1:1,Nk(1)+1:-1:1));
    ctemp = transpose(ctemp);
    eplotplus = [atemp btemp ctemp];
    %non-interacting band structure for band 1
    neplus = real(((ek(:,:,1)+ek(:,:,2)))-(mu*ones(size(S0(:,:,1)))));
    atemp = neplus(1,1:1:Nk(1)+1);
    btemp = neplus(1:1:Nk(2)+1,Nk(1)+1);
    btemp = transpose(btemp);
    ctemp = diag(neplus(Nk(2)+1:-1:1,Nk(1)+1:-1:1));
    ctemp = transpose(ctemp);
    neplotplus = [atemp btemp ctemp];
    %interacting band structure for band 2
    eminus = real(((ek(:,:,1)-ek(:,:,2))-(mu*ones(size(S0(:,:,1))))+(S0(:,:,1)-S0(:,:,2))));
    atemp = eminus(1,1:1:Nk(1)+1);
    btemp = eminus(1:1:Nk(2)+1,Nk(1)+1);
    btemp = transpose(btemp);
    ctemp = diag(eminus(Nk(2)+1:-1:1,Nk(1)+1:-1:1));
    ctemp = transpose(ctemp);
    eplotminus = [atemp btemp ctemp];
    %non-interacting band structure for band 2
    neminus = real(((ek(:,:,1)-ek(:,:,2)))-(mu*ones(size(S0(:,:,1)))));
    atemp = neminus(1,1:1:Nk(1)+1);
    btemp = neminus(1:1:Nk(2)+1,Nk(1)+1);
    btemp = transpose(btemp);
    ctemp = diag(neminus(Nk(2)+1:-1:1,Nk(1)+1:-1:1));
    ctemp = transpose(ctemp);
    neplotminus = [atemp btemp ctemp];
    %superconducting gap as function of momentum for both bands
    sgapplus = real(P0pp)./Zpp;
    atemp = sgapplus(1,1:1:Nk(1)+1);
    btemp = sgapplus(1:1:Nk(2)+1,Nk(1)+1);
    btemp = transpose(btemp);
    ctemp = diag(sgapplus(Nk(2)+1:-1:1,Nk(1)+1:-1:1));
    ctemp = transpose(ctemp);
    sgapplusplot = 5*[atemp btemp ctemp];
    sgapminus = real(P0mm)./Zmm;
    atemp = sgapminus(1,1:1:Nk(1)+1);
    btemp = sgapminus(1:1:Nk(2)+1,Nk(1)+1);
    btemp = transpose(btemp);
    ctemp = diag(sgapminus(Nk(2)+1:-1:1,Nk(1)+1:-1:1));
    ctemp = transpose(ctemp);
    sgapminusplot = 5*[atemp btemp ctemp];
    
    if min(min(sgapminusplot)) < 0
        sgapminusplot = sgapminusplot * -1;
        sgapplusplot = sgapplusplot * -1;
        sgapminus = sgapminus * -1;
        sgapplus = sgapplus * -1;
    end
    %BANDSTRUCTURE PLOTTING BELOW
    
    egap = 0 - max(eplotplus);
    %fprintf('min of electron band=%g meV,max of hole band=%g meV,egap = %g meV for U=%gt lambda=%g T=%gK\n',min(75* eplotminus),max(75* eplotplus),75* egap,U,lammy,T*conv)
    figure(ff);hold on;
    set(gcf,'color','white','Position',get(0,'Screensize'));
    plot(75*eplotplus,'b','LineWidth',2.0);
    plot(75*eplotminus,'r','LineWidth',2.0);
    plot(75*neplotplus,'--b','LineWidth',2.0);
    plot(75*neplotminus,'--r','LineWidth',2.0);
    plot(75*sgapplusplot,'ob','LineWidth',2.0);
    plot(75*sgapminusplot,'or','LineWidth',2.0);
    if figs(iiu) == 'a'
        text(5,600,['\textbf{a}'],'FontSize',35);
    elseif figs(iiu) == 'b'
        text(5,600,['\textbf{b}'],'FontSize',35);
    elseif figs(iiu) == 'c'
        text(5,600,['\textbf{c}'],'FontSize',35);
    end
    text(5,450,['$\lambda=$ ' num2str(lambda_m)],'FontSize',35);
    text(5,350,['$U=$ ' num2str(U) 't'],'FontSize',35);
    text(5,250,['$t_\perp=$ ' num2str(tz) 't'],'FontSize',35);
    %text(5,600,['T = ' num2str(T*conv) 'K'],'FontSize',35);
    %title(['U=' num2str(U) 't'],'FontSize',10,'FontName','Latin Modern Roman','Units', 'normalized','Position', [-0.1, -0.15, 0]);
    ylabel('$\epsilon(\textbf{k}), 5\times\Delta(\textbf{k},$i$\pi$T$)$ [meV]','interpreter','latex','FontSize',35,'FontName','Latin Modern Roman');
    set(gca,'XTick',[1,33,66,99],'YTick',[-750:150:750],'FontSize',35,'FontName','cmr12');
    axis([1,99,-750,750]);
    xticklabels({'$\Gamma$','X','M','$\Gamma$'});
    %xticklabels({' ',' ',' ',' '});
    yticklabels({' ','-600',' ','-300',' ','0',' ','300',' ','600',' '})
    box on;
    yy = 75*[-20:1:20];
    xx = ones(size(yy));
    plot(33*xx,yy,'k');
    plot(66*xx,yy,'k');
    plot(99*xx,yy,'k');
    yline(0,'--b');
    %yline(min(75*eplotminus),'c');
    yline(max(75*eplotplus),'g','LineWidth',1.0);
    plotarray =[plot(33*xx,yy,'k') plot(66*xx,yy,'k') plot(99*xx,yy,'k') yline(0,'--b')];
    set(get(get(plotarray(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(plotarray(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(plotarray(3),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(plotarray(4),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    legend('E$_+(\textbf{k})$','E$_-(\textbf{k})$','$\epsilon_{+}(\textbf{k})$','$\epsilon_{-}(\textbf{k})$','$\Delta_{+}(\textbf{k})$','$\Delta_{-}(\textbf{k})$','FontSize',28,'FontName','Latin Modern Roman','Location','northeast');
    lgd = legend;
    lgd.NumColumns = 3;
    legend('boxoff')
    hold off;
end