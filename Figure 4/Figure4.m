%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Run this script to plot Figures 4 from the paper. 
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

q0 = 0.1;
lamlist =[0.0 0.05 0.10 0.14 0.19 0.23];% list of lambda_m values, the
%control parameter for strength of
%e-ph coupling.
Ulist = [7]; %U value in units of t
tzlist = [2.36]; %corresponding t_\perp value in units of t

%conversion of T for t = 75meV -> [Kelvin]/conv = [t]
conv = 870.3393754624;
%Size of momentum grid
Nk = [32 32];
for ii = 1:2
    K{ii} = (0:1:2*Nk(ii)-1)*pi/Nk(ii);
end
[KX,KY] = meshgrid(K{2},K{1});      %row: y, column: x

gg = 4;
ss = 40;

maxgapp = [];
maxgapm = [];

myfit = fittype('a*tanh(b*sqrt(max(1-x/c,0)))','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b','c'});

%begin lambda loop
for iiu=1:numel(lamlist)
    ss = ss+1;
    
    lambda_m = lamlist(iiu);
    U = Ulist(1);
    tz = tzlist(1);
    dispnm = ['$\lambda_m=$ ' num2str(lambda_m)];
    dispnm1 = ['$\mid\Delta_{+}\mid$'];
    dispnm2 = ['$\mid\Delta_{-}\mid$'];
    %Tc values for U=7 withlambda_m = 0, 0.05, 0.10, 0.14, 0.19, and 0.23,
    %respectively. In units of [t]
    tclist = [0.0538083 0.077715 0.0950005 0.105015 0.115215 0.122526];
    %temperature for last non-zero data point of the gap. Used in Tc fit.
    neartclist = [0.05375 0.0775 0.095 0.105 0.115 0.1225];
    
    fileDir = ['./tz_' num2str(tz) '/U_' num2str(U) '/lambda_' num2str(lambda_m) '/'];
    
    %load temperature list
    gaps = load([fileDir 'gaps_Nk=' num2str(Nk(1)) '_' num2str(Nk(2)) '_U=' num2str(U) '_lam=' num2str(lambda_m) '_q0=' num2str(q0) '.dat']);
    Tlist = gaps(:,1);
    
    %begin temperature loop
    sgapvalplus =[];
    sgapvalminus =[];
    
    for nt = 1:numel(Tlist)
        T = Tlist(nt);
        %load data
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
        eplus = real(((ek(:,:,1)+ek(:,:,2))-(mu*ones(size(S0(:,:,1))))+(S0(:,:,1)+S0(:,:,2))));
        eminus = real(((ek(:,:,1)-ek(:,:,2))-(mu*ones(size(S0(:,:,1))))+(S0(:,:,1)-S0(:,:,2))));
        sgapplus = real(P0pp)./Zpp;
        sgapminus = real(P0mm)./Zmm;
        
        if min(min(sgapminus)) < 0
            sgapminus = sgapminus * -1;
            sgapplus = sgapplus * -1;
        end
        %dE = (2*pi)/(Nk(1)+1);
        %mask_FS = abs(eminus(:,:)) < dE;
        A = nnormgauss(eminus,1);
        sgapvalminus(nt) = max(max((75*A.*sgapminus)));
        if max(max(abs(sgapplus - sgapminus))) < 0.0001
            sgapvalplus(nt) =  max(max(abs(75*A.*sgapplus)));
        else
            sgapvalplus(nt) =  min(min((75*sgapplus)));
        end
        %sgapvalplus(nt) = 75*(max(max(abs(real(P0pp./Zpp)))));
        %sgapvalminus(nt) = 75*(max(max(abs(real(P0mm./Zmm)))));
        %%%%
        %Zee = 1 - imag(S0(:,:,1))/WN(numwi+1);
        %zval(nt) = max(max(Zee));
    end%end temperature loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %                       GAP VS TEMPERATURE PLOT
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(gg);hold on;
    set(gcf,'color','white','Position',get(0,'Screensize'));
    if iiu == 1
        temp = Tlist;
        temp(numel(Tlist)+1) = tclist(iiu);
        sgapp = sgapvalplus;
        sgapm = sgapvalminus;
        sgapp(numel(sgapvalplus)+1) = 0;
        sgapm(numel(sgapvalminus)+1) = 0;
        p2 = plot(temp, sgapp,'LineStyle','-','Marker','none','MarkerSize',3,'Color','#9400D3','LineWidth',3,'DisplayName',dispnm);
        p1 = plot(temp, sgapm,'LineStyle','-','Marker','none','MarkerSize',3,'Color','#9400D3','LineWidth',3,'DisplayName',dispnm);
        set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        set(p1, 'markerfacecolor', get(p1, 'color'));
        set(p2, 'markerfacecolor', get(p2, 'color'));
        maxgapp(iiu) = sgapvalplus(1);
        maxgapm(iiu) = sgapvalminus(1);
        fprintf('max value of Delta minus =%g [meV]\n',maxgapm(iiu))
        
        figure(ss);hold on;
        set(gcf,'color','white','Position',get(0,'Screensize'));
        Tcfit{iiu} = fit(Tlist,transpose(sgapvalminus),myfit,'StartPoint',[sgapvalminus(1),10,neartclist(iiu)]);
        temperature=linspace(0,0.5,1000);
        temp_gap{iiu} = (Tcfit{iiu}.a)*tanh((Tcfit{iiu}.b)*sqrt( 1-temperature./(Tcfit{iiu}.c)) );
        fprintf('Tc fit value for U = %g and lambda=%g is %g --> %g K\n',U,lambda_m,Tcfit{iiu}.c,Tcfit{iiu}.c*conv)
        fprintf('Gap/Tc =%g\n',((maxgapm(iiu)/75)/tclist(iiu)))
        plot(temperature,real(temp_gap{iiu}),'--b','LineWidth',2,'DisplayName',['$T_\textrm{c} \approx$' num2str((Tcfit{iiu}.c)*conv) 'K'])
        plot(Tlist, abs(sgapvalplus), 'Marker','o','Color','#9400D3','LineWidth',1.5,'DisplayName',dispnm1);
        plot(Tlist, abs(sgapvalminus), 'Marker','x','Color','#9400D3','LineWidth',1.5,'DisplayName',dispnm2);
        text(10/conv,2,['$\lambda_m = $' num2str(lambda_m)],'FontSize',35);
        xlabel('T[K]');
        ylabel('$\mid\Delta_{\pm}(\textbf{k}^\textrm{max}_{\pm},$i$\pi$T$)\mid$ [meV]','Fontsize',35);
        axis([0.0 120/conv 0 35]);
        set(gca,'XTick',[0,10/conv,20/conv,30/conv,40/conv,50/conv,60/conv,70/conv,80/conv,90/conv,100/conv,110/conv,120/conv,130/conv,140/conv,150/conv],'YTick',[0:5:50],'FontSize',35,'FontName','Times');
        xticklabels({'0','10','20','30','40','50','60','70','80','90','100','110','120','130','140','150'});
        legend('FontSize',28,'FontName','Latin Modern Roman','Location','northeast')
        lgd = legend;
        lgd.NumColumns = 1;
        legend('boxoff')
        box on;
        
        figure(gg);
    elseif iiu == 2
        temp = Tlist;
        temp(numel(Tlist)+1) = tclist(iiu);
        sgapp = sgapvalplus;
        sgapm = sgapvalminus;
        sgapp(numel(sgapvalplus)+1) = 0;
        sgapm(numel(sgapvalminus)+1) = 0;
        p2 = plot(temp, sgapp,'LineStyle','-','Marker','none','MarkerSize',3,'Color','#007BA7','LineWidth',3,'DisplayName',dispnm);
        p1 = plot(temp, sgapm,'LineStyle','-','Marker','none','MarkerSize',3,'Color','#007BA7','LineWidth',3,'DisplayName',dispnm);
        set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        set(p1, 'markerfacecolor', get(p1, 'color'));
        set(p2, 'markerfacecolor', get(p2, 'color'));
        maxgapp(iiu) = sgapvalplus(1);
        maxgapm(iiu) = sgapvalminus(1);
        fprintf('max value of Delta minus =%g [meV]\n',maxgapm(iiu))
        
        figure(ss);hold on;
        set(gcf,'color','white','Position',get(0,'Screensize'));
        Tcfit{iiu} = fit(Tlist,transpose(sgapvalminus),myfit,'StartPoint',[sgapvalminus(1),10,neartclist(iiu)]);
        temperature=linspace(0,0.5,1000);
        temp_gap{iiu} = (Tcfit{iiu}.a)*tanh((Tcfit{iiu}.b)*sqrt( 1-temperature./(Tcfit{iiu}.c)) );
        fprintf('Tc fit value for U = %g and lambda=%g is %g --> %g K\n',U,lambda_m,Tcfit{iiu}.c,Tcfit{iiu}.c*conv)
        fprintf('Gap/Tc =%g\n',((maxgapm(iiu)/75)/tclist(iiu)))
        plot(temperature,real(temp_gap{iiu}),'--b','LineWidth',2,'DisplayName',['$T_\textrm{c} \approx$' num2str((Tcfit{iiu}.c)*conv) 'K'])
        plot(Tlist, abs(sgapvalplus), 'Marker','o','Color','#FF8C00','LineWidth',2,'DisplayName',dispnm1);
        plot(Tlist, abs(sgapvalminus), 'Marker','x','Color','#FF8C00','LineWidth',2,'DisplayName',dispnm2);
        text(10/conv,2,['$\lambda_m = $' num2str(lambda_m)],'FontSize',35);
        xlabel('T[K]');
        ylabel('$\mid\Delta_{\pm}(\textbf{k}^\textrm{max}_{\pm},$i$\pi$T$)\mid$ [meV]','Fontsize',35);
        axis([0.0 120/conv 0 35]);
        set(gca,'XTick',[0,10/conv,20/conv,30/conv,40/conv,50/conv,60/conv,70/conv,80/conv,90/conv,100/conv,110/conv,120/conv,130/conv,140/conv,150/conv],'YTick',[0:5:50],'FontSize',35,'FontName','Times');
        xticklabels({'0','10','20','30','40','50','60','70','80','90','100','110','120','130','140','150'});
        legend('FontSize',28,'FontName','Latin Modern Roman','Location','northeast')
        lgd = legend;
        lgd.NumColumns = 1;
        legend('boxoff')
        box on;
        
        figure(gg);
    elseif iiu == 3
        temp = Tlist;
        temp(numel(Tlist)+1) = tclist(iiu);
        sgapp = sgapvalplus;
        sgapm = sgapvalminus;
        sgapp(numel(sgapvalplus)+1) = 0;
        sgapm(numel(sgapvalminus)+1) = 0;
        p2 = plot(temp, sgapp,'LineStyle','-','Marker','none','MarkerSize',2,'Color','#3B7A57','LineWidth',3,'DisplayName',dispnm);
        p1 = plot(temp, sgapm,'LineStyle','-','Marker','none','MarkerSize',2,'Color','#3B7A57','LineWidth',3,'DisplayName',dispnm);
        set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        set(p1, 'markerfacecolor', get(p1, 'color'));
        set(p2, 'markerfacecolor', get(p2, 'color'));
        maxgapp(iiu) = sgapvalplus(1);
        maxgapm(iiu) = sgapvalminus(1);
        fprintf('max value of Delta minus =%g [meV]\n',maxgapm(iiu))
        
        figure(ss);hold on;
        set(gcf,'color','white','Position',get(0,'Screensize'));
        Tcfit{iiu} = fit(Tlist,transpose(sgapvalminus),myfit,'StartPoint',[sgapvalminus(1),10,neartclist(iiu)]);
        temperature=linspace(0,0.5,1000);
        temp_gap{iiu} = (Tcfit{iiu}.a)*tanh((Tcfit{iiu}.b)*sqrt( 1-temperature./(Tcfit{iiu}.c)) );
        fprintf('Tc fit value for U = %g and lambda=%g is %g --> %g K\n',U,lambda_m,Tcfit{iiu}.c,Tcfit{iiu}.c*conv)
        fprintf('Gap/Tc =%g\n',((maxgapm(iiu)/75)/tclist(iiu)))
        plot(temperature,real(temp_gap{iiu}),'--b','LineWidth',2,'DisplayName',['$T_\textrm{c} \approx$' num2str((Tcfit{iiu}.c)*conv) 'K'])
        plot(Tlist, abs(sgapvalplus), 'Marker','o','Color','#DC143C','LineWidth',2,'DisplayName',dispnm1);
        plot(Tlist, abs(sgapvalminus), 'Marker','x','Color','#DC143C','LineWidth',2,'DisplayName',dispnm2);
        text(10/conv,2,['$\lambda_m = $' num2str(lambda_m)],'FontSize',35);
        xlabel('T[K]');
        ylabel('$\mid\Delta_{\pm}(\textbf{k}^\textrm{max}_{\pm},$i$\pi$T$)\mid$ [meV]','Fontsize',35);
        axis([0.0 120/conv 0 35]);
        set(gca,'XTick',[0,10/conv,20/conv,30/conv,40/conv,50/conv,60/conv,70/conv,80/conv,90/conv,100/conv,110/conv,120/conv,130/conv,140/conv,150,conv],'YTick',[0:5:50],'FontSize',35,'FontName','Times');
        xticklabels({'0','10','20','30','40','50','60','70','80','90','100','110','120','130','140','150'});
        legend('FontSize',28,'FontName','Latin Modern Roman','Location','northeast')
        lgd = legend;
        lgd.NumColumns = 1;
        legend('boxoff')
        box on;
        
        
        figure(gg);
    elseif iiu == 4
        temp = Tlist;
        temp(numel(Tlist)+1) = tclist(iiu);
        sgapp = sgapvalplus;
        sgapm = sgapvalminus;
        sgapp(numel(sgapvalplus)+1) = 0;
        sgapm(numel(sgapvalminus)+1) = 0;
        p2 = plot(temp, sgapp,'LineStyle','-','Marker','none','MarkerSize',2,'Color','#B0BF1A','LineWidth',3,'DisplayName',dispnm);
        p1 = plot(temp, sgapm,'LineStyle','-','Marker','none','MarkerSize',2,'Color','#B0BF1A','LineWidth',3,'DisplayName',dispnm);
        set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        set(p1, 'markerfacecolor', get(p1, 'color'));
        set(p2, 'markerfacecolor', get(p2, 'color'));
        maxgapp(iiu) = sgapvalplus(1);
        maxgapm(iiu) = sgapvalminus(1);
        fprintf('max value of Delta minus =%g [meV]\n',maxgapm(iiu))
        
        figure(ss);hold on;
        set(gcf,'color','white','Position',get(0,'Screensize'));
        Tcfit{iiu} = fit(Tlist,transpose(sgapvalminus),myfit,'StartPoint',[sgapvalminus(1),10,neartclist(iiu)]);
        temperature=linspace(0,0.5,1000);
        temp_gap{iiu} = (Tcfit{iiu}.a)*tanh((Tcfit{iiu}.b)*sqrt( 1-temperature./(Tcfit{iiu}.c)) );
        fprintf('Tc fit value for U = %g and lambda=%g is %g --> %g K\n',U,lambda_m,Tcfit{iiu}.c,Tcfit{iiu}.c*conv)
        fprintf('Gap/Tc =%g\n',((maxgapm(iiu)/75)/tclist(iiu)))
        plot(temperature,real(temp_gap{iiu}),'--b','LineWidth',2,'DisplayName',['$T_\textrm{c} \approx$' num2str((Tcfit{iiu}.c)*conv) 'K'])
        plot(Tlist, abs(sgapvalplus), 'Marker','o','Color','#03C03C','LineWidth',2,'DisplayName',dispnm1);
        plot(Tlist,abs(sgapvalminus), 'Marker','x','Color','#03C03C','LineWidth',2,'DisplayName',dispnm2);
        text(10/conv,2,['$\lambda_m = $' num2str(lambda_m)],'FontSize',35);
        xlabel('T[K]','Fontsize',35);
        ylabel('$\mid\Delta_{\pm}(\textbf{k}^\textrm{max}_{\pm},$i$\pi$T$)\mid$ [meV]','Fontsize',35);
        axis([0.0 120/conv 0 35]);
        set(gca,'XTick',[0,10/conv,20/conv,30/conv,40/conv,50/conv,60/conv,70/conv,80/conv,90/conv,100/conv,110/conv,120/conv,130/conv,140/conv,150/conv],'YTick',[0:5:50],'FontSize',35,'FontName','Times');
        xticklabels({'0','10','20','30','40','50','60','70','80','90','100','110','120','130','140','150'});
        legend('FontSize',28,'FontName','Latin Modern Roman','Location','northeast')
        lgd = legend;
        lgd.NumColumns = 1;
        legend('boxoff')
        box on;
        
        
        figure(gg);
    elseif iiu == 5
        temp = Tlist;
        temp(numel(Tlist)+1) = tclist(iiu);
        sgapp = sgapvalplus;
        sgapm = sgapvalminus;
        sgapp(numel(sgapvalplus)+1) = 0;
        sgapm(numel(sgapvalminus)+1) = 0;
        p2 = plot(temp, sgapp,'LineStyle','-','Marker','none','MarkerSize',2,'Color','#FF8C00','LineWidth',3,'DisplayName',dispnm);
        p1 = plot(temp, sgapm,'LineStyle','-','Marker','none','MarkerSize',2,'Color','#FF8C00','LineWidth',3,'DisplayName',dispnm);
        set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        set(p1, 'markerfacecolor', get(p1, 'color'));
        set(p2, 'markerfacecolor', get(p2, 'color'));
        maxgapp(iiu) = sgapvalplus(1);
        maxgapm(iiu) = sgapvalminus(1);
        fprintf('max value of Delta minus =%g [meV]\n',maxgapm(iiu))
        
        figure(ss);hold on;
        set(gcf,'color','white','Position',get(0,'Screensize'));
        Tcfit{iiu} = fit(Tlist,transpose(sgapvalminus),myfit,'StartPoint',[sgapvalminus(1),10,neartclist(iiu)]);
        temperature=linspace(0,0.5,1000);
        temp_gap{iiu} = (Tcfit{iiu}.a)*tanh((Tcfit{iiu}.b)*sqrt( 1-temperature./(Tcfit{iiu}.c)) );
        fprintf('Tc fit value for U = %g and lambda=%g is %g --> %g K\n',U,lambda_m,Tcfit{iiu}.c,Tcfit{iiu}.c*conv)
        fprintf('Gap/Tc =%g\n',((maxgapm(iiu)/75)/tclist(iiu)))
        plot(temperature,real(temp_gap{iiu}),'--b','LineWidth',2,'DisplayName',['$T_\textrm{c} \approx$' num2str((Tcfit{iiu}.c)*conv) 'K'])
        plot(Tlist,abs(sgapvalplus), 'Marker','o','Color','#007BA7','LineWidth',2,'DisplayName',dispnm1);
        plot(Tlist,abs(sgapvalminus), 'Marker','x','Color','#007BA7','LineWidth',2,'DisplayName',dispnm2);
        text(10/conv,2,['$\lambda_m = $' num2str(lambda_m)],'FontSize',35);
        xlabel('T[K]','Fontsize',35);
        ylabel('$\mid\Delta_{\pm}(\textbf{k}^\textrm{max}_{\pm},$i$\pi$T$)\mid$ [meV]','Fontsize',35);
        axis([0.0 120/conv 0 35]);
        set(gca,'XTick',[0,10/conv,20/conv,30/conv,40/conv,50/conv,60/conv,70/conv,80/conv,90/conv,100/conv,110/conv,120/conv,130/conv,140/conv,150/conv],'YTick',[0:5:50],'FontSize',35,'FontName','Times');
        xticklabels({'0','10','20','30','40','50','60','70','80','90','100','110','120','130','140','150'});
        legend('FontSize',28,'FontName','Latin Modern Roman','Location','northeast')
        lgd = legend;
        lgd.NumColumns = 1;
        legend('boxoff')
        box on;
        
        
        figure(gg);
    elseif iiu == 6
        temp = Tlist;
        temp(numel(Tlist)+1) = tclist(iiu);
        sgapp = sgapvalplus;
        sgapm = sgapvalminus;
        sgapp(numel(sgapvalplus)+1) = 0;
        sgapm(numel(sgapvalminus)+1) = 0;
        p2 = plot(temp, sgapp,'LineStyle','-','Marker','none','MarkerSize',2,'Color','#DC143C','LineWidth',3,'DisplayName',dispnm);
        p1 = plot(temp, sgapm,'LineStyle','-','Marker','none','MarkerSize',2,'Color','#DC143C','LineWidth',3,'DisplayName',dispnm);
        set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        set(p1, 'markerfacecolor', get(p1, 'color'));
        set(p2, 'markerfacecolor', get(p2, 'color'));
        maxgapp(iiu) = sgapvalplus(1);
        maxgapm(iiu) = sgapvalminus(1);
        fprintf('max value of Delta minus =%g [meV]\n',maxgapm(iiu))
        
        figure(ss);hold on;
        set(gcf,'color','white','Position',get(0,'Screensize'));
        Tcfit{iiu} = fit(Tlist,transpose(sgapvalminus),myfit,'StartPoint',[sgapvalminus(1),10,neartclist(iiu)]);
        temperature=linspace(0,0.5,1000);
        temp_gap{iiu} = (Tcfit{iiu}.a)*tanh((Tcfit{iiu}.b)*sqrt( 1-temperature./(Tcfit{iiu}.c)) );
        fprintf('Tc fit value for U = %g and lambda=%g is %g --> %g K\n',U,lambda_m,Tcfit{iiu}.c,Tcfit{iiu}.c*conv)
        fprintf('Gap/Tc =%g\n',((maxgapm(iiu)/75)/tclist(iiu)))
        plot(temperature,real(temp_gap{iiu}),'--b','LineWidth',2,'DisplayName',['$T_\textrm{c} \approx$' num2str((Tcfit{iiu}.c)*conv) 'K'])
        plot(Tlist,abs(sgapvalplus), 'Marker','o','Color','#3D0C02','LineWidth',2,'DisplayName',dispnm1);
        plot(Tlist,abs(sgapvalminus), 'Marker','x','Color','#3D0C02','LineWidth',2,'DisplayName',dispnm2);
        text(10/conv,2,['$\lambda_m = $' num2str(lambda_m)],'FontSize',35);
        xlabel('T[K]','Fontsize',35);
        ylabel('$\mid\Delta_{\pm}(\textbf{k}^\textrm{max}_{\pm},$i$\pi$T$)\mid$ [meV]','Fontsize',35);
        axis([0.0 120/conv 0 35]);
        set(gca,'XTick',[0,10/conv,20/conv,30/conv,40/conv,50/conv,60/conv,70/conv,80/conv,90/conv,100/conv,110/conv,120/conv,130/conv,140/conv,150/conv],'YTick',[0:5:50],'FontSize',35,'FontName','Times');
        xticklabels({'0','10','20','30','40','50','60','70','80','90','100','110','120','130','140','150'});
        legend('FontSize',28,'FontName','Latin Modern Roman','Location','northeast')
        lgd = legend;
        lgd.NumColumns = 1;
        legend('boxoff')
        box on;
        
        
        figure(gg);
    end%end if statement
    
    %%%%%%%%%%
    %FIGURE 4%
    %%%%%%%%%%
    
    figure(gg);
    xlabel('T [K]','FontSize',30);
    ylabel('$\Delta_{\pm}(\textbf{k}^\textrm{max}_{\pm},$i$\pi$T$)$ [meV]','Fontsize',35);
    axis([0.0 140/conv -35 25]);
    text(100/conv,15,'$\Delta_{-}$','FontSize',35);
    text(10/conv,-15,'$\Delta_{+}$','FontSize',35);
    yp = yline(0,'k','LineWidth',2,'DisplayName','U=7t');
    set(get(get(yp,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(gca,'XTick',[0,10/conv,20/conv,30/conv,40/conv,50/conv,60/conv,70/conv,80/conv,90/conv,100/conv,110/conv,120/conv,130/conv,140/conv,150/conv,160/conv],'YTick',[-35:1:35],'FontSize',35,'FontName','Latin Modern Roman');
    xticklabels({'0',' ','20',' ','40',' ','60',' ','80',' ','100',' ','120',' ','140',' ','160'});
    yticklabels({'-35',' ',' ',' ',' ','-30',' ',' ',' ',' ','-25',' ',' ',' ',' ','-20',' ',' ',' ',' ','-15',' ',' ',' ',' ','-10',' ',' ',' ',' ','-5',' ',' ',' ',' ','0',...
        ' ',' ',' ',' ','5',' ',' ',' ',' ','10',' ',' ',' ',' ','15',' ',' ',' ',' ','20',' ',' ',' ',' ','25',' ',' ',' ',' ','30',' ',' ',' ',' ','35'});
    legend('FontSize',35,'FontName','Latin Modern Roman','Location','southeast')
    lgd = legend;
    lgd.NumColumns = 1;
    legend('boxoff')
    box on;
end%end lambda loop