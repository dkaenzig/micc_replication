

%%% figures for welfare effects and SCC

% graphical parameters
mdefblue  = [0,      0.4470, 0.7410] ;
mdefred   = [0.8500, 0.3250, 0.0980] ;
mdefgreen = [ 53  98  26  ] / 255 ;


% define indices for baseline cases
    % with expectation
irho   = (1:Nrho)' ;
iT     = ( (Nrho+1):(Nrho+NTG) )' ;
iCS    = ( (Nrho+NTG+1):(Nrho+NTG+NCS) )' ;

    % without expectation (surprises)
irho_noexp   = Nparam0 + (1:Nrho)' ;
iT_noexp     = Nparam0 + ( (Nrho+1):(Nrho+NTG) )' ;    
iCS_noexp    = Nparam0 + ( (Nrho+NTG+1):(Nrho+NTG+NCS) )' ;

    % BU data
irho_BU   = 2*Nparam0 + (1:Nrho)' ;
iT_BU     = 2*Nparam0 + ( (Nrho+1):(Nrho+NTG) )' ;    
iCS_BU    = 2*Nparam0 + ( (Nrho+NTG+1):(Nrho+NTG+NCS) )' ;

% graphical parameters
fontsize = 0.95*p.fontsize ;
lw = 4 ;
xlimw = -80 ;

i = 40 ;
xp = 305 ;
yp = 270 ;

% welfare vs. discount rate
i = i + 1 ;
figure(i)
set(gcf,'Position',[100 100 xp yp])

plot(rhoG,100*welfareG(irho,1),'LineWidth',lw); hold on;
plot(rhoG,100*welfareG(irho_BU,1),'-.','LineWidth',lw,'Color',mdefgreen); hold on;
plot(rhoG,100*welfareG(irho_noexp,1),':','LineWidth',lw,'Color',mdefblue); hold on;
plot(rhoG,100*welfareG(irho,2),'--','LineWidth',lw,'Color',mdefred); hold on
plot(0.02,100*GlobalWelfare, 'o', 'MarkerSize', 3*lw, 'MarkerEdgeColor', mdefblue, 'MarkerFaceColor', mdefblue, 'LineWidth', 2);
ylim([xlimw,0])
xlim([min(rhoG),max(rhoG)])
xlabel('Rate of time pref. ($\rho$)','interpreter','latex')
ylabel('Percent','interpreter','latex')
legend('Global $T$, PWT', ...
       'Global $T$, BU', ...
       'Global $T$, unexp.', ...
       'Local $T$','interpreter','latex','Location','SouthEast')
legend boxoff
ax = gca ; ax.XTick = [0.005 0.01 0.02 0.03 0.04 0.05 ] ;
AxisFonts(fontsize,fontsize,fontsize,fontsize,0.85*fontsize,fontsize)
set(gca,'TickLabelInterpreter','latex')

if optns.SaveFigures
    filename = "figures/" + optns.Date + "/WelfareSCC_W_rho_" + optns.Date ;
    printPDF

end

hold off


% welfare vs. warming
i = i + 1 ;
figure(i)
set(gcf,'Position',[100 100 xp yp])

plot(T2100G,100*welfareG(iT,1),'LineWidth',lw); hold on
plot(T2100G,100*welfareG(iT_BU,1),'-.','LineWidth',lw,'Color',mdefgreen); hold on
plot(T2100G,100*welfareG(iT_noexp,1),':','LineWidth',lw,'Color',mdefblue); hold on
plot(T2100G,100*welfareG(iT,2),'--','LineWidth',lw,'Color',mdefred); hold on
plot(3,100*GlobalWelfare, 'o', 'MarkerSize', 3*lw, 'MarkerEdgeColor', mdefblue, 'MarkerFaceColor', mdefblue, 'LineWidth', 2);
ylim([xlimw,0])
xlim([min(T2100G(:)),max(T2100G(:))])
xlabel('2100 temp. vs. pre-ind. ($^\circ$C)','interpreter','latex') ;
ylabel('Percent','interpreter','latex')
ax = gca ; ax.XTick = [ 1.5 2 3 4 5 6 ] ;
AxisFonts(fontsize,fontsize,fontsize,fontsize,fontsize,fontsize)
legend off
set(gca,'TickLabelInterpreter','latex')

if optns.SaveFigures
    filename = "figures/" + optns.Date + "/WelfareSCC_W_T_" + optns.Date ;
    printPDF

end

hold off


% welfare vs. climate sensitivity
i = i + 1 ;
figure(i)
set(gcf,'Position',[100 100 xp yp])

plot(CSG,100*welfareG(iCS,1),'LineWidth',lw); hold on
plot(CSG,100*welfareG(iCS_BU,1),'-.','LineWidth',lw,'Color',mdefgreen); hold on
plot(CSG,100*welfareG(iCS_noexp,1),':','LineWidth',lw,'Color',mdefblue); hold on
plot(CSG,100*welfareG(iCS,2),'--','LineWidth',lw,'Color',mdefred); hold on
plot(adjBaseline,100*GlobalWelfare, 'o', 'MarkerSize', 3*lw, 'MarkerEdgeColor', mdefblue, 'MarkerFaceColor', mdefblue, 'LineWidth', 2);
ylim([xlimw,0])
xlim([min(CSG(:)),max(CSG(:))])
xlabel('Climate sens. vs. median','interpreter','latex')
ylabel('Percent','interpreter','latex')
ax = gca ; ax.XTick = [ 0.5 1 1.5 2 ] ;
AxisFonts(fontsize,fontsize,fontsize,fontsize,fontsize,fontsize)
legend off
set(gca,'TickLabelInterpreter','latex')

if optns.SaveFigures
    filename = "figures/" + optns.Date + "/WelfareSCC_W_sigma_" + optns.Date ;
    printPDF

end

hold off


% SCC vs. discount rate
i = i + 1 ;
figure(i)
set(gcf,'Position',[100 100 xp yp])

plot(rhoG,SCCG(irho,1),'LineWidth',lw); hold on;
plot(rhoG,SCCG(irho_BU,1),'-.','LineWidth',lw,'Color',mdefgreen); hold on;
plot(rhoG,SCCG(irho_noexp,1),':','LineWidth',lw,'Color',mdefblue); hold on;
plot(rhoG,SCCG(irho,2),'--','LineWidth',lw,'Color',mdefred); hold on
plot(0.02,GlobalSCC, 'o', 'MarkerSize', 3*lw, 'MarkerEdgeColor', mdefblue, 'MarkerFaceColor', mdefblue, 'LineWidth', 2);
ylim([0,ylimSCC])
xlim([min(rhoG),max(rhoG)])
xlabel('Rate of time pref.  ($\rho$)','interpreter','latex')
ylabel('Dollars','interpreter','latex')
ax = gca ; ax.XTick = [0.005 0.01 0.02 0.03 0.04 0.05 ] ;
AxisFonts(fontsize,fontsize,fontsize,fontsize,fontsize,fontsize)
legend off
set(gca,'TickLabelInterpreter','latex')

if optns.SaveFigures
    filename = "figures/" + optns.Date + "/WelfareSCC_SCC_rho_" + optns.Date ;
    printPDF

end

hold off


% SCC vs. temperature
i = i + 1 ;
figure(i)
set(gcf,'Position',[100 100 xp yp])

plot(T2100G,SCCG(iT,1),'LineWidth',lw); hold on
plot(T2100G,SCCG(iT_BU,1),'-.','LineWidth',lw,'Color',mdefgreen); hold on
plot(T2100G,SCCG(iT_noexp,1),':','LineWidth',lw,'Color',mdefblue); hold on
plot(T2100G,SCCG(iT,2),'--','LineWidth',lw,'Color',mdefred); hold on
plot(3,GlobalSCC, 'o', 'MarkerSize', 3*lw, 'MarkerEdgeColor', mdefblue, 'MarkerFaceColor', mdefblue, 'LineWidth', 2);
ylim([0,ylimSCC])
xlim([min(T2100G(:)),max(T2100G(:))])
xlabel('2100 temp. vs. pre-ind. ($^\circ$C)','interpreter','latex') ;
ylabel('Dollars','interpreter','latex')
ax = gca ; ax.XTick = [ 1.5 2 3 4 5 6 ] ;
AxisFonts(fontsize,fontsize,fontsize,fontsize,fontsize,fontsize)
legend off
set(gca,'TickLabelInterpreter','latex')

if optns.SaveFigures
    filename = "figures/" + optns.Date + "/WelfareSCC_SCC_T_" + optns.Date ;
    printPDF

end

hold off


% welfare vs. climate sensitivity
i = i + 1 ;
figure(i)
set(gcf,'Position',[100 100 xp yp])

plot(CSG,SCCG(iCS,1),'LineWidth',lw); hold on
plot(CSG,SCCG(iCS_BU,1),'-.','LineWidth',lw,'Color',mdefgreen); hold on
plot(CSG,SCCG(iCS_noexp,1),':','LineWidth',lw,'Color',mdefblue); hold on
plot(CSG,SCCG(iCS,2),'--','LineWidth',lw,'Color',mdefred); hold on
plot(adjBaseline,GlobalSCC, 'o', 'MarkerSize', 3*lw, 'MarkerEdgeColor', mdefblue, 'MarkerFaceColor', mdefblue, 'LineWidth', 2);
ylim([0,ylimSCC])
xlim([min(CSG(:)),max(CSG(:))])
xlabel('Climate sens. vs. median','interpreter','latex')
ylabel('Dollars','interpreter','latex')
ax = gca ; ax.XTick = [ 0.5 1 1.5  2 ] ;
AxisFonts(fontsize,fontsize,fontsize,fontsize,fontsize,fontsize)
legend off
set(gca,'TickLabelInterpreter','latex')

if optns.SaveFigures
    filename = "figures/" + optns.Date + "/WelfareSCC_SCC_sigma_" + optns.Date ;
    printPDF

end

hold off




