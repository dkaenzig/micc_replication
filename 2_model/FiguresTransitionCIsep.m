

%%% plot transitional dynamics with confidence bands

% graphical parameters
mdefblue = [0, 0.4470, 0.7410] ;
mdefred  = [0.8500 0.3250 0.0980] ;
alpha90 = 0.15 ;
alpha95 = 0.075 ;

% choice of CIs to plot
up1 = "up90" ;
lo1 = "lo90" ;
up2 = "up95" ;
lo2 = "lo95" ;

% figure counter
i = 20 ;

% common figure info
tmax = 76 + 50 ;
Tmax = 2024 + tmax ;
fontsize = p.fontsize ; 
ts = 2024 + p.times ;
xticks = [2024 2050 2075 2100 2125 2150] ;
lw = 4 ;
nx = 2 ;
ny = 3 ;

xp = 293 ;
yp = 234 ;


% temperature

i = i + 1 ;
figure(i)
set(gcf,'Position',[100 100 xp yp+23])

plot(ts,irf.global.TD,'LineWidth',lw); hold on;
plot(ts,irf.local.TD,'--','LineWidth',lw,'Color',mdefred); hold on;
ylabel('$^\circ$C','interpreter','latex')
xlabel('Years','interpreter','latex')
xlim([2024,Tmax])
ylim([0,2.2])
ax = gca ; ax.XTick = xticks ;
AxisFonts(fontsize,fontsize,fontsize,fontsize,0.9*fontsize,fontsize)
legend('Global $T$',...
       'Local $T$','interpreter','latex','Location','SouthEast')
legend boxoff
set(gca,'TickLabelInterpreter','latex')

if optns.SaveFigures
    filename = "figures/" + optns.Date + "/TransitionGlobalCI_T" + "_" + optns.Date ;
    printPDF
end

hold off



% output

i = i + 1 ;
figure(i)
set(gcf,'Position',[100 100 xp yp])

plot(ts,zeros(p.Nt,1),'k-','LineWidth',0.5); hold on;

plot(  ts, 100*irf.global.Yp, 'LineWidth',lw, 'Color', mdefblue); hold on;
plotCI(ts, 100*CI.global.Yp.(up1), 100*CI.global.Yp.(lo1), mdefblue, alpha90) ;
plotCI(ts, 100*CI.global.Yp.(up2), 100*CI.global.Yp.(lo2), mdefblue, alpha95) ;

plot(  ts, 100*irf.local.Yp, '--','LineWidth',lw,'Color',mdefred); hold on;
plotCI(ts, 100*CI.local.Yp.(up1), 100*CI.local.Yp.(lo1), mdefred, alpha90) ;
plotCI(ts, 100*CI.local.Yp.(up2), 100*CI.local.Yp.(lo2), mdefred, alpha95) ;

%title('(b) Output','interpreter','latex')
ylabel('Percent','interpreter','latex')
xlim([2024,Tmax])
ylim([-70,10])
ax = gca ; ax.XTick = xticks ;
AxisFonts(fontsize,fontsize,fontsize,fontsize,fontsize,fontsize)
legend off
set(gca,'TickLabelInterpreter','latex')

if optns.SaveFigures
    filename = "figures/" + optns.Date + "/TransitionGlobalCI_Y" + "_" + optns.Date ;
    printPDF
end

hold off


% capital

i = i + 1 ;
figure(i)
set(gcf,'Position',[100 100 xp yp])

plot(ts,zeros(p.Nt,1),'k-','LineWidth',0.5); hold on;

plot(  ts, 100*irf.global.Kp, 'LineWidth', lw, 'Color', mdefblue); hold on;
plotCI(ts, 100*CI.global.Kp.(up1), 100*CI.global.Kp.(lo1), mdefblue, alpha90) ;
plotCI(ts, 100*CI.global.Kp.(up2), 100*CI.global.Kp.(lo2), mdefblue, alpha95) ;

plot(  ts, 100*irf.local.Kp, '--','LineWidth', lw, 'Color', mdefred); hold on;
plotCI(ts, 100*CI.local.Kp.(up1), 100*CI.local.Kp.(lo1), mdefred, alpha90) ;
plotCI(ts, 100*CI.local.Kp.(up2), 100*CI.local.Kp.(lo2), mdefred, alpha95) ;

ylabel('Percent','interpreter','latex')
xlim([2024,Tmax])
ylim([-70,10])
ax = gca ; ax.XTick = xticks ;
AxisFonts(fontsize,fontsize,fontsize,fontsize,fontsize,fontsize)
legend off
set(gca,'TickLabelInterpreter','latex')

if optns.SaveFigures
    filename = "figures/" + optns.Date + "/TransitionGlobalCI_K" + "_" + optns.Date ;
    printPDF
end

hold off



% consumption

i = i + 1 ;
figure(i)
set(gcf,'Position',[100 100 xp yp])


plot(ts,zeros(p.Nt,1),'k-','LineWidth',0.5); hold on;

plot(  ts, 100*irf.global.Ep, 'LineWidth', lw,'Color', mdefblue); hold on;
plotCI(ts, 100*CI.global.Ep.(up1), 100*CI.global.Ep.(lo1), mdefblue, alpha90) ;
plotCI(ts, 100*CI.global.Ep.(up2), 100*CI.global.Ep.(lo2), mdefblue, alpha95) ;

plot(  ts, 100*irf.local.Ep, '--', 'LineWidth', lw, 'Color', mdefred); hold on;
plotCI(ts, 100*CI.local.Ep.(up1), 100*CI.local.Ep.(lo1), mdefred, alpha90) ;
plotCI(ts, 100*CI.local.Ep.(up2), 100*CI.local.Ep.(lo2), mdefred, alpha95) ;

%title('(d) Consumption','interpreter','latex')
ylabel('Percent','interpreter','latex')
xlim([2024,Tmax])
ylim([-70,10])
ax = gca ; ax.XTick = xticks ;
AxisFonts(fontsize,fontsize,fontsize,fontsize,fontsize,fontsize)
legend off
set(gca,'TickLabelInterpreter','latex')

if optns.SaveFigures
    filename = "figures/" + optns.Date + "/TransitionGlobalCI_C" + "_" + optns.Date ;
    printPDF
end

hold off



% welfare

i = i + 1 ;
figure(i)
set(gcf,'Position',[100 100 xp yp])

plot(ts,zeros(p.Nt,1),'k-','LineWidth',0.5); hold on;

plot(  ts, 100*irf.global.welfare, 'LineWidth', lw,'Color', mdefblue); hold on;
plotCI(ts, 100*CI.global.welfare.(up1), 100*CI.global.welfare.(lo1), mdefblue, alpha90) ;
plotCI(ts, 100*CI.global.welfare.(up2), 100*CI.global.welfare.(lo2), mdefblue, alpha95) ;

plot(   ts, 100*irf.local.welfare,'--','LineWidth',lw,'Color',mdefred); hold on;
plotCI(ts, 100*CI.local.welfare.(up1), 100*CI.local.welfare.(lo1), mdefred, alpha90) ;
plotCI(ts, 100*CI.local.welfare.(up2), 100*CI.local.welfare.(lo2), mdefred, alpha95) ;

ylabel('Percent','interpreter','latex')
xlim([2024,Tmax])
ylim([-70,10])
ax = gca ; ax.XTick = xticks ;
AxisFonts(fontsize,fontsize,fontsize,fontsize,fontsize,fontsize)
legend off
set(gca,'TickLabelInterpreter','latex')

% to make the graph plot below to be of same size as here
LI = ax.LooseInset; 

if optns.SaveFigures
    filename = "figures/" + optns.Date + "/TransitionGlobalCI_W" + "_" + optns.Date ;
    printPDF
end

hold off



% SCC

i = i + 1 ;
figure(i)
set(gcf,'Position',[100 100 1.045*xp 0.89*yp])

plot(ts,zeros(p.Nt,1),'k-','LineWidth',0.5); hold on;

plot(  ts, irf.global.SCC*ones(size(ts)),'LineWidth',lw,'Color', mdefblue); hold on;
plotCI(ts, CI.global.SCC.(lo1)*ones(size(ts))', CI.global.SCC.(up1)*ones(size(ts))', mdefblue, alpha90) ;
plotCI(ts, CI.global.SCC.(lo2)*ones(size(ts))', CI.global.SCC.(up2)*ones(size(ts))', mdefblue, alpha95) ;
    % switch order for SCC
    
plot(   ts, irf.local.SCC*ones(size(ts)),'--','LineWidth',lw,'Color', mdefred); hold on;
plotCI(ts, CI.local.SCC.(lo1)*ones(size(ts))', CI.local.SCC.(up1)*ones(size(ts))', mdefred, alpha90) ;
plotCI(ts, CI.local.SCC.(lo2)*ones(size(ts))', CI.local.SCC.(up2)*ones(size(ts))', mdefred, alpha95) ;
    % switch order for SCC

ylabel('Dollars','interpreter','latex')
xlim([2024,Tmax])
ylim([ylimSCC0,ylimSCC1])
ax = gca ; ax.XTick = {} ; 

AxisFonts(fontsize,fontsize,fontsize,fontsize,fontsize,fontsize)
legend off
set(gca,'TickLabelInterpreter','latex')

% to force the graph plot to be of same size as above
ax.LooseInset = LI;

if optns.SaveFigures
    filename = "figures/" + optns.Date + "/TransitionGlobalCI_SCC" + "_" + optns.Date ;
    printPDF
end

hold off

