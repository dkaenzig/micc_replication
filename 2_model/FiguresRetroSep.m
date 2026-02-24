


%%% retrospective analysis

% graphical parameters
xp = 320 ;
yp = 300 ;
i  = 30 ;

mdefblue = [0, 0.4470, 0.7410] ;
mdefred  = [0.8500 0.3250 0.0980] ;
fontsize = p.fontsize ;
lw = 4 ;
xticks = [ 1960 1980 2000 2019 ] ;


% plot results

i = i + 1 ;
figure(i)
set(gcf,'Position',[100 100 xp yp])

plot(p.times((1:length(Ydata)-1))+1960,GYdata,'LineWidth',lw) ; hold on ;
plot(p.times((1:length(Ydata)-1))+1960,GYpotential,'--','Color',mdefred,'LineWidth',lw); hold on;
yline(mean(GYdata),'Color',mdefblue,'LineWidth',1.5*lw) ; hold on;
yline(mean(GYpotential),'--','Color',mdefred,'LineWidth',1.5*lw)
xlim([1960,2019])
ylim([-0.05,0.10])
ax = gca ; ax.XTick = xticks ;
AxisFonts(fontsize,fontsize,fontsize,fontsize,0.9*fontsize,fontsize)
legend('With climate change','Without climate change',...
       'interpreter','latex','Location','NorthWest')
legend boxoff
set(gca,'TickLabelInterpreter','latex')

if optns.SaveFigures
    filename = "figures/" + optns.Date + "/Retro_g" + optns.Date ;
    printPDF
end

hold off

i = i + 1 ;
figure(i)
set(gcf,'Position',[100 100 xp yp])

plot(p.times((1:length(Ydata)-1))+1960,fractionLost,'LineWidth',lw); hold on;
yline(mean(fractionLost),'Color',mdefblue,'LineWidth',1.5*lw); hold on;
plot(years,fitted,':','Color',mdefblue,'LineWidth',1.5*lw); hold on;
xlim([1960,2019])
ylim([-0.6,0.3])
ax = gca ; ax.XTick = xticks ;
AxisFonts(fontsize,fontsize,fontsize,fontsize,0.9*fontsize,fontsize)
legend off
set(gca,'TickLabelInterpreter','latex')

if optns.SaveFigures
    filename = "figures/" + optns.Date + "/Retro_gLost" + optns.Date ;
    printPDF
end

hold off

i = i + 1 ;
figure(i)
set(gcf,'Position',[100 100 xp yp])

plot(p.times(1:length(Ydata))+1960,Ydata/Ydata(1),'LineWidth',lw) ; hold on ;
plot(p.times(1:length(Ydata))+1960,Ypotential/Ypotential(1),'--','Color',mdefred,'LineWidth',lw)
xlim([1960,2019])
ax = gca ; ax.XTick = xticks ;
AxisFonts(fontsize,fontsize,fontsize,fontsize,0.9*fontsize,fontsize)
legend('With climate change','Without climate change',...
       'interpreter','latex','Location','NorthWest')
legend boxoff
set(gca,'TickLabelInterpreter','latex')

if optns.SaveFigures
    filename = "figures/" + optns.Date + "/Retro_Ylevel" + optns.Date ;
    printPDF
end

hold off
