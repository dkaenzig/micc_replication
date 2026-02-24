

%%% model fit with confidence intervals

% transparency of CIs
alpha90 = 0.15 ;
alpha95 = 0.075 ;

% choice of CIs to plot
up1 = "up90" ;
lo1 = "lo90" ;
up2 = "up95" ;
lo2 = "lo95" ;

% y axis limits
y_zeta.global = 100*[-0.08,0.02] ;
y_gdp.global  = 100*[-0.20,0.06] ;
y_cap.global  = 100*[-0.15,0.06] ;
y_tem.global  = [0,1] ;

y_zeta.local = 100*[-0.02,0.01] ;
y_gdp.local  = 100*[-0.03,0.02] ;
y_cap.local  = 100*[-0.02,0.01] ;
y_tem.local  = [0,1] ;

% other graphical parameters
fontsize      = p.fontsize ;
lw            = 4 ;
xsize         = 450 ;
ysize         = 300 ;
colors.global = mdefblue ;
colors.local  = mdefred  ;
horizon       = Data.PWT.global.horizon ;


for temp = [ "global" , "local"]
    

    % temperature
        
    figure()
    set(gcf,'Position',[100 100 xsize ysize])
    plot(  horizon, Data.PWT.(temp).temperature, '--', 'LineWidth',lw, 'Color', colors.(temp) ); hold on;
    
    xlim([0,10])
    ylim(y_tem.(temp))
    xlabel('Years','interpreter','latex')
    ylabel('$^\circ$C','interpreter','latex')
    AxisFonts(fontsize,fontsize,fontsize,fontsize,fontsize,fontsize)
    set(gca,'TickLabelInterpreter','latex')
    legend off
    
    if optns.SaveFigures
        filename = "figures/" + optns.Date + "/Est"+temp+"_Temp_OnlyProdCI_" + optns.Date ;
        printPDF
    end
    
    hold off
    
    
    
    % productivity
     
    figure()
    set(gcf,'Position',[100 100 xsize ysize])
    
    plot(  p.times,zeros(p.Nt,1),'k-','LineWidth',0.5); hold on;
    plot(  p.times, 100*fit.(temp).zetas, 'LineWidth',lw,'Color',colors.(temp)); hold on;
    plotCI(p.times, 100*CI.(temp).zetas.(up1), 100*CI.(temp).zetas.(lo1), colors.(temp), alpha90) ;
    plotCI(p.times, 100*CI.(temp).zetas.(up2), 100*CI.(temp).zetas.(lo2), colors.(temp), alpha95) ;
    xlim([0,10])
    ylim(y_zeta.(temp))
    xlabel('Years','interpreter','latex')
    ylabel('Percent','interpreter','latex')
    AxisFonts(fontsize,fontsize,fontsize,fontsize,fontsize,fontsize)
    set(gca,'TickLabelInterpreter','latex')
    legend off
    
    if optns.SaveFigures
        filename = "figures/" + optns.Date + "/Est"+temp+"_Zeta_OnlyProdCI_" + optns.Date ;
        printPDF
    end
    
    hold off


    
    
    
    % output
        
    figure()
    set(gcf,'Position',[100 100 xsize ysize])

    plot(  horizon, Data.PWT.(temp).gdp ,'--', 'LineWidth',lw, 'Color', colors.(temp)); hold on;
    plot(  p.times, fit.(temp).Yp_fit,         'LineWidth',lw, 'Color', colors.(temp)); hold on;
    plot(  p.times, zeros(p.Nt,1),       'k-', 'LineWidth',0.5); hold on;
    plotCI(horizon, CI.(temp).gdp.(up1), CI.(temp).gdp.(lo1), colors.(temp), alpha90) ;
    plotCI(horizon, CI.(temp).gdp.(up2), CI.(temp).gdp.(lo2), colors.(temp), alpha95) ;
    legend('Data','Model','interpreter','latex','Location','NorthWest')
    legend boxoff
    xlim([0,10])
    ylim(y_gdp.(temp))
    xlabel('Years','interpreter','latex')
    ylabel('Percent','interpreter','latex')
    AxisFonts(fontsize,fontsize,fontsize,fontsize,fontsize,fontsize)
    set(gca,'TickLabelInterpreter','latex')
    
    if optns.SaveFigures
        filename = "figures/" + optns.Date + "/Est"+temp+"_GDP_OnlyProdCI_" + optns.Date ;
        printPDF
    end
    
    hold off



    % capital
        
    figure()
    set(gcf,'Position',[100 100 xsize ysize])
    
    plot(  horizon, Data.PWT.(temp).capital, '--', 'LineWidth',lw, 'Color', colors.(temp)); hold on;
    plot(  p.times, fit.(temp).Kp_fit,             'LineWidth',lw, 'Color', colors.(temp)); hold on;
    plot(  p.times,      zeros(p.Nt,1),      'k-', 'LineWidth',0.5); hold on;
    plotCI(horizon, CI.(temp).capital.(up1), CI.(temp).capital.(lo1), colors.(temp), alpha90) ;
    plotCI(horizon, CI.(temp).capital.(up2), CI.(temp).capital.(lo2), colors.(temp), alpha95) ;
    legend('Data','Model','interpreter','latex','Location','SouthWest')
    legend boxoff
    xlim([0,10])
    ylim(y_cap.(temp))
    xlabel('Years','interpreter','latex')
    ylabel('Percent','interpreter','latex')
    AxisFonts(fontsize,fontsize,fontsize,fontsize,fontsize,fontsize)
    set(gca,'TickLabelInterpreter','latex')
    
    if optns.SaveFigures
        filename = "figures/" + optns.Date + "/Est"+temp+"_Cap_OnlyProdCI_" + optns.Date ;
        printPDF
    end
    
    hold off
    
        
end

