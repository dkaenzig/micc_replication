

%%% plots and exports damage functions and transitional dynamics region by region

% graphical parameters
mdefblue = [0, 0.4470, 0.7410] ;
mdefred  = [0.8500 0.3250 0.0980] ;

% graphical inputs
tmax = 76 + 50 ; 
Tmax = 2024 + tmax ;
fontsize = 14 ; 
ts = 2024 + p.times ;
xticks = [2024 2050 2075 2100 2125 2150] ;
lw = 4 ;
nx = 3 ;
ny = 3 ;

% plot warming transition by looping over regions
figure()
set(gcf,'Position',[100 100 1000 1000])

for i = 1:9

    is   = num2str(i) ;
    ri   = "r"+is ;

    name = regionnames.group_name(i) ;
        
    subplot(nx,ny,i)
    plot(ts,100*irf.(ri).global.Yp,'LineWidth',lw); hold on;
    yline(0,'k-','LineWidth',0.5)
    title(name{1},'interpreter','latex')
    xlabel('Years','interpreter','latex')
    ylabel("Percent",'interpreter','latex')
    xlim([2024,Tmax])
    ylim([-80,40])
    ax = gca ; ax.XTick = xticks ;
    AxisFonts(fontsize,fontsize,fontsize,fontsize,0.9*fontsize,fontsize)
    set(gca,'TickLabelInterpreter','latex')
    legend off
    hold off

end


if optns.SaveFigures
    filename = "figures/" + optns.Date + "/TransitionRegions" + optns.Date ;
    printPDF
end

hold off




% plot damage functions by looping over regions
figure()
set(gcf,'Position',[100 100 1000 1000])

for i = 1:9

    is   = num2str(i) ;
    ri   = "r"+is ;

    name = regionnames.group_name(i) ;
        
    subplot(nx,ny,i)
    plot(p.times,100*fit.(ri).global.zetas,'LineWidth',lw); hold on;
    yline(0,'k-','LineWidth',0.5)
    title(name{1},'interpreter','latex')
    xlabel('Years','interpreter','latex')
    ylabel("Percent",'interpreter','latex')
    xlim([0,10])
    ylim([-10,5])
    AxisFonts(fontsize,fontsize,fontsize,fontsize,0.9*fontsize,fontsize)
    set(gca,'TickLabelInterpreter','latex')
    legend off
    hold off

end

if optns.SaveFigures
    filename = "figures/" + optns.Date + "/DamageFunctionRegions" + optns.Date ;
    printPDF
end

hold off


%%% export damage functions, its parameters, and counterfactuals

% initializes
damage_params = table() ;
damage_func   = table() ;
Yp_path       = table() ;

% world
damage_params.("Parameter") = ["A";"B";"C"] ;
xs    = fit.global.xs ; 
xs(3) = xs(2) + xs(3) ;
damage_params.("World") = xs ;

hs = (0:1:20)' ;
damage_func.("Horizon") = hs ;
damage_func.("World")   = fit.global.zetas(1+hs) ;

ts = (1:1:tmax)' ;
Yp_path.("Years") = 2025+ts ;
Yp_path.("World") = irf.global.Yp(ts) ;


% regions
for i = 1:9

    is    = num2str(i) ;
    ri    = "r"+is ;
    xs    = fit.(ri).global.xs ;
    xs(3) = xs(2) + xs(3) ;
    name  = regionnames.group_name(i) ;
    name  = name{1} ;

    damage_params.(name) = xs ;
    damage_func.(name)   = fit.(ri).global.zetas(1+hs) ;
    Yp_path.(name)       = irf.(ri).global.Yp(ts) ;

    if optns.SaveFigures
        writetable(damage_params, "tables/" + optns.Date + "/DamageParamers"  + optns.Date + ".csv") ;
        writetable(damage_func,   "tables/" + optns.Date + "/DamageFunctions" + optns.Date+ ".csv") ;
        writetable(Yp_path,       "tables/" + optns.Date + "/OutputPath"      + optns.Date+ ".csv") ;
    end


end