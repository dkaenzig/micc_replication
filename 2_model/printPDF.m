
% prints current figure to pdf with right size
% h is current figure handle (get with h = gcf ; before running function)
% filename must be a string defined before running the script
% cannot make it into a function because of how gcf works

h = gcf ;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(h,filename,'-dpdf','-r0')
