close all
t = datenum(1995,7,27):1:datenum(1995,7,27)+length(y1)-1;
y1 = sperobj(1).counts;

ax1=subplot(4,1,4),plot(t,y1,'k')
datetick('x','keeplimits')
set(ax1,'YLim',[0 2000],'XLim',[min(t) max(t)]) % max is 3221
hold on
ylabel('VT')

y2 = sperobj(2).counts;
ax2=subplot(4,1,3),plot(t,y2,'k')
datetick('x','keeplimits')
set(ax2,'YLim',[0 1200],'XLim',[min(t) max(t)]) % max is 4089
ylabel('HYBRID')

y3 = sperobj(3).counts;
ax3=subplot(4,1,2),plot(t,y3,'k')
datetick('x','keeplimits')
ylabel('LONG-PERIOD')
set(ax3,'XLim',[min(t) max(t)])

y4 = sperobj(4).counts;
ax4=subplot(4,1,1),plot(t,y4,'k')
ylabel('ROCKFALL')
datetick('x','keeplimits')
set(ax4,'XLim',[min(t) max(t)])

nansum(y1)+nansum(y2)+nansum(y3)+nansum(y4)
%%
figure
pt = datenum(1995,7,27):1:datenum(1995,7,27)+length(y1)-1;
y1 = sperobj(1).counts;
datestr(min(t))
datestr(max(t))
y1(isnan(y1))=0;
plot(t,cumsum(y1),'b','LineWidth',2)
hold on
y2(isnan(y2))=0;
y3(isnan(y3))=0;
y4(isnan(y4))=0;
plot(t,cumsum(y2),'r','LineWidth',2)
plot(t,cumsum(y3),'g','LineWidth',2)
plot(t,cumsum(y4),'k','LineWidth',2)
legend({'VT';'HYBRID';'LP';'ROCKFALL'},'Location','northwest')
datetick('x')
set(gca,'XLim',[min(t) max(t)])
ylabel('Cumulative event counts')
print -dpng ~/Desktop/data_mastering_paper/figures/ASNE_EVENT_COUNTS.PNG
