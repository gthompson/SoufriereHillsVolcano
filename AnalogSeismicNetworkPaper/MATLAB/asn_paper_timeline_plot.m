
% make a sort of helicorder plot with 1 row = 1 year
startyear=1988.9;
endyear=2011;

extrusion_color = [.8 .8 .8];


 % define phases of extrusion
phase1.snum = datenum(1995,11,15);
phase1.enum=datenum(1998,3,10);
phase1.col = extrusion_color;
phase1.label='Phase 1';
phase1.plotas = 'patch';
phase1.level = 0;

phase2.snum = datenum(1999,11,27);
phase2.enum = datenum(2003,8,1);
phase2.col = extrusion_color;
phase2.label='Phase 2';
phase2.plotas = 'patch';
phase2.level = 0;

phase3.snum = datenum(2005,8,1);
phase3.enum = datenum(2007,4,20);
phase3.col = extrusion_color;
phase3.label='3';
phase3.plotas = 'patch';
phase3.level = 0;

phase4a.snum=datenum(2008,8,8);
phase4a.enum=datenum(2008,10,8);
phase4a.col= extrusion_color;
phase4a.label='';
phase4a.plotas = 'patch';
phase4a.level = 0;

phase4b.snum=datenum(2008,12,2);
phase4b.enum=datenum(2009,1,3);
phase4b.col= extrusion_color;
phase4b.label='';
phase4b.plotas = 'patch';
phase4b.level = 0;

phase5.snum = datenum(2009,10,8);
phase5.enum=datenum(2010,2,11);
phase5.col = extrusion_color;
phase5.label='5';
phase5.plotas = 'patch';
phase5.level = 0;


% define events
% volcanic/seismic/natural below line, network-related above line
event1.dnum = datenum(1995,7,18);
event1.label = sprintf('1st\nphreatic\nexplosion');
event1.pos = 'below';
event1.labelheight = 1;


event2.dnum = datenum(1995,7,27);
event2.label = sprintf('ASN &\nVDAP\nsystem\ninstalled');
event2.pos = 'above';
event2.labelheight = 1;


event3.dnum = datenum(1996,10,19);
event3.label = sprintf('DSN\ninstalled');
event3.pos = 'above';
event3.labelheight = 0.2;


event4.dnum = datenum(2001,3,2);
event4.label = sprintf('Seislog\nreplaced\nXdetect');
event4.pos = 'above';
event4.labelheight = 0.9;


event5.dnum = datenum(2004,12,16);
event5.label = sprintf('ASN\nphased\nout\n(DSN\nupgrade)');
event5.pos = 'above';
event5.labelheight = 0.9;


event6.dnum = datenum(2004,1,17);
event6.label = sprintf('last\nRSAM\ndata\ndownload');
event6.pos = 'above';
event6.labelheight = 0.2;

event7.dnum = (phase4a.snum + phase4a.enum)/2;
event7.label = '4a';
event7.pos = 'below';
event7.labelheight = 0.5;

event8.dnum = (phase4b.snum + phase4b.enum)/2;
event8.label = '4b';
event8.pos = 'below';
event8.labelheight = 0.2;

event9.dnum = datenum(1992,6,30);
event9.label = sprintf('GH, BE\n & SP\nrestored');
event9.pos = 'above';
event9.labelheight = 0.2;

event10.dnum = datenum(1989,9,17);
event10.label = sprintf('hurricane\nHugo');
event10.pos = 'below';
event10.labelheight = 0.2;

event11.dnum = datenum(1992,8,9);
event11.label = sprintf('swarm');
event11.pos = 'below';
event11.labelheight = 0.2;

event12.dnum = datenum(1994,2,3);
event12.label = sprintf('swarm');
event12.pos = 'below';
event12.labelheight = 0.2;

event13.dnum = datenum(1994,2,3);
event13.label = sprintf('WH\nadded');
event13.pos = 'above';
event13.labelheight = 0.2;

event14.dnum = datenum(1994,11,27);
event14.label = sprintf('major\nswarm');
event14.pos = 'below';
event14.labelheight = 0.4;

event15.dnum = datenum(1989,9,17);
event15.label = sprintf('local\nnetwork\ndestroyed');
event15.pos = 'above';
event15.labelheight = 0.2;


% makesinglelineplot
close all
figure
figbox=get(gcf,'Position');
set(gcf,'Position',[100 100 figbox(3)*1.8 figbox(4)]);
pheight = 0.13;
lspace = 0.02 + pheight; % space between line and 'arrow' line
startdnum=datenum(startyear,1,1); 
enddnum= datenum(endyear,1,1);
plot([startdnum enddnum], [0 0], 'k','LineWidth',3)
hold on

% add phases of extrusion
addsinglelinephase(phase1, pheight)
addsinglelinephase(phase2, pheight)
addsinglelinephase(phase3, pheight)
addsinglelinephase(phase4a, pheight)
addsinglelinephase(phase4b, pheight)
addsinglelinephase(phase5, pheight)

% add timeline events
addsinglelineevent(event1, pheight/2)
addsinglelineevent(event2, pheight/2)
addsinglelineevent(event3, lspace)
addsinglelineevent(event4, lspace)
addsinglelineevent(event5, pheight/2)
addsinglelineevent(event6, pheight/2)
addsinglelineevent(event7, lspace)
addsinglelineevent(event8, lspace)
addsinglelineevent(event9, pheight/2)
addsinglelineevent(event10, pheight/2)
addsinglelineevent(event11, pheight/2)
addsinglelineevent(event12, pheight/2)
addsinglelineevent(event13, pheight/2)
addsinglelineevent(event14, pheight/2)
addsinglelineevent(event15, pheight/2)
%patch([phase1.snum phase1.enum phase1.enum phase1.snum],[phase1.level-0.5 phase1.level-0.5 phase1.level+0.5 phase1.level+0.5],phase1.col)
set(gca,'YLim',[-1.8 1.8])
set(gca,'YTick',[])
set(gca,'XLim',[startdnum enddnum]);
datetick('x','keeplimits')
% set(gca,'XTick',[startyear:1:endyear]);
% a = get(gca,'XTickLabel');  
% set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
xtickangle(45)
set(gca,'FontSize',12)
hold off
grid on
%%
print -deps ~/src/AnalogSeismicNetworkPaper/FIGURES/asn_paper_timeline.eps
print -dpng ~/src/AnalogSeismicNetworkPaper/FIGURES/asn_paper_timeline.png
% cannot get a good PDF render of same figure
orient(gcf,'landscape')
print -dpdf -bestfit ~/src/AnalogSeismicNetworkPaper/FIGURES/asn_paper_timeline.pdf
%%
function addsinglelinephase(p,pheight)
    patch([p.snum p.enum p.enum p.snum],[p.level-pheight p.level-pheight p.level+pheight p.level+pheight],p.col)
    text(mean([p.snum p.enum]), p.level, p.label, 'HorizontalAlignment', 'center', ...
          'VerticalAlignment', 'middle');
end

function addsinglelineevent(e,lspace)
    
    if strcmp(e.pos,'above')
        %scatter(e.dnum, lspace, 'v')
        plot([e.dnum e.dnum],[lspace e.labelheight+lspace],'k','LineWidth',1.5)
        
        text(e.dnum, e.labelheight+lspace, e.label,  'HorizontalAlignment', 'center', ...
          'VerticalAlignment', 'bottom');
    else
        %scatter(e.dnum, -0.1, '^')
        plot([e.dnum e.dnum],[-lspace -e.labelheight-lspace],'k','LineWidth',1.5)
        text(e.dnum, -e.labelheight-lspace, e.label,  'HorizontalAlignment', 'center', ...
          'VerticalAlignment', 'top');
    end

end

% %%  multi-line plot
% close all
% figure
% makehelichronplot(startyear,endyear)
% makepatch(phase1)
% makepatch(phase2)
% makepatch(phase3)
% makepatch(phase4a)
% makepatch(phase4b)
% makepatch(phase5)
% addevent(event1);
% addevent(event2);
% addevent(event3);
% addevent(event4);
% addevent(event5);
% addevent(event6);
% addevent(event7);
% addevent(event8);
% set(gca,'YLim', [startyear-0.3  endyear-1+0.3])

% function addevent(thisp)
%     sy = dnum2year(thisp.dnum);
%     sdv = floor(sy);
%     hold on
%     if strcmp(thisp.pos,'above')
%         scatter(mod(sy,1),floor(sy)-0.1,70,'kv');
%         text(mod(sy,1),floor(sy)-0.25,thisp.label,  'HorizontalAlignment', 'center', ...
%       'VerticalAlignment', 'bottom');
%     else
%         scatter(mod(sy,1),floor(sy)+0.1,70,'k^');
%         text(mod(sy,1),floor(sy)+0.25,thisp.label,  'HorizontalAlignment', 'center', ...
%       'VerticalAlignment', 'top');   
%     end
%     hold off
% end
% 
% function makepatch(thisp)
%     sy = dnum2year(thisp.snum);
%     ey = dnum2year(thisp.enum);
%     sdv = floor(sy);
%     edv = floor(ey);
%     hold on
%     thisp.plotas
%     if strcmp(thisp.plotas,'patch')
%         %patch([thisp.snum thisp.enum thisp.enum thisp.snum thisp.snum], [thisp.level thisp.level thisp.level+1 thisp.level+1 thisp.level], thisp.col);
%         %datetick('x')
%         for thisy = sdv:edv
%             thissy = mod(max([sy thisy]), 1);
%             thisey = mod(min([ey thisy+1]), 1); 
%             lvl = thisp.level;
%             if thisey==0
%                 thisey=1;
%             end
%             patch([thissy thisey thisey thissy thissy], [thisy-0.1+lvl thisy-0.1+lvl thisy+0.1+lvl thisy+0.1+lvl thisy-0.1+lvl], thisp.col)
% 
%         end
%     else
%         scatter(mod(sy,1),floor(sy)-0.15,'kv');
%         scatter(mod(ey,1),floor(ey)-0.15,'kv');
%     end
%     hold off
% end
% 
% function makehelichronplot(starty,endy)
%     figure
%     subplot(1,1,1)
%     for y = starty:endy-1
%         line([0 1],[y y],'LineWidth',1,'Color','k');
%     end
%     set(gca,'YTick',starty:1:endy);
%     %, ['YTickLabel',[endy-1:-1:starty])
%     set(gca,'Ydir','reverse')
%     dnum = [];
%     for c=1:12
%         dnum = [dnum datenum(1995,c,1)];
%     end
%     xticks = dnum2year(dnum) - 1995;
%     xticklabels = {'J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'};
%     set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
% end
% 
% function y=dnum2year(dnum)
%     for c=1:length(dnum)
%         dv = datevec(dnum(c));
%         y(c) = dv(1) + (dnum(c)-datenum(dv(1),1,1))/365;
%     end
% end
