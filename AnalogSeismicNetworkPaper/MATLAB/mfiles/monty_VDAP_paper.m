%% Hypo71 data
summaryfile = '/Volumes/GoogleDrive/My Drive/data/MVO/seismicdata/Hypocentres/summaryall.txt'
cobj = Catalog.retrieve('vdap',summaryfile);

%%
cobj2=cobj;
ind = find(abs(cobj.lat-median(cobj.lat)) < 0.15);
cobj=cobj.subset('indices',ind);
ind = find(abs(cobj.lon-median(cobj.lon)) < 0.15);
cobj=cobj.subset('indices',ind);
%%
erobj=cobj.eventrate('binsize',1)
erobj.plot('metric',{'counts';'energy'})
print -dpng hypo71counts.png
%% RSAM_1 data
close all
rsamdir = 'MVOdata/RSAM_1/';
d=dir(fullfile(rsamdir,'M*.DAT'));
stations = {};
for c=1:numel(d)
    thissta = d(c).name(1:4);
    stations{c} = thissta;
end
stations = unique(stations);
snum=datenum(1995,1,1)
enum=datenum(2004,12,31)
fout = fopen('RSAM_data_captured.csv','w+');
fprintf(fout,'station,ondate,offdate,days,days captured,days missing\n')
for stanum=1:numel(stations)
    r=rsam.read_bob_file(fullfile(rsamdir, 'SSSSYYYY.DAT'),'snum',snum,'enum',enum,'sta',stations{stanum});
    okdnum=r.dnum(r.data>0);
    dmin = min(okdnum);
    dmax = max(okdnum);
    dlen = numel(okdnum)/1440;
    ddiff = (dmax-dmin);
    dmissing = (ddiff-dlen);
    fprintf(fout,'%s,%s,%s,%10.4f,%10.4f,%10.4f\n',stations{stanum},datestr(dmin,31),datestr(dmax,31),ddiff,dlen,dmissing);
    fdnum=floor(okdnum);
    fsta = fopen(sprintf('RSAM_%s.txt',stations{stanum}),'w+');
    fprintf(fsta,'date\t%s\n',stations{stanum});
    for daynum=datenum(1995,7,28):datenum(2004,1,10)
        l=length(find(fdnum==daynum))/1440;
        fprintf(fsta,'%s\t%.2f\n',datestr(daynum,'yyyymmdd'),l);
    end
    fclose(fsta);
end
fclose(fout)

%% TILT data
clc
snum=datenum(1995,7,27);
enum=datenum(2000,10,13);
tiltdir = 'MVOdata/TILT/';
d=dir(tiltdir);
stations = {};
fout = fopen('TILT_data_captured.csv','w+');
fprintf(fout,'station,ondate,offdate,days,days captured,days missing\n');
for c=1:numel(d)
    thissta = d(c).name;
    fprintf('%s %s\n',d(c).folder, d(c).name)
    t = rsam.read_bob_file(fullfile(d(c).folder, thissta, 'XAXSYYYY.DAT'),'snum',snum,'enum',enum);
    okdnum=t.dnum(abs(t.data)>0);
    dmin = min(okdnum);
    dmax = max(okdnum);
    dlen = numel(okdnum)/144;
    ddiff = (dmax-dmin);
    dmissing = (ddiff-dlen);
    %t.plot()
    if dlen<1
        continue
    end
    fprintf(fout,'%s,%s,%s,%10.4f,%10.4f,%10.4f\n',thissta,datestr(dmin,31),datestr(dmax,31),ddiff,dlen,dmissing);
    fdnum=floor(okdnum);
    fsta = fopen(sprintf('TILT_%s.txt',thissta),'w+');
    fprintf(fsta,'date\t%s\n',thissta);
    for daynum=datenum(1995,7,28):datenum(2004,1,10)
        l=length(find(fdnum==daynum))/144;
        fprintf(fsta,'%s\t%.2f\n',datestr(daynum,'yyyymmdd'),l);
    end
    fclose(fsta);
    
end
fclose(fout);

%% RSAM EVENT10
close all, clc
%rsameventsdir = 'MVOdata/RSAM/EVNT_10/';
rsameventsdir = 'MVOdata/VDAP/RSAMTRIGGERS/';
d=dir(fullfile(rsameventsdir,'M*.DAT'));
stations = {};
for c=1:numel(d)
    thissta = d(c).name(1:4);
    stations{c} = thissta;
end
stations = unique(stations);
disp(stations)
snum=datenum(1995,7,27);
enum=datenum(2001,12,31,23,59,59);
fout = fopen('EVENTS_RSAM_data_captured.csv','w+');
fprintf(fout,'station,ondate,offdate,days,days captured,days missing\n');
for stanum=1:numel(stations)
    r=rsam.read_bob_file(fullfile(rsameventsdir, 'SSSSYYYY.DAT'),'snum',snum,'enum',enum,'sta',stations{stanum});
    okdnum=r.dnum(r.data>0);
    dmin = min(okdnum);
    dmax = max(okdnum);
    dlen = numel(okdnum)/144;
    ddiff = (dmax-dmin);
    dmissing = (ddiff-dlen);
    fprintf(fout,'%s,%s,%s,%10.4f,%10.4f,%10.4f\n',stations{stanum},datestr(dmin,31),datestr(dmax,31),ddiff,dlen,dmissing);
    fdnum=floor(okdnum);
    fsta = fopen(sprintf('EVENTS_RSAM_%s.txt',stations{stanum}),'w+');
    fprintf(fsta,'date\t%s\n',stations{stanum});
    for daynum=datenum(1995,7,28):datenum(2004,1,10)
        l=length(find(fdnum==daynum))/144;
        fprintf(fsta,'%s\t%.2f\n',datestr(daynum,'yyyymmdd'),l);
    end
    fclose(fsta);
    if dlen>365
        r.plot()
    end
end
fclose(fout);
%%
close all
r1=rsam.read_bob_file(fullfile(rsameventsdir, 'SSSSYYYY.DAT'),'snum',snum,'enum',enum,'sta','MWHT');
r2=rsam.read_bob_file(fullfile(rsameventsdir, 'SSSSYYYY.DAT'),'snum',snum,'enum',enum,'sta','MLGT');
r3=rsam.read_bob_file(fullfile(rsameventsdir, 'SSSSYYYY.DAT'),'snum',snum,'enum',enum,'sta','MGAT');
r=r1;
y=[r1.data; r2.data; r3.data];
size(y)
datetick('x')
set(gca,'XLim',[snum enum])
ylabel('RSAM Events per 10 minutes (sum of MWHT, MLGT, MGAT)')
r.data = nansum(y,1);
size(r.data)
size(r.dnum)
figure
plot(r.dnum, r.data)
%%
% Author: Glenn Thompson 2001
% Description: Plots a histogram of rsam event and tremor alarms
% Usage: histogram_rsam_alarms(snum,enum);
% binsize = bin size in number of days

load rsam_alarms.mat

% event alarms
% i=find(dnum_e>snum & dnum_e<enum);
% dnum_e=dnum_e(i);
%y=ones(length(dnum_e),1);
%[t,y]=bindata(dnum_e,y,binsize);
%bar(t,y);
counts_e=1:length(dnum_e);
snum=floor(min(dnum_e));
enum=ceil(max(dnum_e));
figure
stairs([snum dnum_e],[0 counts_e],'k-','LineWidth',2);
datetick('x')
hold on
set(gca,'XLim',[snum enum]);
grid;
%title(['RSAM event alarms']);
ylabel('cumulative # of alarms');

        
%         fout1=fopen('C:\earthworm\html\alarms\latest_RSAM_EVENT_alarms.html','w');
%         fprintf(fout1,'<html><title>Latest RSAM EVENT alarms</title>\n');
%         fprintf(fout1,'<table border=1>\n<tr><td>Time</td></tr>');
%         for alarmnum=1:length(dnum_e)
%             fprintf(fout1,'<tr><td>%s</td></tr>\n',datestr(dnum_e(alarmnum),31));
%         end
%         fprintf(fout1,'</table></html>\n');
%         fclose(fout1);

% tremor alarms
j=find(dnum_t>snum & dnum_t<enum);
dnum_t=dnum_t(j);
%y=ones(length(dnum_t),1);
%[t,y]=bindata(dnum_t,y,binsize);
%bar(t,y);
counts_t=1:length(dnum_t);
%subplot(2,1,2),
stairs([snum dnum_t],[0 counts_t],'r-','LineWidth',2);
datetick('x')
set(gca,'XLim',[snum enum]);
grid on;
%title(['RSAM tremor alarms']);
%xlabel(sprintf('Last updated at %s',datestr(now,31)));
ylabel('cumulative # of alarms');
legend({'event';'tremor'})

%         fout1=fopen('C:\earthworm\html\alarms\latest_RSAM_TREMOR_alarms.html','w');
%         fprintf(fout1,'<html><title>Latest RSAM TREMOR alarms</title>\n');
%         fprintf(fout1,'<table border=1>\n<tr><td>Time</td></tr>');
%         for alarmnum=1:length(dnum_t)
%             fprintf(fout1,'<tr><td>%s</td></tr>\n',datestr(dnum_t(alarmnum),31));
%         end
%         fprintf(fout1,'</table></html>\n');
%         fclose(fout1);
%%
print -dpng  rsam_alarms.png

%% SP event counts
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

%%
close all, clc
hoblittdir = '~/Desktop/data_mastering_paper/Hoblitt_MRAT_DATA';
snum = datenum(1995,7,1);
enum = datenum(1995,12,31,23,59,59);

tx = rsam.read_bob_file(fullfile(hoblittdir, 'XAXS95.DAT'),'snum',snum,'enum',enum);
ty = rsam.read_bob_file(fullfile(hoblittdir, 'YAXS95.DAT'),'snum',snum,'enum',enum);
temp = rsam.read_bob_file(fullfile(hoblittdir, 'TEMP95.DAT'),'snum',snum,'enum',enum);
batt = rsam.read_bob_file(fullfile(hoblittdir, 'BATT95.DAT'),'snum',snum,'enum',enum);
tx.data = (tx.data - 4096) / 10; % microradians
tx.units = 'microradians';
ty.data = (ty.data - 4096) / 10; % microradians
ty.units = 'microradians';
temp.data = (temp.data - 4096) * 0.122; % Degrees - but is it F or C?
temp.units = 'degrees';
batt.data = (batt.data - 4096) * 0.00122 * 4; % Volts
batt.units = 'Volts';

location = Position(16 + 43.51/60, -(62 + 9.74/60));

info = {'LONG GROUND TILTMETER SITE -  Station installed 8/03/95'

'Tiltmeter is an Applied Geomechanics 701-1 bi-axial tiltmeter'
'installed approximately 1.57 meters below the surface.  The tiltmeter is about'
'7 feet east of the Long Ground Health Center front porch and 5 feet south of '
'the sidewalk.  Location is 16 deg. 43.51 min. North and 62 deg. 09.74 min. West.'
'Station was installed with Filter=out and Gain=high.'

'The following channels are transmitted and recorded.'

'BATT	Battery voltage for the field site.  Volts= sub 4096; mul 0.00122; mul 4'
'TEMP	Temperature at the buried sensor. Temp= sub 4096; mul 0.122'
'XAXS	X-axis reading. To get microradians: sub 4096 and div 10'
'YAXS	Y-axis reading. To get microradians: sub 4096 and div 10' }

tx.plot()
ty.plot()
temp.plot()
batt.plot()
