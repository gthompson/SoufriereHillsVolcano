function histogram_rsam_alarms();
% Author: Glenn Thompson 2001
% Description: Plots a histogram of rsam event and tremor alarms
% Usage: histogram_rsam_alarms(snum,enum);
% binsize = bin size in number of days

load rsam_alarms.mat
close all;
figure;

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
print -dpng -f2 rsam_alarms.png
        
        