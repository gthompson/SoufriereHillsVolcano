save wavfilesperday2019.mat counts2 dnum2

counts2(dnum>datenum(2000,4,1) & dnum<datenum(2002,1,1)) = 72;
close all
bar(dnum2,counts2/72*100)
datetick('x')
axis tight
set(gca,'YLim',[0 100])
title('DSN continuous data')
ylabel('Percent of daily files captured')
xlabel('Date')
%%
figure
bar(dnum2, smooth(counts2/72*100,30))
datetick('x')
axis tight
set(gca,'YLim',[0 100])
title('DSN continuous data')
ylabel('Percent of monthly files captured')
xlabel('Date')


%%
d = mvoewavfilesperday;
yyyy=d.VarName2;
mm=d.VarName3;
dd=d.VarName4;
counts=d.VarName5;
dnum=datenum(yyyy,mm,dd);
size(counts)
size(dnum)
%%


indices = find(dnum>datenum(1996,10,21) & dnum<datenum(2008,9,1));
counts2 = counts(indices);
dnum2=dnum(indices);

close all
bar(dnum2,counts2)
datetick('x')
axis tight
title('DSN event waveform data')
ylabel('Event WAV files per day')
xlabel('Date')

figure
bar(dnum2, smooth(counts2,30))
datetick('x')
axis tight
title('DSN event waveform data')
ylabel('Event WAV files per month')
xlabel('Date')