function text_rsam_alarms(snum,enum);
% Author: Glenn Thompson
% Description: Writes a list of RSAM alarms between given dates/times
% Usage: text_rsam_alarms(snum,enum);

load rsam_alarms.mat

% user select output file
[newfile,newpath]=uiputfile(['*.dat'],'Choose text file');
newfile=[newfile,'.dat'];
fullpath=[newpath,'/',newfile];
fout=fopen(fullpath,'w');

% Event alarms
fprintf(fout,'\nEvent alarms between %s and %s\n',datestr(snum),datestr(enum));
dnumlast=0;
i=find(dnum_e > snum);
dnum_temp1=dnum_e(i);
i=find(dnum_temp1 < enum);
dnum_temp2=dnum_temp1(i);

dnum_last = 0;
count = 0;
for c=1:length(dnum_temp2)
   if dnum_temp2(c) > dnum_last+10/1440
      count = count + 1;
      fprintf(fout,'%d:\t%s\n',count,datestr(dnum_temp2(c)));
      dnum_last=dnum_temp2(c);
   end
end
disp(' ');

% Tremor alarms
fprintf(fout,'\nTremor alarms between %s and %s\n',datestr(snum),datestr(enum));
dnumlast=0;
i=find(dnum_t > snum);
dnum_temp1=dnum_t(i);
i=find(dnum_temp1 < enum);
dnum_temp2=dnum_temp1(i);

dnum_last = 0;
count = 0;
for c=1:length(dnum_temp2)
   if dnum_temp2(c) > dnum_last+9/1440
      count = count + 1;
      fprintf(fout,'%d:\t%s\n',count,datestr(dnum_temp2(c)));
      dnum_last=dnum_temp2(c);
   end
end

fclose(fout);
