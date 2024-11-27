function rsam_alarm2excel();
% Author: Glenn Thompson 2001
% Usage: rsam_alarm2excel;
% Description: Converts rsam_alarms.mat file to Excel format
% Output file is chosen by user

% load the input file
load rsam_alarms.mat

% user select output file
[newfile,newpath]=uiputfile(['*.xls'],'Export data to Excel file');
newfile=[newfile,'.xls'];
fullpath=[newpath,'/',newfile];
fout=fopen(fullpath,'w');

% output event alarm data in Excel format
l=length(data_e);
fprintf(fout,'EVENT ALARMS\t');
for cc=1:length(stations)
   fprintf(fout,'%s\t\t',stations{cc});
end
fprintf(fout,'\n\t');
for cc=1:length(stations)
   fprintf(fout,'THRESH\tDATA\t');
end
for c=1:l
   fprintf(fout,'\n%s',datestr(dnum_e(c)));
   for cc=1:length(stations)
      fprintf(fout,'\t%4.0f',threshold_e(c,cc));
      fprintf(fout,'\t%4.0f',data_e(c,cc));
   end
end
fclose(fout);
