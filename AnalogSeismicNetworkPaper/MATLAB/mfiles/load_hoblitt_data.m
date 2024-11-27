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


%%


'The following channels are transmitted and recorded.'

'BATT	Battery voltage for the field site.  Volts= sub 4096; mul 0.00122; mul 4'
'TEMP	Temperature at the buried sensor. Temp= sub 4096; mul 0.122'
'XAXS	X-axis reading. To get microradians: sub 4096 and div 10'
'YAXS	Y-axis reading. To get microradians: sub 4096 and div 10' }

tx.plot()
ty.plot()
temp.plot()
batt.plot()


%%
close all
rsamdir = 'MVOdata/RSAM_1/';
d=dir(fullfile(rsamdir,'M*.DAT'));
stations = {};
for c=1:numel(d)
    thissta = d(c).name(1:4);
    stations{c} = thissta;
end
stations = unique(stations)
snum=datenum(1995,1,1)
enum=datenum(2004,12,31)
%     dp = 'INETER_DATA/RSAM/SSSSYYYY.DAT';
%     s = rsam.read_bob_file('file', dp, 'snum', datenum(2015,1,1), ...
%           'enum', datenum(2015,2,1), 'sta', 'MOMN', 'units', 'Counts')
rcell={}
stations = stations(1:2)
fout = fopen('RSAM_data_captured.txt','w+');
for stanum=1:numel(stations)
    r=rsam.read_bob_file(fullfile(rsamdir, 'SSSSYYYY.DAT'),'snum',snum,'enum',enum,'sta',stations{stanum});
    rcell{stanum} = r;
    r.plot()
    okdnum=r.dnum(r.data>0);
    dmin = min(okdnum);
    dmax = max(okdnum);
    dlen = numel(okdnum)/1440;
    ddiff = (dmax-dmin);
    dmissing = (ddiff-dlen);
    fprintf(fout,'%s\t%s\t%s\t%10.4f\t%10.4f\n',stations{stanum},datestr(dmin,31),datestr(dmax,31),dlen,dmissing);
end
fclose(fout)
%%

