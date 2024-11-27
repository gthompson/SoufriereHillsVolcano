README.TXT 12-Jan-1995 RB

This file documents three utility programs for making line data files 
for HypoMap.  These programs are; LIN2MAP, MAP2LIN, and VPFGet.  LIN2MAP 
and MAP2LIN simply convert ASCII format line data files (.LIN files) to 
binary (.MAP files) and binary files to ASCII, respectively.  These 
programs can be useful when editing or converting line data from other 
sources.

VPFGet is a small program that I threw together to extract line data 
from the DCW (The Digital Chart of the World) package and create binary 
line data files for HypoMap.  This program is a quick hack and may have 
some problems.  Let me know if you have any trouble or just need help 
getting it to work.

The idea is to first run the VPFView software that comes with DCW and 
select your study area in location and extent.  Once this is done you 
should generate a "feature location list" on disk for the features that 
you want.  This file is then processed by VPFGet to extract these 
features from the DCW and write a binary line data file for HypoMap.  
This process is best illustrated by example.

Lets say we want to make a map of the western US.  We want coastlines, 
international borders and administrative (state) borders on the map in 
various colors.  The following steps will do the job:

1) Start VPFView and choose the NOAMER database.

2) Select the center of your study area somewhere near the center of the 
   western US.

3) Choose the Zoom option from the Graphics menu.  Note that you don't
   need any features to be selected or displayed in order to do this.

4) Choose to Zoom Out (the default study area is 5 degrees across, we 
   want more).  You should now see the map of North America again.

5) Press the right mouse button until the area is large enough and 
   center over the area that we want.  Press the left mouse button to
   select it.

6) Now have VPFView generate feature location lists for each feature 
   that we want.  Select Feature location list from the List menu, select
   Disk when prompted for the destination, and then select
   Political/Oceans, Coastlines and then wait for the filename prompt.
   Type COAST.FLL at the filename prompt.

7) Repeat 6 for each desired feature.  Select Political/Oceans, 
   International Boundaries.  Name this file NATION.FLL.  Select
   Political/Oceans, Administrative Boundaries.  Name this file ADMIN.FLL.

8) Quit VPFView.

Now we'll use VPFGet to retrieve the line data from the database using 
the feature location lists that we just generated.  VPFGet expects 4 
arguments on the command line as follows (type VPFGet ? to see this help):
   
   Usage: VPFGET FLLFile MAPFile Color LineType

   FFLFile  = Feature location list from VPFVIEW
   MAPFile  = Output map linedata file for HYPOMAP
   Color    = Color of these lines (0-255)
   LineType = Line type of these lines (0-7)

Here's a list of the first 16 available colors that are available:

0 = Black, 1 = Blue, 2 = Green, 3 = Cyan
4 = Red, 5 = Magenta, 6 = Brown, 7 = White
8 = Gray, 9 = Light Blue, 10 = Light Green, 11 = Light Cyan
12 = Light Red, 13 = Light Magenta, 14 = Yellow, 15 = Bright White

Linetypes are: 0 = solid, 1-6 = dash-dot lines, 7 = dotted.

Now, the feature location list (.FLL) contains pathnames relative to the
main feature subdirectory on the CD-ROM.  So in the following 
discussion, I'll assume that your CD-ROM drive is drive S:.  If this is 
not the case, substitute you drive letter for S: below.  The VPFGet 
program must be in a directory that is on your path at this point or you 
will have to use the fully qualified filespec for VPFGet in the commands 
below.  I will also assume that VPFView is installed in C:\VPFVIEW and 
that the map files that we will create are to go into C:\MAPS.

1) Log into drive S: (type S: and press return).

2) Change directory to S:\DCW\PO (type CD \DCW\PO and press return). This
   is the Political/Ocean feature directory on the CD-ROM.

3) Extract the coastlines by typing:

     VPFGET C:\VPFVIEW\COAST.FLL C:\MAPS\COAST.MAP 7 0

   This will create the map file, C:\MAPS\COAST.MAP which contains the
   coastlines for the western US as white, solid lines.

4) Repeat 3 for NATION.FLL and ADMIN.FLL setting the color and linetype 
   as you wish.

5) Now log back into the C: drive and change directory to C:\MAPS.

6) Concatenate these map files into one file with the following command:

      COPY /B COAST.MAP+NATION.MAP+ADMIN.MAP WEST_US.MAP

   Note that the /B switch in the command above is very important!

Now you have the western US in a map file for use with HypoMap.  Edit your
HYPOMAP.INI and set: 

   Line_Data_File = C:\MAPS\WEST_US.MAP
   Origin_Latitude = 37.0
   Origin_Longitude = -120.0
   Extent = 20.0.

Now run HypoMap and you should see the western United States.

Robert Banfill
Small Systems Support
2 Boston Harbor Place
Big Water,  UT  84741-2205  USA
(801) 675-5827 Voice, (801) 675-5820 Fax
banfill@gldvxa.cr.usgs.gov
