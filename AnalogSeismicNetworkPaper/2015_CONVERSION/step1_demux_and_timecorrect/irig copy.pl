#!C:\Perl\bin\perl.exe
# This program just applied IRIG.EXE to all DMX files
chdir("d:\\DMX_PHA\\");
@alldirs=glob("????");
foreach $YYYY (@alldirs) {
	chdir("d:\\DMX_PHA\\$YYYY");
	print "\n\nDirectory now e:\\DMX_PHA\\$YYYY\n";
	@alldmxfiles=glob("*.DMX");
	foreach $dmxfile (@alldmxfiles) {
		system("d:\\PROGRAMS\\SUDS\\IRIG.EXE $dmxfile");
	}	
}
