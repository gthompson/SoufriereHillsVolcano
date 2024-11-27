#!C:\Perl64\bin\perl.exe
use Cwd;
chdir("y:\\suds\\DMX\\");
@yeardirs=glob("????");
foreach $yeardir (@yeardirs) {
	chdir("$yeardir");
	@monthdirs = glob("??");
	foreach $monthdir (@monthdirs) {
		chdir($monthdir);
		print "Directory: ".cwd()."\n";
		@alldmxfiles=glob("*.DMX");
		foreach $dmxfile (@alldmxfiles) {
#			$MM=substr($dmxfile,2,2);
			print "Processing $dmxfile\n";
			system("d:\\PROGRAMS\\sudsei_mvo.exe $dmxfile");
#        	        if (!(-e "d:\\SEISAN\\WAV\\IRIG_\\$YYYY")) {
#				mkdir("d:\\SEISAN\\WAV\\IRIG_\\$YYYY",777);
#			}
#        	        if (!(-e "e:\\SEISAN\\WAV\\IRIG_\\$YYYY\\$MM")) {
#				mkdir("d:\\SEISAN\\WAV\\IRIG_\\$YYYY\\$MM",777);
#			}
#        	        #print "\n\nMoving *.IRIG* to d:\\SEISAN\\WAV\\IRIG_\\$YYYY\\$MM\n\n";
#        	        system("move *.IRIG* d:\\SEISAN\\WAV\\IRIG_\\$YYYY\\$MM");
	}	
}
