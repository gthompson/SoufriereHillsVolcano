#!C:\Perl\bin\perl.exe
# This program applies IRIG.EXE to all DMX files
# IRIG.EXE is a PC_SUDS utility that uses the IRIG time trace in a DMX file to correct the time of the SUDS file
# This was needed because computer clocks drifted. The IRIG time trace contained encoded the seconds, minutes, hours and days from a 1-s GPS time signal
#
# The environment variables SUDS_TOP_DIR and SUDS_BIN need to be set
# DMX files are demultiplexed SUDS files stored in year directories under SUDS_TOP_DIR
# PC-SUDS utilities are executables stored under SUDS_BIN

# Glenn Thompson 2015 University of South Florida

use Cwd;
use Env;
$SUDS_TOP_DIR = $ENV{SUDS_TOP_DIR};
$SUDS_BIN = $ENV{SUDS_BIN};
die("$SUDS_TOP_DIR does not exist\n") unless (-e $SUDS_TOP_DIR);
die("$SUDS_BIN does not exist\n") unless (-e $SUDS_BIN);
chdir($SUDS_TOP_DIR);
@yeardirs=glob("????");
foreach $YYYY (@yeardirs) {
	chdir("$YYYY");
	$pwd = getcwd();
	print "\nDirectory now $pwd\n";
	@alldmxfiles=glob("*.DMX");
	foreach $dmxfile (@alldmxfiles) {
		system("$SUDS_BIN\\IRIG.EXE $dmxfile");
	}	
}
