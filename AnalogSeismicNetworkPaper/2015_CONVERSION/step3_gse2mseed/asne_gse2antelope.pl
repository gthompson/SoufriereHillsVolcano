#!/usr/bin/perl
# This program loops over GSE files in YYYY/MM directories, creates a Seisan-style filenr.lis, and then runs
# the program fix_asne_archive.py on this list.

# GSE files were created by converting DMX files to GSE using the PC-SUDS program SUDS2GSE,
# after time correcting and demultiplexing.
#
# The environment variables GSE_TOP_DIR and VOLCANO_OBSPY need to be set

# Glenn Thompson 2015 University of South Florida
$filenrlis = "gse.lis";
use Cwd;
chdir($GSE_TOP_DIR);
'/raid/data/suds/WVM');
foreach $YYYY (glob("????")) {
	chdir($YYYY);
	print cwd(),"\n";
	foreach $MM (glob("??")) {
		chdir($MM);
		print cwd(),"\n";
		unlink($filenrlis) if (-e $filenrlis);
		@allgsefiles = glob("*.gse");
		printf "%d GSE files\n",($#allgsefiles + 1);
		if ($#allgsefiles > 0) {
			system("python ~/src/volcanoObsPy/seisan_scripts/dirf.py . gse > $filenrlis");
			system("python -W ignore ~/bin/fix_asne_archive.py $filenrlis ASNE_WVM");
			unlink($filenrlis) if (-e $filenrlis);
		}
		chdir("..");
		
	}
	chdir("..");
}
