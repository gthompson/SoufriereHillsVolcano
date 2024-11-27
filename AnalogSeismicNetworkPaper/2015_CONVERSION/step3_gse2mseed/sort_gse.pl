#!/usr/bin/perl
use Cwd;
use File::Copy;
use File::Path;
#chdir('/home/t/thompsong/data/SUDS_DATA');
chdir('/raid/data/suds/DMX');
foreach $YYYY (glob("????")) {
	chdir($YYYY);
	print cwd(),"\n";
	foreach $MM (glob("??")) {
		chdir($MM);
		$targetdir='/raid/data/suds/GSE/$YYYY/$MM';
		mkpath($targetdir);

                foreach $gsefile (glob("*.gse")) {
			#if (-e "/home/t/thompsong/data/gse/$gsefile") {
			#	unlink($gsefile);
			#} else {
				copy($gsefile,  $targetdir);
			#}
		}

		chdir("..");
		
	}
	chdir("..");
}
