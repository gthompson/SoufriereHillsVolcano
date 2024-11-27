#!/usr/bin/perl
$filenrlis = "gse.lis";
use Cwd;
chdir('/raid/data/suds/WVM');
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
