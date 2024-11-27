#!/usr/bin/perl
use Cwd;
use File::Copy;
chdir("/cygdrive/d/roxio/seisan");
@yeardirs=glob("????");
foreach $yeardir (@yeardirs) {
	chdir($yeardir);
	print "Directory: ".cwd()."\n";
	@monthdirs = glob("??");
	foreach $monthdir (@monthdirs) {
		chdir($monthdir);
		print "Directory: ".cwd()."\n";
		@allfiles = glob("*S*");
		foreach $file (@allfiles) {
			($newfile = $file) =~s/ /0/g;
			print "move $file -> $newfile\n" if ($newfile ne $file);
			move($file, $newfile) if ($newfile ne $file);
		}
		chdir("..");
	}
	chdir("..");
}