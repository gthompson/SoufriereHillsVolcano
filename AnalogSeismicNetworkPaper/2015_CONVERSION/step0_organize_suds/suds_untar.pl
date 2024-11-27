#!/usr/bin/perl
use Cwd;
print "Directory: ".cwd()."\n";
@yeardirs=glob("????");
foreach $yeardir (@yeardirs) {
	chdir("$yeardir");
	print "Directory: ".cwd()."\n";
	@monthdirs = glob("??");
	foreach $monthdir (@monthdirs) {
		chdir($monthdir);
		print "Directory: ".cwd()."\n";
		#system("tar -xvf *.tar");
		@alltarfiles = glob("*.tar");
		foreach $tarfile (@alltarfiles) {
			system("tar -xf $tarfile");
		}
		chdir("..");

	}
	chdir("..");	
}
