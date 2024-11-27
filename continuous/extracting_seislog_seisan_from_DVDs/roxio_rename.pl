#!/usr/bin/perl
use Cwd;
use File::Copy;
use File::Path;
print "Directory: ".cwd()."\n";
@yeardirs=glob("????");
foreach $yeardir (@yeardirs) {
	chdir("$yeardir");
	print "Directory: ".cwd()."\n";
	@monthdirs = glob("??");
	foreach $monthdir (@monthdirs) {
		chdir($monthdir);
		print "Directory: ".cwd()."\n";

		# unzip
		@allfiles = glob("*.BZ2");
		printf(cwd().": There are %d BZ2 files\n", ($#allfiles + 1));
		if ($#allfiles > -1) {
			foreach $file (@allfiles) {
				system("bunzip2 $file");
			}
		}

		# rename
		@allfiles = glob("*.out");
		printf(cwd().": There are %d OUT files\n", ($#allfiles + 1));
		if ($#allfiles > -1) {
			foreach $file (@allfiles) {
				$seislogfile = `head -1 $file | cut -c897-921`;
				if ($seislogfile =~ /T.MVO/) {
					$seislogdir = "/cygdrive/d/roxio/seislog/$yeardir/$monthdir";
					mkpath($seislogdir) unless (-e $seislogdir);
					move($file, "$seislogdir/$seislogfile");
				}
			}
		}		

		chdir("..");

	}
	chdir("..");	
}
