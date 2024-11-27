#!/usr/bin/perl
use Cwd;
use File::Copy;
#use File::Basename;
@yeardirs=glob("????");
foreach $yeardir (@yeardirs) {
	chdir($yeardir);
	print "Directory: ".cwd()."\n";
	@monthdirs = glob("??");
	foreach $monthdir (@monthdirs) {
		chdir($monthdir);
		print "Directory: ".cwd()."\n";
		@seislogfiles = glob("*T.*");
		foreach $seislogfile (@seislogfiles) {
			@fileparts = split /T\./, $seislogfile;
			$basename = $fileparts[0];
			#($basename,$path,$suffix) = fileparse($seislogfile);
			#print "$basename\n";
			if (length($basename) == 15) {
				$YY = substr($basename, 0, 2);
				$MM = substr($basename, 2, 2);
				$DD = substr($basename, 5, 2);
				$HH = substr($basename, 8, 2);
				$II = substr($basename, 10, 2);
				$SS = substr($basename, 13, 2);
				if ($YY > 90) {
					$YYYY = "19$YY";
				} else {
					$YYYY = "20$YY";
				}

			}

			if (length($basename) == 18) {
				$YYYY = substr($basename, 0, 4);
				$MM = substr($basename, 5, 2);
				$DD = substr($basename, 8, 2);
				$HH = substr($basename, 11, 2);
				$II = substr($basename, 13, 2);
				$SS = substr($basename, 16, 2);
			}

			$seisanbase = $YYYY."-".$MM."-".$DD."-".$HH.$II."-".$SS."S";
			#print "Looking for $seisanbase\n";
			@seisanfiles = glob($seisanbase."*");
			if ($#seisanfiles == -1) {
				# No Seisan file for this Seislog file
				print "\nqnxsei $seislogfile\n";
				system("echo $seislogfile > qnxsei.in");
				system("/raid/apps/seisan/PRO/qnxsei < qnxsei.in");
			} else {
				mkdir("./done") unless (-d "./done");
				#print "$seislogfile ";
				#move($seislogfile, "./done");
			}
			@seisanfiles = glob($seisanbase."*");
			if ($#seisanfiles == -1) {
				# qnxsei must have crashed
				$movedir = "corrupt";
			}
			else
			{
				$movedir = "done";
			}
			mkdir($movedir) unless (-d $movedir);
                        #print "$seislogfile ";
                        move($seislogfile, $movedir);
                                
		}
		chdir("..");
	}
	chdir("..");
}
