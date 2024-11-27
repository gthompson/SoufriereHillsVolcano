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
				$seislogfile = `head -1 $file | cut -c896-922`;
				$seislogfile =~ s/\W//g;
				$seislogfile =~ s/MVO/\.MVO/;
				print "Extracted filename $seislogfile from $file\n";
				#chomp($seislogfile);	
				if ($seislogfile =~ /T.MVO/) {
					$seislogdir = "/cygdrive/d/roxio/seislog/$yeardir/$monthdir";
					mkpath($seislogdir) unless (-e $seislogdir);
					unless (-e "$seislogdir/$seislogfile") {
						print "Move $file -> $seislogdir/$seislogfile\n";
						move($file, "$seislogdir/$seislogfile");
					}
				} else { # is it a seisan file?
					$timetext = `head -1 $file | cut -c36-59`;
					$numchans = substr($timetext,0,2);
					$year = substr($timetext,3,2);
					$year = "19$year" if (substr($year,0,1) eq "9");
					$year = "20$year" if (substr($year,0,1) eq "0");
					$year = ($year + 8) if ($year < 1995);
					($month = substr($timetext,10,2)) =~ s/\W//g;
					($day = substr($timetext,13,2)) =~ s/\W//g;
					($hour = substr($timetext,16,2)) =~ s/\W//g;
					($minute = substr($timetext,19,2)) =~ s/\W//g;
					($second = substr($timetext,22,2)) =~ s/\W//g;
					$seisanfile = sprintf("$year-%02d-%02d-%02d%02d-%02d",$month,$day,$hour,$minute,$second)."S.MVO___".$numchans;
					print "Extracted filename $seisanfile from $file\n";	
					if ($year > 1994 && $year < 2015) {
						$seisandir = "/cygdrive/d/roxio/seisan/$yeardir/$monthdir";
						mkpath($seisandir) unless (-e $seisandir);
						unless (-e "$seisandir/$seisanfile") {
							print "Move $file -> $seisandir/$seisanfile\n";
							move($file, "$seisandir/$seisanfile");
						}
					} else { # file type unknown
						$unknowndir = "/cygdrive/d/roxio/unknown/$yeardir/$monthdir";
						mkpath($unknowndir) unless (-e $unknowndir);
						unless (-e "$unknowndir/$file") {
							print "Move $file -> $unknowndir/$file\n";
							move($file, "$unknowndir/$file");
						}						
					}
									
				}
				if (-e $file) { # must be a duplicate
					$dupdir = "/cygdrive/d/roxio/duplicates/$yeardir/$monthdir";
					mkpath($dupdir) unless (-e $dupdir);
					unless (-e "$dupdir/$file") {
						print "Move $file -> $dupdir/$file\n";
						move($file, "$dupdir/$file");
					}					
				}

			}
		}		

		chdir("..");

	}
	chdir("..");	
}
