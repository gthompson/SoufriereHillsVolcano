#!C:\Perl64\bin\perl.exe
use Cwd;
use File::Copy;
chdir("y:\\suds\\WVM\\");
@yeardirs=glob("????");
foreach $yeardir (@yeardirs) {
	chdir($yeardir);
	print "Directory: ".cwd()."\n";
	@monthdirs = glob("??");
	foreach $monthdir (@monthdirs) {
		chdir($monthdir);
		print "Directory: ".cwd()."\n";
		for ($day=1; $day < 32; $day++) {
			$YYMMDD = sprintf("%s%s%02d",substr($yeardir,2,2), $monthdir, $day);
			@allwvmfiles=glob($YYMMDD."*.WVM");
			print "\nYYMMDD = $YYMMDD\n" if ($#allwvmfiles > -1);
			foreach $wvmfile (@allwvmfiles) {
				($DMX = $wvmfile) =~ s/WVM/DMX/;
				($DMXDIR = cwd()) =~ s/WVM/DMX/;
				($dmx = $wvmfile) =~ s/WVM/dmx/;
				unless ((-e $dmx) || (-e "$DMXDIR/$DMX") ) {
					# demux WVM -> dmx unless dmx or DMX exist
					print "demux $wvmfile\n";
					system("demux.exe $wvmfile");
				}
				if (-e $dmx) {
					# if dmx exists and DMX doesn't, move dmx -> DMX
					unless (-e "$DMXDIR/$DMX") {
						print "move $dmx -> $DMXDIR/$dmx\n";
						move($dmx, "$DMXDIR/$DMX");
					} else {
						# DMX exists, so just rename dmx
						print "move $dmx -> done$dmx\n";
						move($dmx, "done".$dmx);
					}
				}
				# we've processed this wvm file
				print "move $wvmfile -> done$wvmfile\n";
				move($wvmfile, "done".$wvmfile);
			}
		}
		chdir("..");
	}
	chdir("..");
}
