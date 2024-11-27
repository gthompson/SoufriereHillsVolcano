#!C:\Perl64\bin\perl.exe
use Cwd;
use File::Copy;
print "Directory: ".cwd()."\n";
@yeardirs=glob("????");
foreach $yeardir (@yeardirs) {
	chdir("$yeardir");
	print "Directory: ".cwd()."\n";
	@monthdirs = glob("??");
	foreach $month (@monthdirs) {
		chdir($month);
		#@months = qw(01 02 03 04 05 06 07 08 09 10 11 12);
		#foreach $month (@months) {
		#print "$month\n";
		$YYMM = sprintf("%s%s",substr($yeardir,2,2), $month);
		@allfiles = glob($YYMM."*");
		printf("There are %d files matching %s\n", ($#allfiles + 1), $YYMM);
		if ($#allfiles > -1) {
			mkdir($month);
			#print "Directory: ".cwd()."\n";
			for ($day=1; $day < 32; $day++) {
				$YYMMDD = sprintf("%s%s%02d",substr($yeardir,2,2), $month, $day);
				@allwvmfiles=glob($YYMMDD."*.WVM");
				printf("\n\WVM: $YYMMDD %d files\n",($#allwvmfiles+1)) if ($#allwvmfiles > -1);
				foreach $wvmfile (@allwvmfiles) {
					($dmx = $wvmfile) =~ s/WVM/dmx/;
					($gse = $dmxfile) =~ s/dmx/gse/;
					unless (-e $dmx || -e $gse) {
						print "Processing $wvmfile\n";
						eval {
							system("demux.exe $wvmfile");
						}
					}
					unlink($wvmfile) if (-e $dmx);
				}
	
				@alldmxfiles=glob($YYMMDD."*.dmx");
				printf("\ndmx $YYMMDD: %d files\n",(@alldmxfiles+1)) if ($#alldmxfiles > -1);
				foreach $dmxfile (@alldmxfiles) {
					($gse = $dmxfile) =~ s/dmx/gse/;
					unless (-e $gse) {
						print "Processing $dmxfile\n";
						eval {
							system("sud2gse.exe $dmxfile");
						}
					}
					unlink($dmxfile) if (-e $gse);
				}
				
				@allDMXfiles=glob($YYMMDD."*.DMX");
				print "\nDMX $YYMMDD:\n" if ($#allDMXfiles > -1);
				foreach $DMXfile (@allDMXfiles) {
					($gse = $DMXfile) =~ s/DMX/gse/;
					unless (-e $gse) {
						print "Processing $DMXfile\n";
						eval {
							system("sud2gse.exe $DMXfile");
						}
					}
					unlink($DMXfile) if (-e $gse);
				}
			}

			#@allgsefiles = glob($YYMM."*.gse");
			#foreach $gsefile (@allgsefiles) {
			#	move($gsefile,"$month\\$gsefile");
			#}
		}
		chdir("..");

	}
	chdir("..");	
}
