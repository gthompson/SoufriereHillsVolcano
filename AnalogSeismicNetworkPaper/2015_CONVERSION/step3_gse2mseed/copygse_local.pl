#!C:\Perl64\bin\perl.exe
use Cwd;
use File::Copy;
use File::Path;
#chdir("c:\\users\\thompsong\\Desktop\\SUDS_DATA\\DMX_PHA\\");
chdir("d:\\SUDS_DATA\\");
print "Directory: ".cwd()."\n";
@yeardirs=glob("????");
foreach $yeardir (@yeardirs) {
	chdir("$yeardir");
	print "Directory: ".cwd()."\n";
	@monthdirs = glob("??");
	foreach $monthdir (@monthdirs) {
		($targetdir = cwd()) =~ s/d:/m:/;
		chdir($monthdir);
		print "Directory: ".cwd()."\n";
		for ($day=1; $day < 32; $day++) {
			$YYMMDD = sprintf("%s%s%02d",substr($yeardir,2,2), $monthdir, $day);
			@allgsefiles=glob($YYMMDD."*.gse");
			printf "Found %d GSE files\n", ($#allgsefiles + 1); 
			if ($#allgsefiles > -1) {
					unless (-d $targetdir) {
						print "mkpath $targetdir\n";
						mkpath($targetdir);
					}
			}
			print "\nYYMMDD = $YYMMDD\n" if ($#allgsefiles > -1);
			foreach $gsefile (@allgsefiles) {
#print "got here\n";
				#unless (-e "$targetdir/$gse") {
					print "cp $gsefile $targetdir\n";
					copy($gsefile, "$targetdir/");
				#}
			}
		}
		chdir("..");

	}
	chdir("..");	
}
