#!C:\Perl\bin\perl.exe
chdir("d:\\DMX_PHA\\");
@alldirs=glob("????");
print "Found directories @alldirs\n";
open(FOUT,">D:\\ELIANA\\no_wav_file_found.txt");
foreach $dir (@alldirs) {
	chdir("d:\\DMX_PHA\\$dir");
	print "\n\nDirectory now d:\\WVM\\$dir\n";
	@allphafiles=glob("*.PHA");
	print "Found $#allphafiles for $dir\n";
	foreach $phafile (@allphafiles) {
		#print "Processing $phafile - ";
		$dmxfile = $phafile;
		$dmxfile =~ s/PHA/DMX/;
		if (-e $dmxfile) {
			print "copying $dmxfile to D:\\ELIANA";
                	system("copy $dmxfile D:\\ELIANA");
                	system("copy $phafile D:\\ELIANA");
                } else {
                	print "$dmxfile not found\n";
                	print OUT "$phafile: $dmxfile not found\n";
                }
	}
}
close(FOUT);
