#!/usr/bin/perl
open($FH, "<", "badDMXfileList.txt") or die("Could not open file\n");
while ($file = <$FH>) {
	chomp($file);
	$file =~ s/\\/\//g;

	$fileroot = substr($file,0,-5);

	#print("<ENTER> to continue\n");
	#$key = <STDIN>;
	unless (-e "/c/Seismo/WAV/ASNE_/$fileroot.mseed") {
	    if (-e "/c/Seismo/WAV/ASNE_/$fileroot.gse") {
		    print("Need to process this GSE file $fileroot\n");
		    system("echo $fileroot >> GSEfilesToProcess2.txt");
	    } else {

	        print("\n$file:\n");
	        system("ls -l /z/Montserrat/MASTERING/VDAP/SUDS/$fileroot".".*");
	        system("ls -l /z/Montserrat/MASTERING/VDAP/GSE/gse_all/$fileroot".".*");
	        system("ls -l $fileroot".".*");
	        system("ls -l /c/Seismo/WAV/ASNE_/$fileroot".".*");
	    };
	};
		
}
close($FH);
# Use this to check for bad files after running this:
#     grep -i "^-.*DMX$" actuallyBadDMXfiles.txt | grep -v "/z/" 
# Then also examine GSEfilestoProcess.txt -> these can be directly added to the ASNE_ database.
# Seisan can handle GSE!
