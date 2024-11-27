#!/usr/bin/perl
use Cwd;
use File::Copy;
#chdir('/home/t/thompsong/data/SUDS_DATA');
chdir('/raid/data/suds');
$obsoletedir = "/raid/data/suds/obsolete/";
foreach $YYYY (glob("????")) {
	chdir($YYYY);
	print cwd(),"\n";
	#foreach $MM (glob("??")) {
	#	chdir($MM);
		foreach $tarfile (glob("*.tar")) {
			system("tar -xf $tarfile");
			unlink($tarfile);
		}
		foreach $wvmfile (glob("*.WVM")) {
			($dmxfile = $wvmfile) =~ s/WVM/dmx/;
			($DMXfile = $wvmfile) =~ s/WVM/DMX/;
                        ($gsefile = $dmxfile) =~ s/WVM/gse/;
			if (-e $dmxfile || -e $DMXfile || -e $gsefile) {
				move($wvmfile, $obsoletedir);
			};
		}
                foreach $dmxfile (glob("*.dmx")) {
                        ($gsefile = $dmxfile) =~ s/dmx/gse/;
			($DMXfile = $dmxfile) =~ s/dmx/DMX/;
                        if (-e $gsefile | -e $DMXfile) {
				move($dmxfile, $obsoletedir);
                        }
                }
                foreach $DMXfile (glob("*.DMX")) {
                        ($gsefile = $DMXfile) =~ s/DMX/gse/;
                        if (-e $gsefile) {
				move($DMXfile, $obsoletedir);
                        }
                }

		foreach $wvmfile (glob("*.WVM")) {
			$gsefile = substr($wvmfile,0,length($wvmfile)-3)."gse";
                        #($gsefile = $dmxfile) =~ s/WVM/gse/;
			if ( -e $gsefile) {
				move($wvmfile, $obsoletedir);
			};
		}
	#	chdir("..");
		
	#}
	chdir("..");
}
