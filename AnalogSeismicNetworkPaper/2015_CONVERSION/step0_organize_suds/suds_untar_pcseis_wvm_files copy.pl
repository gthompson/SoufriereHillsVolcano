#!/usr/bin/perl
chdir("/raid/data/suds/WVM");
@yeardirs = <????>;
foreach $YYYY (@yeardirs) {
	chdir("$YYYY");
	@monthdirs = glob("??");
	foreach $MM (@monthdirs) {
		if (-d $MM) {
			chdir($MM);
			foreach $tarfile (glob("*.tar")) {
				system("tar -xvf $tarfile");
			}
			chdir("..");
		}
	}
	chdir("..");
}
