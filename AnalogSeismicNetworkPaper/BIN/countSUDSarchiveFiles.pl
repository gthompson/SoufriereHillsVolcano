#!/usr/bin/perl
foreach my $EXT ("PUN", "PHA", "WVM", "gse", "mseed", "DMX") {
	my $count1 = 0;
	my $count2 = 0;
	my $countu = 0;
	my @yyyydirs = glob("????");
	my $startdir = "/c/Users/thompsong/SUDS";
	chdir($startdir);
	foreach my $yyyydir (@yyyydirs) {
		my @mmdirs = glob("$yyyydir/??");
		foreach my $mmdir1 (@mmdirs) {
			chdir($mmdir1);
			my @files1;
			my $files2;
			my $filesu;
			if ($EXT eq "DMX") {
				my $ext = lc $EXT;
				my @files1u = glob("*.$EXT");
				my @files1l = glob("*.$ext");
				@files1 = find_unique(@files1u, @files1l);
			} else {
				@files1 = glob("*.$EXT");
			}
	
			my $mmdir2 = "/c/seismo/WAV/ASNE_/$mmdir1";
			chdir($mmdir2);
			if ($EXT eq "DMX") {
				my $ext = lc $EXT;
				my @files2u = glob("*.$EXT");
				my @files2l = glob("*.$ext");
				@files2 = find_unique(@files2u, @files2l);
			} else {
				@files2 = glob("*.$EXT");
			}
	
			my @filesu = find_unique(@files1, @files2);
	
			#print("$mmdir1\t$#files1\t$#files2\t$#filesu\n");
			$count1 = $count1 + 1 + $#files1;
			$count2 = $count2 + 1 + $#files2;
			$countu = $countu + 1 + $#filesu;
			chdir($startdir);
		}
	}
	print("TOTAL $EXT: $count1\t$count2\t$countu\n");
}
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

sub find_unique {
	my @all = @_;
	my @roots = ();
	foreach my $file (@all) {
		my $fileroot = substr($file,0,8);
		push @roots, $fileroot;
	}
	my @uniqroots = uniq(@roots);
	return @uniqroots;
}
