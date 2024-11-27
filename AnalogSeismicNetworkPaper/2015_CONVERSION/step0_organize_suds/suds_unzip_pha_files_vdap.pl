#!/usr/bin/perl
@zips = glob("/raid/data/suds/DMX_PHA/DONE_PHA/*.ZIP");
foreach $zip (@zips) {
	system("unzip $zip");
}
