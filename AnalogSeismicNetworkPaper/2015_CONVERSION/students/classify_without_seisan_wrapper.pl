#!/usr/bin/perl
#use Term::Menus;
#use Switch;
#use Time::Piece;
#use Time::Seconds;
##use DateTime;
our ($year,$month,$day,$hour,$minute,$sec);
$year = 1996;
$month = 4;
$day = 1;
$hour = 0;
$minute = 0;
$sec = 0;
while (1) {

	system("clear");
	print("*** enter all times as Montserrat local time, as they are listed in the books ***\n");
	print("*** this program will convert to UTC ***\n");
	printf("Year [%4d]: ?",$year);
	$newyear = <STDIN>;
	chomp($newyear);
	$newyear = $year if ($newyear eq "");
	if ($newyear >= 1995 and $newyear <= 2008) 
	{
		$year = $newyear;
	} else {
		print "Invalid year\n";
		next;
	}

	printf("Month [%02d]: ?",$month);
	$newmonth = <STDIN>;
	chomp($newmonth);
	$newmonth = $month if ($newmonth eq "");
	if ($newmonth > 0  and $newmonth <= 12) 
	{
		$month = $newmonth;
	} else {
		print "Invalid month\n";
		next;
	}

	printf("Day [%02d]: ?",$day);
	$newday = <STDIN>;
	chomp($newday);
	$newday = $day if ($newday eq "");
	if ($newday > 0  and $newday <= 31) 
	{
		$day = $newday;
	} else {
		print "Invalid day\n";
		next;
	}

	while (1) {

	printf("Hour [%02d]: ?",$hour);
	$newhour = <STDIN>;
	chomp($newhour);
	$newhour = $hour if ($newhour eq "");
	if ($newhour >= 0  and $newhour <= 23) 
	{
		$hour = $newhour;
	} else {
		print "Invalid hour\n";
		next;
	}

	printf("Minute [%02d]: ?",$minute);
	$newminute = <STDIN>;
	chomp($newminute);
	$newminute = $minute if ($newminute eq "");
	if ($newminute >= 0  and $newminute <= 59) 
	{
		$minute = $newminute;
	} else {
		print "Invalid minute\n";
		next;
	}

	# now add 4 hours to time from local to UTC
	$dtyear = $year; $dtmonth  = $month; $dtday = $day; $dthour = $hour; $dtminute = $minute;
	$dthour = $dthour + 4;
	@dayspermonth = qw(31 28 31 30 31 30 31 31 30 31 30 31);
	$dayspermonth[1] = 29 if (($dtyear % 4) == 0);	
	if ($dthour>23) {
		$dtday++;
		$dthour-=24;
		if ($dtday>$dayspermonth[$dtmonth-1]) {
			$dtmonth++;
			$dtday=1;
			if ($dtmonth>12) {
				$dtyear++;
				$dtmonth=1;
			}
		}
	}

	#my $dt = DateTime->new( year => $year, month => $month, day => $day, hour => $hour, minute => $minute, second => 0);
	#$dt->add( hours => 4);
	printf("\n *** Converting to UTC (adding 4 hours) ***\n");

	#printf("Date: ".$dt->year."/%02d/%02d Time: %02d:%02d ",$dt->month,$dt->day,$dt->hour,$dt->minute);
	printf("Date: ".$dtyear."/%02d/%02d Time: %02d:%02d ",$dtmonth,$dtday,$dthour,$dtminute);
	printf("\n\nClassification\n");

	print "1. Rockfall (LV r)\n" .
	      "2. Long-period rockfall (LV e)\n" .
	      "3. Long-period (LV l)\n" .
	      "4. Hybrid (LV h)\n" .
	      "5. Volcano-tectonic (LV t)\n" .
	      "6. Local (L)\n" .
	      "7. Regional (R)\n" .
	      "8. Distal/Teleseism (D)\n" .
	      "9. unknown (u)\n" .
	      "10. noise (n)\n" .
	      "11. OTHER \n";

	my $input = <STDIN>;
	chomp($input);
	
	if($input eq '1') { $class = "LV"; $subclass = "r"; }
	elsif($input eq '2') { $class = "LV"; $subclass = "e"; }
	elsif($input eq '3') { $class = "LV"; $subclass = "l"; }
	elsif($input eq '4') { $class = "LV"; $subclass = "h"; }
	elsif($input eq '5') { $class = "LV"; $subclass = "t"; }
	elsif($input eq '6') { $class = "L"; $subclass = ""; }
	elsif($input eq '7') { $class = "R"; $subclass = ""; }
	elsif($input eq '8') { $class = "D"; $subclass = ""; }
	elsif($input eq '9') { $class = "L"; $subclass = "u"; }
	elsif($input eq '10') { $class = "L"; $subclass = "n"; }
	elsif($input eq '11') { 
		print("Class (e.g. LV) ? ");
		$class = <STDIN>;
		chomp($class); 
		print("subclass (leave blank if none) ? ");
		$subclass = <STDIN>;
		chomp($subclass);
	}

	

	#printf("\nDate: ".$dt->year."/%02d/%02d Time: %02d:%02d ",$dt->month,$dt->day,$dt->hour,$dt->minute);
	printf("\nDate: ".$dtyear."/%02d/%02d Time: %02d:%02d ",$dtmonth,$dtday,$dthour,$dtminute);
	printf("Class: $class Subclass: $subclass \n\n");

	if(prompt_yn("Create S-file ?")) {
		#print("~/bin/event_without_waveform.pl ".$dt->year." ".$dt->month." ".$dt->day." ".$dt->hour." ".$dt->minute." 0 $class $subclass\n");
		print("classify_without_seisan.pl ".$dtyear." ".$dtmonth." ".$dtday." ".$dthour." ".$dtminute." 0 $class $subclass\n");
		#system("~/bin/event_without_waveform.pl ".$dt->year." ".$dt->month." ".$dt->day." ".$dt->hour." ".$dt->minute." 0 $class $subclass");
		system("classify_without_seisan.pl ".$dtyear." ".$dtmonth." ".$dtday." ".$dthour." ".$dtminute." 0 $class $subclass");
	}	

	unless(prompt_yn("Another event for this day ?")) {
		last;
	}

	}


	unless(prompt_yn("Another day ?")) {
		die("Quitting\n");
	}
	
	printf("\n\n\n");
}

sub prompt {
  my ($query) = @_; # take a prompt string as argument
  local $| = 1; # activate autoflush to immediately show the prompt
  print $query;
  chomp(my $answer = <STDIN>);
  return $answer;
}

sub prompt_yn {
  my ($query) = @_;
  my $answer = prompt("$query (Y/N): ");
  return lc($answer) eq 'y';
}
