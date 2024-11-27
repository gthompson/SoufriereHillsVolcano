#!/usr/bin/perl 
#use Time::gmtime();
use File::Path qw(make_path);
if ($#ARGV<6) {
print("Usage:\n$0 YYYY MM DD hour minute seconds mainclass subclass\n");
print("e.g:\n$0 1996 04 01 23 34 22.0 LV r\n");
$argsgot = $#ARGV + 1;
die("Expect 7 or 8 arguments, but got $argsgot\n");
}
use Env;
$SEISAN_TOP=$ENV{SEISAN_TOP};
my ($yyyy,$mm,$dd,$hh,$mi,$ss,$class,$subclass) = @ARGV;
#my $sfiledir = sprintf("/Users/tom/Desktop/seismo/REA/ASNE_/%4d/%02d",$yyyy,$mm); # e.g.: Desktop/seismo/REA/ASNE_/1996/04/01-2334-21L.S199604
my $sfiledir = sprintf("$SEISAN_TOP/REA/ASNE_/%4d/%02d",$yyyy,$mm); # e.g.: Desktop/seismo/REA/ASNE_/1996/04/01-2334-21L.S199604
make_path($sfiledir) unless (-e $sfiledir);
my $sfile = sprintf("$sfiledir/%02d-%02d%02d-%02dL.S%04d%02d",$dd,$hh,$mi,$ss,$yyyy,$mm); # e.g.: Desktop/seismo/REA/ASNE_/1996/04/01-2334-21L.S199604
printf("Sfile = $sfile\n");
if (-e $sfile) {
	printf("$sfile already exists. Contents are:\n\n");
	system("cat $sfile");
	printf("\n");
        if (prompt_yn("Do you want to overwrite $sfile")) {
       		# do what you do for "yes"
		printf("Deleting existing Sfile\n");
		unlink($sfile);
   	} else {
       		# do what you do for "no"
		die("Not overwriting - Quitting\n");
	}
}

open SFILE, ">$sfile" or die "ERROR: could not open $sfile - $!\n";
$class = "$class " if(length($class)==1);
die("Class must be 1 or 2 characters") if (length($class)<1 or length($class)>2);
printf(SFILE " %4d %2d%2d %02d%02d %4.1f %s                                                        1\n",$yyyy,$mm,$dd,$hh,$mi,$ss,$class);
($gsec,$gmin,$ghour,$gmday,$gmon,$gyear,$gwday,$gyday,$gisdst) = gmtime();
$gyear = $gyear + 1900;
printf(SFILE " ACTION:REG %2d-%02d-%02d %2d:%02d OP:ts   STATUS:               ",substr($gyear,2,2),$gmon+1,$gmday,$ghour,$gmin);
printf(SFILE "ID:%04d%02d%02d%02d%02d%02d     I\n",$yyyy,$mm,$dd,$hh,$mi,$ss);
if ($class eq "LV") {
	$subclass = "$subclass " if (length($subclass)==1);
	printf(SFILE " VOLC MAIN %s                                                                  3\n",$subclass);
}
printf(SFILE " STAT SP IPHASW D HRMM SECON CODA AMPLIT PERI AZIMU VELO AIN AR TRES W  DIS CAZ7\n");
printf("\n");
close(SFILE);
printf("new Sfile $sfile contents are:\n\n");
system("cat $sfile");
printf("\n");
1;

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
