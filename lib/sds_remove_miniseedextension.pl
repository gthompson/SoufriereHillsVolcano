#!/usr/bin/perl
$startdir = $ARGV[0];
if (-e $startdir) {
    print("Changing to $startdir\n");
    chdir($startdir);
}
if ($#ARGV==1) {
    $dry_run = 1;
} else {
    $dry_run = 0;
}
#use glob;
@allyears = glob("[12][0-9][0-9][0-9]");
foreach $year (@allyears) {
    print("$year\n");
    chdir($year);
    @allnets = glob("??");
    foreach $net (@allnets) {
        chdir($net);
        @allstations = glob("*");
        foreach $station (@allstations) {
            chdir($station);
            @allchannels = glob("*.D");
            foreach $channel (@allchannels) {
                system("pwd");
                chdir($channel);
                @allbadnames = glob("*.miniseed");
                foreach $badname (@allbadnames) {
                    $goodname = $badname;
                    $goodname =~s/\.miniseed//;
                    if (-e $goodname) {
                        print("Caution moving $badname to $goodname. Already exists\n");
                        $badsize = -s $badname;
                        $goodsize = -s $goodname;
                        if ($badsize > $goodsize) {
                            print("Moving $badname to $goodname, because $badsize > $goodsize\n");
                            system("mv $badname $goodname") unless $dry_run;
                        } else {
                            print("Not moving $badname to $goodname, because $badsize < $goodsize\n");
                            if ($badsize == 0) {
                                print("Removing empty file $badname\n");
                                system("rm $badname") unless $dry_run;
                            }
                            if ($badsize == $goodsize) {
                                print("Removing $badname because it appears to be same file as $goodname\n");
                                system("rm $badname") unless $dry_run;
                            }
                        }
                    } else {
                        print("Moving $badname to $goodname\n");
                        system("mv $badname $goodname") unless $dry_run;
                    }
                }
                chdir("..");
            }
            chdir("..");
        } 
        chdir("..");
    }
    chdir("..");
}
