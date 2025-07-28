#!/usr/bin/perl -w

use strict;
use File::Basename;
use File::Copy;
use Getopt::Std; 
#use Switch;
#use threads;
#use threads::shared;

my ($rs,@files)=@ARGV;
my $nparse=10;
my $nchilds=20;

#my @files=`ls $name*.root`;

unless ($#files>=0) {
  print "NO  files - quitting\n";
  exit;
}
RunDivision($rs,@files);

sub RunDivision{
    my($res,@file)=@_;
    chomp $res;
    print " entering runDivision: ${res} \n";# \n @files \n";   

    if ($#file==0){
        my $fl=$file[0];
        chomp $fl;
	my $cmd="cp $fl ${res}";
	print "$cmd \n";
	`$cmd`;
	return 0;
    }

    my $outName="${res}_t";
    my $count=0;
    my @childs;
    while (@file){
	my $n=0;
	my $list=" ";
	while( @file && ($n<$nparse)){
	    my $f=pop @file;
	    chomp $f;
	    $list="$list $f";
	    $n++;
	}
	my $cmd="";
	if ($n==1){
	    $cmd="cp $list ${outName}$count.root";
	}
	else { 
        #start new thread
        $cmd="hadd -f ${outName}$count.root $list";
       }
	my $pid = fork();
        if ($pid) {
	    # parent
	    #print "pid is $pid, parent $$\n";
	    push(@childs, $pid);
        } elsif ($pid == 0) {
	    # child
	    runCmd($cmd);
	    exit 0;
        } else {
	    die "couldnt fork: $!\n";
        }

	if ($#childs == $nchilds){
	    foreach (@childs) {
		my $tmp = waitpid($_, 0);
		print "done with pid $tmp\n";
	    }
	    @childs=();
	}

	$count++;
    }
    
    foreach (@childs) {
	my $tmp = waitpid($_, 0);
	print "done with pid $tmp\n";
    }
     my @newlist=`ls ${outName}*.root`;
    RunDivision("${outName}.root",@newlist);
    my $cmd="mv ${outName}.root ${res}";
    print "$cmd \n";
    `${cmd}`;
    `rm ${outName}*.root`;
}

sub runCmd{
    my $cmd=shift;
    print "$cmd \n";
    my $tmp=`$cmd`;
    print "$tmp \n";
}
