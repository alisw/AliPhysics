#!/usr/bin/perl -w
# runEMCALQA.pl v.0.5 A.SHABETAI alexandre.shabetai@cern.ch

use Cwd 'abs_path';
use File::Basename;
use File::Path qw(mkpath );
use File::Copy;
use Getopt::Long;

use strict;
use warnings;

##########################################################
# What this script is doing and how to use it

sub usage()
{
 print STDERR << "EOF";
 This script runs the EMCAL QA macros on selected runs.

 usage: $0 --period <PERIOD> --pass <PASS> [ --gridprefix <GRIDPREFIX> ] [--extraGridPath <EXTRAGRIDPATH> ] [ --runlist <runlist.txt> ] [ --col <pp||PbPb> ] 

 -h --help       : shows this (help) message
 
The following parameters are mandatory: 
 --period        : the running period to consider 				                                            				  (ex: LHC1XX                  )
 --pass          : the production pass to use (use "simu" for simulated data)                                                				  (ex: pass1                   )    

 You may specify the following options if you want to change their default values: 
 --runlist       : a text file with the list of runs to consider - one run by line   (default: <period>-<pass>-runlist.txt                              ) (ex: MyLHC13arunlist.txt     )
 --extraGridPath : extra grid path to use                                            (default: none                                                     ) (ex: ESDs/QA68               )
 --localpath     : path where the existing local QAresults files are stored          (default: <period>/<pass>                                          ) (ex: <a_prefix>              )
 --savePlots     : save the plots as image files (0: no ; 1: png ; 2: png + pdf)     (default: 1 (png)                                                  ) (ex: 0                       )   
 --debug         : debug flag (0 or 1)                                               (default; 0                                                        ) (ex: 1                       )
 --gridQAname    : basename QA filename on the grid (without the .root suffix)	     (default: QAresults                                                ) (ex: QAresults               )    
 --gridTimeOut   : timeOut for grid acces				             (default: 10 s                                                     ) (ex: 10                      )
 --gridprefix    : the beginning of the grid path to use                             (default: /alice/data/<year>/<period> - or .../sim/... if pass=simu) (ex: /alice/data/201X/LHC1XX )
 --filter        : filter the input file                                             (default: 0 (no)                                                   ) (ex: 1		       )
 --treeSelections: filename of file containing additional cuts on the trending tree  (default: none                                                     ) (ex: MySelection.C           )
 --trigger       : force to use the specified trigger name                           (default: trigger names are automaticly read from the input files  ) (ex: trigEMC                 )
 --col	         : force the type of collisions pp or PbPb                           (default: pp - or PbPb if period ends with "h"                     ) (ex: PbPb                    )      
                         
 NB: all switches can be abbreviated, ex: --period is valid as well as -pe

  example: $0  --period LHC13a --pass pass1 
EOF

exit;
}

######################################################################################
# the values below should not be changed (please use the command line options instead)

my $help;
my $period;
my $pass;
my $gridprefix;
my $localpath;
my $runlist;
my $col;
my $extragridPath;
my $savePlots = 0;
my $debug = 0;
my $filter = 0;
my $trigger = ""; 
my $GridTimeOut = 10;
my $gridQAname = "QAresults";
my $treeSelections = ""; 

######################################################################################
# get options 

GetOptions ( 
 "help!" => \$help,
 "period=s"  => \$period,                #string
 "pass=s" => \$pass,                     #string
 "gridprefix=s" => \$gridprefix,         #string
 "localpath=s" => \$localpath,           #string
 "runlist=s" => \$runlist,               #string
 "col=s" => \$col,                       #string
 "extraGridPath=s" => \$extragridPath,   #string
 "savePlots=i" => \$savePlots,           #int
 "debug=i" => \$debug,                   #int
 "filter=i" => \$filter,                 #bool
 "triggerForced=s" => \$trigger,
 "gridTimeOut=i" => \$GridTimeOut,       #int
 "gridQAname=s" => \$gridQAname,         #string
 "treeSelections=s" => \$treeSelections) #string
 or die("Error in command line arguments\n");
 usage() if ($help);
 die "Missing --period! (you can get help with $0 -h)" unless $period;
 die "Missing --pass! (you can get help with $0 -h)" unless $pass;

if ($pass eq "sim") {$pass= "simu";}
if (!(defined($runlist))) {$runlist = "$period-$pass-runlist.txt";}
if (!(defined($col)) && index($period, 'h') != -1) {$col = "PbPb";}
if (!(defined($col))) {$col = "pp";}

if (!(defined($gridprefix))) {
 my $year=sprintf("20%s",substr($period, 3, 2));
 if ($pass ne "simu") {$gridprefix = "/alice/data/$year/$period";} else {$gridprefix = "/alice/sim/$year/$period";}
}
if (!(defined($extragridPath))) {
  if    (($pass ne "simu") && (index($pass, 'c') != -1)) { $extragridPath = $pass;}
  elsif ($pass ne "simu") { $extragridPath = "ESDs/$pass";}
  else  {$extragridPath = "";}
 }

######################################################################################
# the values below should not be changed

my $macroDir = abs_path(dirname($0));
my $locDir = "$period/$pass";
$ENV{QAPATH} = "$locDir/";
my $processedrunlist = "$locDir/processed.txt";
my $notprocessedrunlist = "$locDir/notprocessed.txt";
my $exit_code = 0;

my $calo = "EMCAL";
my @gridsuffixes = qw(.root);
if (index($pass, 'c') != -1) {@gridsuffixes = qw(_barrel.root _outer.root);}
if (index($pass, 'muon_calo_') != -1) {@gridsuffixes = qw(.root);}
if (!(defined($localpath))) {$localpath = $locDir;} else {$localpath = $localpath."/";}

######################################################################################
# sanity checks

if (! -e $locDir) {mkpath($locDir);}
if ((defined(glob("$locDir/*data.txt")))) {unlink(glob("$locDir/*data.txt"));}
my $includePath = qx{aliroot -b -l << EOF
.include
EOF
};
if (index($includePath,"$ENV{ALICE_ROOT}/include") == -1) 
  {print "Please add \$ALICE_ROOT/include to your include path using a rootlogon.C of .rootrc\n"; exit;}

my $ALIEN_RUNTIME_ROOT = $ENV{'ALIEN_RUNTIME_ROOT'};
if (!(defined($ALIEN_RUNTIME_ROOT))) {die "ALIEN_RUNTIME_ROOT not defined!";}
my $proxyValidity = qx{$ALIEN_RUNTIME_ROOT/bin/xrdgsiproxy info 2>&1 | grep "time left" | cut -d " " -f 6 | cut -d: -f1}; chomp($proxyValidity);
if (!(defined($proxyValidity)) || $proxyValidity eq "0h") {  die "No valid proxy found (or proxy valid for less that one hour). Nothing done!"};
my $isValidToken = qx{$ALIEN_RUNTIME_ROOT/bin/alien-token-info | grep -c "Token is still valid"};
if ($isValidToken==0) { die "No valid token found. Nothing done!";}
######################################################################################
# start working

open(STDOUT, "| tee -i \"$locDir/$period\_$pass\_EMCALQA.log\"");
open(STDERR, '>&STDOUT');
open(FIN, $runlist) or die "Could not open $runlist: $! please specify --runlist of create $runlist";
open(FOUT, ">$processedrunlist") or die "Could not open $processedrunlist: $!";
open(FOUT2, ">$notprocessedrunlist") or die "Could not open $notprocessedrunlist: $!";

print "\n Executing $0";
 foreach (@ARGV)  { print " $_";} 
print " ...\n";

# read the run list
my @runs;
while( my $line = <FIN>)  {
    chomp($line);  
    my @fields = split(/,/, $line);
    push @runs, @fields
}
    @runs = sort { $a <=> $b } @runs;

# create the directory structure and copy the output of the QA train from the grid
    foreach my $run (@runs) {   
    $run =~ tr/ //ds;
    if($pass ne "simu") {$run = sprintf("%09d", $run);}
    print "\nProcessing run $run...\n";
    foreach my $suffix (@gridsuffixes) {
     my $gridFile = "alien://$gridprefix/$run/$extragridPath/$gridQAname$suffix"; 
     my $shortrun = $run; $shortrun =~ s/^0+//; $shortrun =~ s/ +$//; 
     my $locFileName = "$shortrun$suffix"; 
    
     print "\n Processing file $locFileName...\n";
    
     if (! -e "$locDir/$run") {mkdir("$locDir/$run");} #

if ($localpath eq $locDir) {
     for (my $trial = 1; $trial <= 2; $trial++)  
     {  
 
      system("root -b -l -q \'$macroDir/CopyQAFile.C+g\(\"$gridFile\",\"$localpath\",kFALSE,\"$locFileName\",$GridTimeOut\)\'");  
      $exit_code = $?>>8;
    
     if (! $exit_code) { print "\nFile $gridQAname$suffix processed succesfully (run $shortrun, try $trial)\n";last} else { print "\nThe processing of file $gridQAname$suffix (run $shortrun) failed! (try $trial)\n";};
     if ($exit_code eq 252 || $exit_code eq 253) {unlink("$localpath/$locFileName"); }
    }  
 }
# run the run level QA for the runs that were selected (i.e. in the runlist and that do have a valid local QA output file)

 if (!$exit_code) {

  if (-e "$localpath/$locFileName" )  { 

   system("aliroot -b -l -q \'$macroDir/CreateEMCALRunQA.C+g\(\"$localpath/$locFileName\",\"$run\",\"$period\",\"$pass\",$savePlots,$filter,\"$trigger\",\"$col\",\"$calo\"\)\'"); 
   $exit_code = $?>>8;
   
   if (! $exit_code) {  print FOUT "$shortrun\n"; print "\nCreateEMCALRunQA() for file $gridQAname$suffix processed succesfully (run $shortrun)\n";} else { print FOUT2 "$shortrun\n" ; print "\nThe processing of CreateEMCALRunQA() for file $gridQAname$suffix failed! (run $shortrun)\n";};
  
 if ($suffix ne ".root") { move("$locDir/$run/trending.root", "$locDir/$run/trending$suffix"); move("$locDir/$run/${period}_${pass}_${shortrun}_QAplots.root","$locDir/$run/${period}_${pass}_${shortrun}_QAplots$suffix");}

   } else { print " File not found: $localpath/$locFileName";} }
  } 
}
 
close FIN;
close FOUT;
close FOUT2;
$exit_code = 0;

######################################################################################
# merge the QA trees

if ((! $exit_code) && defined(glob("$locDir/*/trendin*.root")) && (-e glob("$locDir/*/trendin*.root")))
{

 print "\nMerging....\n";	
 
 system("hadd -v 1 -f $locDir/trending.root $locDir/*/trendin*.root");
 $exit_code = $?>>8;

 if (! $debug) { system("rm $locDir/*/trendin*.root"); }

} 
######################################################################################
# process and plot the period level QA 

print "\n";

if ((! $exit_code) && (-e "$locDir/trending.root")) { system("root -b -l -q \'$macroDir/PlotEMCALQATrendingTree.C(\"$locDir/trending.root\",$savePlots,\"$treeSelections\",\"$trigger\"\)\'"); }

$exit_code = $?>>8;
if (! $exit_code) { print "\nPlotEMCALQATrendingTree() processed succesfully\n";} else { print "\nThe processing of PlotEMCALQATrendingTree() failed!\n";};
   
######################################################################################
# cleanup

move($processedrunlist,"$locDir/runlist.txt")  or die "Cannot move the runlist to the correct local directory : $!";
move("$locDir/trendingPlots.root","$locDir/$period"."_$pass"."_trendingPlots.root");

######################################################################################
# final merging

if ((! $exit_code) && defined(glob("$locDir/*/${period}_${pass}_*_QAplot*.root $locDir/${period}_${pass}_trendingPlots.root $locDir/trending.root")) && (-e glob("$locDir/*/${period}_${pass}_*_QAplot*.root $locDir/${period}_${pass}_trendingPlots.root $locDir/trending.root")))
{

 print "\n Final Merging....\n";

 system("hadd -v 1 -f $locDir/${period}_${pass}_EMCALQA.root  $locDir/*/${period}_${pass}_*_QAplot*.root $locDir/${period}_${pass}_trendingPlots.root $locDir/trending.root");

 if (! $debug) {system("rm -f  $locDir/*/${period}_${pass}_*_QAplot*.root $locDir/${period}_${pass}_trendingPlots.root"); }

}

######################################################################################
# cleanup some files used for debugging

if (! $debug && (defined(glob("$locDir/*/*.txt")))) {unlink(glob("$locDir/*/*.txt"));}
if (! $debug && (defined(glob("$locDir/*.txt")))) {unlink(glob("$locDir/*.txt"));}
if (! $debug && ! $savePlots && (defined(glob("$locDir/*/"))))  {system("rmdir --ignore-fail-on-non-empty  $locDir/*/");}


######################################################################################
#| make final pdf

if ((! $exit_code) && -e "$locDir/${period}_${pass}_EMCALQA.root") { system("root -b -l -q \'$macroDir/MakeQAPdf.C+g(\"$locDir/${period}_${pass}_EMCALQA.root\"\)\'"); }


######################################################################################
# that's all folks 

print "\n Done! \n";


