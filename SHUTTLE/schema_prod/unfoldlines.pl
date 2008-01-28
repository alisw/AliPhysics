#!/usr/bin/perl
# unfoldlines.pl is a script to unfold lines when they start with a leading
# space. It is meant to correct the output of ldapsearch making it comparable
# to the ldif files 
# Usage: ldapsearch [...args...] | ./unfoldlines.pl > verify.out
# if verify.out is the file to be diffed against the ldif files
# 

my @lines=<>;

$folded=0;
$cnt=1;
while ($lines[$cnt]) {
  $line = $lines[$cnt];
  if ($line =~ /^ /) {
    $folded = 1;
    $line =~ s/^ //;
    chomp $lines[$cnt-1];
    $line = $lines[$cnt-1] . $line;
    print $line;
  } else {
    if($folded==0) {print $lines[$cnt-1];}
    $folded=0;
  }
  $cnt++;
}

close(FILE);

