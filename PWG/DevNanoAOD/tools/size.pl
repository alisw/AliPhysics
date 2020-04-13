#!/usr/bin/perl

use warnings;
use strict;

my $file = $ARGV[0];

my @treeInfo = `aliroot -q treeinfo.C'("$file")'`;

my @keywords = ("tracks", "header", "vertices", "v0s", "cascades", "conversionphotons", "mcparticles", "fHeaders");
my %totalSize = map { $_ => 0 } @keywords;
my %comprSize = map { $_ => 0 } @keywords;
# print %size;

my $current = undef;
for (@treeInfo) 
{
  if (/\*Br/) {
    my @elem = split ":";
    
    $elem[1] =~ s/\s//g;
    $elem[1] =~ s/\*//g;
#     print $elem[1];
    if (exists($totalSize{$elem[1]})) {
      $current = $elem[1];
#       print "Switched to $current\n";
    }
  }
  $totalSize{$current} += $1 if (defined $current && /Total  Size=\s+(\d+) bytes/);
  $comprSize{$current} += $1 if (defined $current && /File Size  =\s+(\d+)/);
}

for (@keywords) {
  $comprSize{"Total"} += $comprSize{$_};
  $totalSize{"Total"} += $totalSize{$_};
}

push @keywords, "Total";

printf "%17s: \t %10s \t         %10s \n", "Branch", "Compressed", "Uncompressed";
for (@keywords) {
  next if $comprSize{$_} == 0;
  printf "%17s: \t %10d (%3d%%) \t %12d \n", $_, $comprSize{$_}, $comprSize{$_} / $comprSize{"Total"} * 100, $totalSize{$_} ;
}

# *Br   10 :fUniqueID : UInt_t                                                 *
# *Entries :      547 : Total  Size=       2761 bytes  File Size  =        124 *
# *Baskets :        1 : Basket Size=      32000 bytes  Compression=  18.28     *
# *............................................................................*
