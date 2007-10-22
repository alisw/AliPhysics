#!/usr/bin/perl

my ($BINDIR, $BINNAME) = $0 =~ m!^(.*)/(.*)!;
if ($#ARGV != 1) {
  print STDERR <<"fnord";
usage: $BINNAME <module> <class-base>
   eg: $BINNAME base     QuadSet
       $BINNAME gl       QuadSet
       $BINNAME ged      QuadSet
       $BINNAME gedsubed QuadSet
Should be in module sub-directory (Reve/ or Alieve/).
Note that GL and Editor suffixes are not present!
fnord
  exit 1;
}

my $MODULE = $ARGV[0];
# Flat structure now.
# die "'$MODULE' not a directory" unless -d $MODULE;

%suff = ( 'gl' => 'GL', 'ged' => 'Editor', 'gedsubed' => 'Editor');
my $STEM  = $ARGV[1];
my $CLASS = $STEM . $suff{$MODULE};

if ($MODULE eq 'gedsubed') {
  $replace_xxclass = 1;
  $XXCLASS = $STEM . 'SubEditor';
}

# Flat structure now.
# my $H_NAME = "$MODULE/$CLASS.h";
# my $C_NAME = "$MODULE/$CLASS.cxx";
my $H_NAME = "$CLASS.h";
my $C_NAME = "$CLASS.cxx";

die "File '$H_NAME' already exists" if -e $H_NAME;
die "File '$C_NAME' already exists" if -e $C_NAME;

sub find_skel {
  my ($stem, $suff) = @_;
  my $file1 = "$stem-$MODULE.$suff";
  return $file1 if -e $file1;
  my $file2 = "$stem.$suff";
  return $file2 if -e $file2;
  die "Skeleton file not found neither as '$file 1' nor '$file2'";
}

my $SKEL_H_NAME = find_skel(".SKEL", "h");
my $SKEL_C_NAME = find_skel(".SKEL", "cxx");

print "Using skeleton files '$SKEL_H_NAME' and '$SKEL_C_NAME'\n";

my ($skel_h, $skel_c);
{
  my $ex_sla = $/; undef $/;
  open H, "$SKEL_H_NAME" or die "can't open $SKEL_H_NAME";
  $skel_h = <H>; close H;
  open C, "$SKEL_C_NAME" or die "can't open $SKEL_C_NAME";
  $skel_c = <C>; close C;
  $/ = $ex_sla;
}

print "Replacing CLASS -> $CLASS, STEM -> $STEM.\n";

for $f ($skel_h, $skel_c) {
  $f =~ s/XXCLASS/$XXCLASS/g if ($replace_xxclass);
  $f =~ s/CLASS/$CLASS/g;
  $f =~ s/STEM/$STEM/g;
}

print "Writing files '$H_NAME', '$C_NAME'.\n";

open H, ">$H_NAME" or
  die "can't open $H_NAME";
print H $skel_h; close H;
open H, ">$C_NAME" or
  die "can't open $C_NAME";
print H $skel_c; close H;

if($DUMPINFO) {
  print "# Now you should also edit the link-def file.\n";
  print "# Remember to run 'make depend'.\n";
}
