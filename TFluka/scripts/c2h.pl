#!/usr/bin/perl

# Written by Endre Futo
# 21-May-2002

# read the filename from the command argument
$inputfn = shift @ARGV;
$len = length($inputfn);
if ($len == 0) {
   die "\nUsage: perl c2h \\(commonname\\)\n";
}
open(in, $inputfn) || die "\nCannot open the input file\n";

# construct the output file name
$outputfn = $inputfn;
# chop last character
$rightpar = chop($outputfn);
# translate to lower case
$commonname = $outputfn;
$commonname =~ tr/(/F/;
$commonname = join("",$commonname,"_H");
$outputfn =~ tr/A-Z/a-z/;
# translate '(' to 'F'
$outputfn =~ tr/(/F/;
# append '.h' to the output file name and prefix it by '>'
$outputfn = join("",">",$outputfn,".h");
open(out,$outputfn);
print out "#ifndef ",$commonname,"\n";
print out "#define ",$commonname," 1 \n";
$newline = join("","#include \"cfortran.h\"\n");
print out $newline;
$newline = join("","#include \"Rtypes.h\" \n");
print out $newline;
$newline = join("","#include \"Fdimpar.h\" \n");
print out $newline;
$newline = join("","extern \"C\" {\n");
print out $newline;

# Print the hints
# ===============
$newline = join ("","\nCheck in the input file",$inputfn,":\n");
print $newline;
print "1. At the end of the input file must be at least one blank line\n";

$newline = $outputfn;
$newline =~ tr/>/ /;
$newline = join("","\nIn the C++ header file",$newline," created from the FORTRAN common file ",$inputfn," always check: \n");
print $newline;
print "1. Mutidimensional arrays - swap the dimensions\n";
print "2. Arrays dimensioned as X(x:y) - should became X[x-y+1]\n";
print "3. CHARACTER* variables - swap the dimension and CHARACTER* size\n";
print "4. LOGICAL variables - all should be integers\n";
print "5. All comment lines\n";
print "6. All continuation lines\n";
print "7. All unions created from EQUIVALENCEs\n";
print "8. All double constants - exponent should be e not d\n";
print "9. All constants with exponent- exponent should not start with 0\n";

$oncommon = 0;
$firstcom = 0;

# Loop over lines read from the input file
# ========================================
while($line = <in>) {

# Translate everything to lower case
# ==================================

   $line =~ tr/A-Z/a-z/;
#---------------------------------------------------------------------------

# Treat comments
# ==============

   if ($line =~ /^\*/) {
      $newline = join("","\/\/",$line);
      print out $newline;
   }
#---------------------------------------------------------------------------

# Treat equivalence (became unions)
# ================================

   elsif ($line =~ /equivalence/) {
# chop last two characters
      $rightpar = chop($line);
      $rightpar = chop($line);
# search for '('
      $pos = index($line,'(');
# shift out everything before '(' inclusive
      $newline1 = substr($line,$pos+1);
# split the rest of the line according to blank
      @vars = split(" ",$newline1);
# join the line again - the line will now be without blanks
      $newline1 = join("",@vars);
# translate '(' to '['
      $newline1 =~ tr/(/[/;
# translate ')' to ']'
      $newline1 =~ tr/)/]/;
# split the rest of the line according to comma
      ($var1,$var2) = split(",",$newline1);
# here may come handling of variables
#
# first variable is an integer
      if ($var1 =~ /^[i-n]/) {
         if ($var2 =~ /^[i-n]/) {
            $newline = join("","union { int    ",$var1,";"," int    ",$var2,";};\n");
         }
         else {
            $newline = join("","union { int    ",$var1,";"," double ",$var2,";};\n");
         }
      }
# first variable is a double precision
      else {
         if ($var2 =~ /^[i-n]/) {
            $newline = join("","union { double ",$var1,";"," int    ",$var2,";};\n");
         }
         else {
            $newline = join("","union { double ",$var1,";"," double ",$var2,";};\n");
         }
      }
      print out $newline;
   }
#---------------------------------------------------------------------------

# Treat parameters (became constants)
# ===================================

   elsif ($line =~ /parameter/) {
# chop last two characters
      $rightpar = chop($line);
      $rightpar = chop($line);
# search for '('
      $pos = index($line,'(');
# shift out everything before '(' inclusive
      $newline1 = substr($line,$pos+1);
# store the variable name
      ($var) = split("=",$newline1);
# split the variable name according to blank
      @vars = split(" ",$var);
# join the line again - the variable name will now be without blanks
      $var = join("",@vars);
# search for '='
      $pos = index($line,'=');
# shift out everything before '=' inclusive
      $newline1 = substr($line,$pos+1);
# split the rest of the line according to blank
      @vars = split(" ",$newline1);
# join the line again - the line will now be without blanks
      $newline1 = join("",@vars);

# parameter is an integer
      if ($var =~ /^[i-n]/) {
         $newline = join("","const int ",$var," = ",$newline1,";\n");
         print out $newline;
      }
# parameter is real in double precision
      else {
         $newline = join("","const double ",$var," = ",,$newline1,";\n");
         print out $newline;
      }
   }
#---------------------------------------------------------------------------

# Treat commons (became struct)
# =============================

# first line of the common
# ------------------------
   elsif ($line =~ /common/) {
# finish the previous common - for the first one there is nothing to finish
      if ($firstcom eq 1) {
         $newline = join("","} ",$comname,"Common;","\n");
         print out $newline;
         $comnamebig = $comname;
# translate common name back to capitals
         $comnamebig =~ tr/a-z/A-Z/;
         $newline = join("","#define ",$comnamebig," COMMON_BLOCK(",$comnamebig,",",$comname,")\n");
         print out $newline;
         $newline = join("","COMMON_BLOCK_DEF(",$comname,"Common",",",$comnamebig,");\n");
         print out $newline;
      }
# from now on the first common is over
      $firstcom = 1;
      $newline = join("","\ntypedef struct {","\n");
      print out $newline;
# we are inside a common
      $oncommon = 1;
# look for the first '/'
      $pos = index($line,'/');
# shift out everything before first '/' inclusive
      $newline1 = substr($line,$pos+1);
      $newline2 = $newline1;
# eliminate the second '/'
      $newline2 =~ tr/\// /;
# determine the common name
      ($comname) = split(" ",$newline2);
# look for the second '/'
      $pos = index($newline1,'/');
# shift out everything before second '/' inclusive
      $newline2 = substr($newline1,$pos+1);
# split the rest of the line according to blank
      @vars = split(" ",$newline2);
# join the line again - the line will now be without blanks
      $newline2 = join("",@vars);
# here may come handling of multidimensional arrays
#
# translate comma to blank
      $newline2 =~ tr/,/ /;
# translate '(' to '['
      $newline2 =~ tr/(/[/;
# translate ')' to ']'
      $newline2 =~ tr/)/]/;
# split the rest of the line according to blank
      @vars = split(" ",$newline2);
# loop over common variables
      for ($i=0; $i < @vars; $i++) {
         if ($vars[$i] =~ /^[i-n]/) {
            $newline = join("","   int    ",$vars[$i],";\n");
            print out $newline;
         }
         else {
            $newline = join("","   double ",$vars[$i],";\n");
            print out $newline;
         }
      }
   }

# continuation line of the common
# -------------------------------
   elsif (($line =~ /&/) && ($oncommon eq 1)) {
# look for '&'
      $pos = index($line,'&');
# shift out everything before '&' inclusive
      $newline1 = substr($line,$pos+1);
# split the rest of the line according to blank
      @vars = split(" ",$newline1);
# join the line again - the line will now be without blanks
      $newline2 = join("",@vars);
# here may come handling of multidimensional arrays
#
# translate comma to blank
      $newline2 =~ tr/,/ /;
# translate '(' to '['
      $newline2 =~ tr/(/[/;
# translate ')' to ']'
      $newline2 =~ tr/)/]/;
# split the rest of the line according to blank
      @vars = split(" ",$newline2);
# loop over common variables
      for ($i=0; $i < @vars; $i++) {
         if ($vars[$i] =~ /^[i-n]/) {
            $newline = join("","   int    ",$vars[$i],";\n");
            print out $newline;
         }  
         else {
            $newline = join("","   double ",$vars[$i],";\n");
            print out $newline;
         }
      }
   }

# the line does not belong to the common (end of common)
# ------------------------------------------------------
   elsif ($oncommon eq 1) {
      $newline = join("","} ",$comname,"Common;","\n");
      print out $newline;
      $comnamebig = $comname;
# translate common name back to capital
      $comnamebig =~ tr/a-z/A-Z/;
      $newline = join("","#define ",$comnamebig," COMMON_BLOCK(",$comnamebig,",",$comname,")\n");
      print out $newline;
      $newline = join("","COMMON_BLOCK_DEF(",$comname,"Common",",",$comnamebig,");\n");
      print out $newline;
      $oncommon = 0;
   }
#---------------------------------------------------------------------------

# all other lines are just repeated as comments
# ---------------------------------------------
   else {
# prefix line with '//'
      $newline = join("","\/\/",$line);
# translate common name back to capitals
      $newline =~ tr/a-z/A-Z/;
      print out $newline;
   }
 

}
# closing curly parenthesis od the extern "C" { statement
print out "}\n";
print out "#endif\n";

close(in);
close(out);


