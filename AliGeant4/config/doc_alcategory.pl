#!/usr/local/bin/perl
# $Id$
# Ivana Hrivnacova 15.2.00
#
# This script generates the class category page
# from all *.h files found in the category directory

# no test of arguments is performed
$CAT = $ARGV[0];
$NAME = $ARGV[1];

# open output file
$output =  $ENV{'AG4_INSTALL'} . "/doc/HTML/" . $NAME . "Category.html";
open(OUTPUT, ">" . $output);

print "Processing class category: " . $NAME . "\n"; 

# print the begining of file
print OUTPUT "<HTML>\n";
print OUTPUT "\n";
print OUTPUT "<HEAD>\n";
print OUTPUT "<TITLE>Class category: ". $NAME . "</TITLE></HEAD>\n";
print OUTPUT "\n";
print OUTPUT "\n";
print OUTPUT "<BODY bgcolor=#FFFFFF>\n";
print OUTPUT "\n";
print OUTPUT "<!-- Header material -->\n";
print OUTPUT "<table border=0   cellpadding=5 cellspacing=0 width=\"100%\">\n";
print OUTPUT "  <tr bgcolor=#d0ffd0>\n";
print OUTPUT "     <td align=left width=30%>\n";
print OUTPUT "     <img alt=\"Alice\"\n";
print OUTPUT "       src=\"http://AliSoft.cern.ch/offline/geant4/gif/AliceLogo.gif\"\n";
print OUTPUT "	    width=\"60\" height=\"60\" align=\"absmiddle\" border=1>\n";
print OUTPUT "     <td align=center width=40%>\n";
print OUTPUT "        <font size=\"+2\">\n";
print OUTPUT "           Alice Geant4 Simulation Code Prototype        </font>\n";
print OUTPUT "     <td align=right width=30% valign=bottom>\n";
print OUTPUT "        <font size=\"-1\">\n";
print OUTPUT "        <script language=\"JavaScript\">\n";
print OUTPUT "        document.write(\"Last modified \"+ document.lastModified)\n";
print OUTPUT "        // end of script -->\n";
print OUTPUT "        </script></font>\n";
print OUTPUT "     </td>\n";
print OUTPUT "  </tr>\n";
print OUTPUT "</table>\n";
print OUTPUT "<CENTER>\n";
print OUTPUT "<H2>Class category: " . $NAME . "</H2>\n";
print OUTPUT "</CENTER>\n";
print OUTPUT "\n";
print OUTPUT "<P><HR SIZE=5><BR>\n";
print OUTPUT "\n";
print OUTPUT "<UL><BR>\n";
print OUTPUT "\n";
print OUTPUT "<LI><STRONG>C++ header files:</STRONG>\n";
print OUTPUT "\n";
print OUTPUT "  <UL>\n";

# print the linked header files
$CATSTRING = "\"Category: " . $CAT . "\"";
@FILELIST = `find . -name \"*.h\" -exec grep -l  $CATSTRING  {} \\;`;

foreach $FILEPATH (@FILELIST) { 
  @TEMP = split('/',$FILEPATH);
  $FILE = @TEMP[@TEMP - 1];
  chop $FILE;
  print "   Linking file " . $FILE . "\n"; 
  print OUTPUT "  <LI><A HREF=\"" . $FILE . ".html\">" . $FILE . "</A>\n";
}

# print the end of file
$today = localtime(time);
$today =~ s/ \d\d:\d\d:\d\d / /;
@list = getpwuid($<);
$user = $list[6];
print OUTPUT "</UL>\n";
print OUTPUT "\n";
print OUTPUT "</UL>\n";
print OUTPUT "\n";
print OUTPUT "<P><HR SIZE=5>\n";
print OUTPUT "\n";
print OUTPUT "<ADDRESS>\n";
print OUTPUT "Created on $today by <B>$user</B> <BR>\n";
print OUTPUT "using the HTML generator\n";
print OUTPUT "<A HREF=\"http://home.cern.ch/~binko/Ddl2Html/Ddl2Html.html\">Ddl2Html description</A>\n";
print OUTPUT " (the source <A HREF=\"http://home.cern.ch/~binko/Ddl2Html/Ddl2Html.code\">Perl5 code</A>)\n";
print OUTPUT "</ADDRESS>\n";
print OUTPUT "\n";
print OUTPUT "</BODY bgcolor=#FFFFFF >\n";
print OUTPUT "\n";
print OUTPUT "</HTML>\n";

# close output file
close(OUTPUT);    
