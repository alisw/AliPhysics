#!/usr/local/bin/perl5
# $Id$
################################################################################
#
# Automatic HTML generator for projects using Objectivity
#
# Author : Pavel Binko
#
# Last update : 09/04/96
#
################################################################################
# Modified for Alice specifics by Ivana Hrivnacova: 22.5.98;
# some more changes 12.2.99;
#
# ------------------------------------------------------------------------------
# Analyse the command arguments
#
if( $#ARGV == -1 )
{
  $pwd = `pwd`;
  chop($pwd);
  @help = split(/\//,$pwd);
  $help[$#help] =~ tr/a-z/A-Z/;
  $html_title = "The ".$help[$#help]." Project";
  undef $help;
}
elsif( $#ARGV == 0 && $ARGV[0] =~ /^-[?hH]/ )
{
  print "\n";
  print "The Ddl2Html functionality :\n\n";
  print "  Ddl2Html has to be started from the project directory you'd like\n";
  print "           to document (to convert to html)\n";
  print "  Ddl2Html creates (if not existing) the directory HTML in the project\n";
  print "           directory and adds the afs access rights \"cern:nodes rl\"\n";
  print "           to both HTML and project directories\n";
  print "  Ddl2Html creates index.html file and *.ddl.html, *.h.html and\n";
  print "           *.C.html files in the HTML directory\n\n";
  print "  Contents of the project directory remains unchanged.\n\n\n";
  print "The Ddl2Html command syntax :\n\n";
  print "  Ddl2Html [ arguments ] \n\n";
  print "    No argument : The project directory name will be chosen\n";
  print "                  as the project name (in upper case)\n\n";
  print "    One argument : a) If the argument starts with -h or -H or -?,\n";
  print "                      this help text will be printed\n";
  print "                   b) If the argument matches with a file name\n";
  print "                      in the project directory, the first line\n";
  print "                      from the file will be the project title\n";
  print "                   c) Othetwise  the whole text from the command\n";
  print "                      line will be the project title\n\n";
  print "    More than one argument : The whole text from the command line\n";
  print "                             will be the project title\n\n\n";
  print "For further details see the :\n";
  die "    http://home.cern.ch/~binko/Ddl2Html/Ddl2Html.html\n\n";
}
elsif( $#ARGV == 0 && -e $ARGV[0] )
{
  $first = 1;
  while( <> )
  {
    if( $first )
    {
      chop($_);
      $html_title = $_;
      $first = 0;
    }
  }
  undef $first;
}
else
{
  $html_title = join( " ", @ARGV );
}

study();

print "\n\n";
print "*****************************************************************\n";
print "*\n";
print "*  $html_title\n";
print "*\n";
print "*****************************************************************\n\n\n";

#
# ------------------------------------------------------------------------------
# Get some information about the creator
#
# UID is in the variable $<
@list = getpwuid($<);
$user = $list[6];

#
# ==============================================================================
# Get a list of all .ddl files
#
@ddls = <*.ddl>;
print "List of Objectivity DDL files = @ddls\n\n";

#
# ------------------------------------------------------------------------------
# Get a list of .h files, where no .ddl files exist (ignore _ref.h files)
#
while( $hdrname = <*.h> )
{
  if( ! ( $hdrname =~ /_ref.h$/ ) )
  {
    $helpname = $hdrname;
    $helpname =~ s/.h$/.ddl/;
    if( ! -e $helpname )
    { @hdrs = (@hdrs, $hdrname);
    }
  }
}
print "List of C++ header files = @hdrs\n\n";

#
# ------------------------------------------------------------------------------
# Get a list of all C++ programs
#
while( $cppname = <*.C> )
{
  if( ! ( $cppname =~ /_ddl.C$/ ) )
  {
    $helpname = $cppname;
    $helpname =~ s/.C$//;
    if( -e "$helpname.ddl" )
    { @cppddls = (@cppddls, $cppname);
    }
    elsif( -e "$helpname.h" )
    { @cpphdrs = (@cpphdrs, $cppname);
    }
    else
    { @cppmain = (@cppmain, $cppname);
    }
  }
}
print "List of C++ programs to the *.ddl = @cppddls\n\n";
print "List of C++ programs to the *.h = @cpphdrs\n\n";
print "List of main programs = @cppmain\n\n";

#
# ------------------------------------------------------------------------------
# Get a list of all Makefiles, makefiles, s.Makefiles and s.makefiles
#
@maks = <*.gmk>;
print "List of Makefiles = @maks\n\n";

#
# ------------------------------------------------------------------------------
# Get a list of all config.* files
#
#@cnfs = <config.*>;
@cnfs = <setup* *boot*>;
print "List of configuration files = @cnfs\n\n";

#
# ==============================================================================
# Analyse the .ddl and .h files
#
print "Analysing the .ddl and .h files ...\n\n";

foreach $file (@ddls, @hdrs)
{
  open (HF, $file);
  $fileline = 0;

  while (<HF>)
  {
    $fileline += 1;
    chop;
    if ( ! (/^[ \t]*\/\//) )	# ignore C++ comment lines
    {
      s/\/\/.*//;		# ignore all after C++ comment sign //

      if ( /\bclass\b/ && ! ( /template.*<.*class.*>/ ) )
      {
        @words = split();			# split line read in
        for( $i = 0 ; $i < $#words ; $i++ )
        {
          if( $words[$i] eq "class" )
          {
            $words[$i+1] =~ s/:.*//;
            $words[$i+1] =~ s/{.*//;
            if( !($words[$i+1] =~ /;/ ) && !($words[$i+2] =~ /^;/ )
                                        && ($words[$i+1] ne "") )
            {
              $fileclass{$words[$i+1]} = $file."//".$fileline;
              $i = $#words;
            }
          }
        }
      }

      if ( /\benum\b/ )
      {
        @words = split();			# split line read in
        for( $i = 0 ; $i < $#words ; $i++ )
        {
          if( $words[$i] eq "enum" )
          {
            $words[$i+1] =~ s/{.*//;
            if( $words[$i+1] ne "" )
            {
              $fileenum{$words[$i+1]} = $file."//".$fileline;
              $i = $#words;
            }
          }
        }
      }

    }
  }
  close(HF);
}

# debug
# print "\n\n\nDEBUG INFORMATION :\n";
# foreach $class (keys %fileclass)
# { print "DEBUG: class $class found in $fileclass{$class}\n";
# }
# foreach $enum (keys %fileenum)
# { print "DEBUG: enum $enum found in $fileenum{$enum}\n";
# }
# print "\n\n";
# end debug

#
# ==============================================================================
# Some variables needed for html files
#
$today = localtime(time);
$today =~ s/ \d\d:\d\d:\d\d / /;

#
# ------------------------------------------------------------------------------
# Create the directory html (if not existing)
#
if( ! -e "HTML" )
{
  `mkdir HTML`;
}
if( -e "/usr/sue/bin/fs" )
{
  `/usr/sue/bin/fs setacl . cern:nodes rl`;
  `/usr/sue/bin/fs setacl HTML cern:nodes rl`;
}

#
# ==============================================================================
# Create *.ddl.html and *.h.html files
#
foreach $file (@ddls, @hdrs)
{
print "Writing file $file.html ... \n";

# Open the .ddl or .h file and the .html file
open (HDRF, $file);
open (HLPF,">HTML/$file.html");

&html_header($file);
#print HLPF "<HR SIZE=5><P>\n\n";
#&centred_header1($html_title);
#print HLPF "<P><HR SIZE=5><P>\n\n";
&see_source_file;
print HLPF "<P><HR SIZE=5><P>\n\n";
&centred_header1($file);
print HLPF "<P><HR><P>\n\n";

print HLPF "<PRE>\n\n";
&html_ddl_h_cpp_code;
print HLPF "\n</PRE>\n\n";

print HLPF "<P><HR SIZE=5><P>\n\n";
&back_to_index_see_source_file;
&create_address;

close(HLPF);
close(HDRF);
}
print "\nHTML for Objectivity DDL files and C++ header files created ... \n\n";

#
# ------------------------------------------------------------------------------
# Create *.C.html files
#
foreach $file (@cppddls, @cpphdrs, @cppmain)
{
print "Writing file $file.html ... \n";

# Open the C++ file and the .html file
open (HDRF, $file);
open (HLPF,">HTML/$file.html");

&html_header($file);
#print HLPF "<HR SIZE=5><P>\n\n";
#&centred_header1($html_title);
#print HLPF "<P><HR SIZE=5><P>\n\n";
&see_ddl_h_file;
print HLPF "<P><HR SIZE=5><P>\n\n";
&centred_header1($file);
print HLPF "<P><HR><P>\n\n";

print HLPF "<PRE>\n\n";
&html_ddl_h_cpp_code;
print HLPF "\n</PRE>\n\n";

print HLPF "<P><HR SIZE=5><P>\n\n";
&back_to_index_see_ddl_h_file;
&create_address;

close(HLPF);
close(HDRF);
}
print "\nHTML for C++ files created ... \n\n";

#
# ------------------------------------------------------------------------------
# Copy Make and Config files
#
foreach $file (@maks, @cnfs)
{
print "Copying file $file ... \n";
`cp $file HTML/$file`;
}
print "\nMake and Config files copied ... \n\n";

#
# ------------------------------------------------------------------------------
# Write index.html
#
print "Writing index file ... \n";

open (HLPF,">HTML/index.html");

&html_header($html_title);

#print HLPF "<H1 ALIGN=RIGHT>";
#print HLPF "<TABLE COLSPEC=\"20\" BORDER=5 CELLPADDING=0 CELLSPACING=0>\n";
#print HLPF "<TR><TH>";
#print HLPF "<A HREF=\"index.html\" TARGET=_top>No Frame</A>";
#print HLPF "</TH></TR>\n";
#print HLPF "</TABLE></H1>\n\n";

print HLPF "<!-- Header material -->\n";
print HLPF "<table border=0   cellpadding=5 cellspacing=0 width=\"100%\">\n";
print HLPF "  <tr bgcolor=#d0ffd0>\n";
print HLPF "     <td align=left width=30%>\n";
print HLPF "     <img alt=\"Alice\"\n"; 
print HLPF "        src=\"http://AliSoft.cern.ch/offline/geant4/gif/AliceLogo.gif\"\n"; 
print HLPF "	    width=\"60\" height=\"60\" align=\"absmiddle\" border=1>\n";
print HLPF "     <td align=center width=40%>\n";
print HLPF "        <font size=\"+2\">\n";
print HLPF "           Geant4 Project";
print HLPF "        </font>\n";
print HLPF "     <td align=right width=30% valign=bottom>\n"; 
print HLPF "        <font size=\"-1\">\n"; 
print HLPF "        <script language=\"JavaScript\">\n"; 
print HLPF "        document.write(\"Last modified \"+ document.lastModified)\n"; 
print HLPF "        // end of script -->\n"; 
print HLPF "        </script></font>\n"; 
print HLPF "     </td>\n";
print HLPF "  </tr>\n";
print HLPF "</table>\n";

#&centred_header1($html_title);
#print HLPF "<P><HR SIZE=5><P>\n\n";

print HLPF "<UL><BR>\n\n";

if( @hdrs )
{
print HLPF "<LI><STRONG>C++ header files:</STRONG>\n";
  print HLPF "  <UL>\n";
  foreach $file (@hdrs)
  { print HLPF "  <LI><A HREF=\"$file.html\">$file</A>\n";
  }

if( @cpphdrs )
{
#  &and_cpp_code;
#  foreach $file (@cpphdrs)
#  { print HLPF "  <LI><A HREF=\"$file.html\">$file</A>\n";
#  }
}
print HLPF "  </UL>\n\n";

print HLPF "<P><HR><P>\n\n";
}

if( @ddls )
{
print HLPF "<LI><STRONG>Objectivity DDL files:</STRONG>\n";
  print HLPF "  <UL>\n";
  foreach $file (@ddls)
  { print HLPF "  <LI><A HREF=\"$file.html\">$file</A>\n";
  }

if( @cppddls )
{
#  &and_cpp_code;
#  foreach $file (@cppddls)
#  { print HLPF "  <LI><A HREF=\"$file.html\">$file</A>\n";
#  }
}
print HLPF "  </UL>\n\n";

print HLPF "<P><HR><P>\n\n";
}

if( @cppmain )
{
print HLPF "<LI><STRONG>Main programs/Extern methods:</STRONG>\n";
  print HLPF "  <UL>\n";
  foreach $file (@cppmain)
  { print HLPF "  <LI><A HREF=\"$file.html\">$file</A>\n";
  }
  print HLPF "  </UL>\n\n";

print HLPF "<P><HR><P>\n\n";
}

if( @maks || @cnfs )
{
print HLPF "<LI><STRONG>Makefiles and configuration files:</STRONG>\n";
  print HLPF "  <UL>\n";
  foreach $file (@maks)
#   { print HLPF "  <LI><A HREF=\"../$file\">$file</A>\n";
  { print HLPF "  <LI><A HREF=\"$file\">$file</A>\n";
  }
  foreach $file (@cnfs)
#   { print HLPF "  <LI><A HREF=\"../$file\">$file</A>\n";
  { print HLPF "  <LI><A HREF=\"$file\">$file</A>\n";
  }
  print HLPF "  </UL>\n\n";
}

print HLPF "</UL>\n\n";

&create_address;

close(HLPF);

print "\n\n";

########################  END OF THE PROGRAM  ##################################

#
# ------------------------------------------------------------------------------
# Subroutine create_address
#
sub create_address
{
  print HLPF "<P><HR SIZE=5><BR>\n\n";

  print HLPF "<ADDRESS>\n";
  print HLPF "Created on $today by <B>$user</B> <BR>\n";
  print HLPF "using the HTML generator\n";
  print HLPF "<A HREF=\"http://home.cern.ch/~binko/Ddl2Html/Ddl2Html.html\">Ddl2Html description</A>\n";
  print HLPF " (the source ";
  print HLPF "<A HREF=\"http://home.cern.ch/~binko/Ddl2Html/Ddl2Html.code\">Perl5 code</A>)\n";
  print HLPF "</ADDRESS>\n\n";

  print HLPF "</BODY bgcolor=#FFFFFF >\n\n";

  print HLPF "</HTML>\n";
}

#
# ------------------------------------------------------------------------------
# Subroutine back_to_index
#
sub back_to_index
{
  print HLPF "<H4>Back to: ";
  print HLPF "<A HREF=\"../G4CodePrototype.html\"> Class categories index, </A>";
  print HLPF "<A HREF=\"index.html\">Alphabetical index</A>";
  print HLPF "</H4>\n\n";
}

#
# ------------------------------------------------------------------------------
# Subroutine see_source_file
#
sub see_source_file
{
  # The C++ file corresponding to the .ddl or .h file
  $cfile = $file;
  $cfile =~ s/.h$/.C/;
  $cfile =~ s/.ddl$/.C/;

  if( -e $cfile )
  {
    print HLPF "<H4>See the source file ";
    print HLPF "<A HREF=\"$cfile.html\">$cfile</A>";
    print HLPF "</H4>\n\n";
  }
  undef $cfile;
}

#
# ------------------------------------------------------------------------------
# Subroutine see_ddl_h_file
#
sub see_ddl_h_file
{
  # The .ddl or .h file corresponding to the C++ file
  $ddlfile = $file;
  $ddlfile =~ s/.C$/.ddl/;
  if( -e $ddlfile )
  {
    print HLPF "<H4>See the Objectivity DDL file ";
    print HLPF "<A HREF=\"$ddlfile.html\">$ddlfile</A>";
    print HLPF "</H4>\n\n";
  }
  else
  {
    $hfile = $file;
    $hfile =~ s/.C$/.h/;
    if( -e $hfile )
    {
      print HLPF "<H4>See the C++ header file ";
      print HLPF "<A HREF=\"$hfile.html\">$hfile</A>";
      print HLPF "</H4>\n\n";
    }
  }
  undef $ddlfile;
  undef $hfile;
}

#
# ------------------------------------------------------------------------------
# Subroutine back_to_index_see_source_file
#
sub back_to_index_see_source_file
{
  &back_to_index;
  &see_source_file;
}

#
# ------------------------------------------------------------------------------
# Subroutine back_to_index_see_ddl_h_file
#
sub back_to_index_see_ddl_h_file
{
  &back_to_index;
  &see_ddl_h_file;
}

#
# ------------------------------------------------------------------------------
# Subroutine and_cpp_code
#
sub and_cpp_code
{
  print HLPF "<P>\n";
  print HLPF "<STRONG>... and the corresponding C++ code :</STRONG>\n\n";
}

#
# ------------------------------------------------------------------------------
# Subroutine centred_header1
#
sub centred_header1
{
  local($sometitle) = @_;

  print HLPF "<CENTER>\n";
  print HLPF "<H1>$sometitle</H1>\n";
  print HLPF "</CENTER>\n\n";
}

#
# ------------------------------------------------------------------------------
# Subroutine html_header
#
sub html_header
{
  local($sometitle) = @_;

  print HLPF "<HTML>\n\n";

  print HLPF "<HEAD>\n";
  print HLPF "<TITLE>$sometitle</TITLE>";
  print HLPF "</HEAD>\n\n";

  print HLPF "<BODY bgcolor=#FFFFFF>\n\n";
}

#
# ------------------------------------------------------------------------------
# Subroutine html_ddl_h_cpp_code
#
sub html_ddl_h_cpp_code
{
$fileline = 0;

while (<HDRF>)
{
  $fileline += 1;
  chop;
		
  s/</&lt\;/g;			# convert special characters to html
  s/>/&gt\;/g;

# Write the HTML
  foreach $class (keys %fileclass)		# add links for classes 
  {
    undef $newstr;
    undef $newerstr;

    ( $locfileclass, $locfileline ) = split( /\/\//, $fileclass{$class} );

    if( ($file eq $locfileclass) && ($fileline eq $locfileline) )
    { print HLPF "<A NAME=\"$class\_classdef\">\n";
    }

    if ( /\b$class\b/ )
    {
      if( $file eq $locfileclass )
      { $newstr="<a href=\"\#$class\_classdef\">$class</a>";
      }
      else
      { $newstr="<a href=\"$locfileclass.html\#$class\_classdef\">$class</a>";
      }

undef $hotovo;

while( $_ ne "" )
{
  if( /\b$class\b/ )
  {
    $hotovo = $hotovo.$`;
    $zbytek = $';

    if( !( ($hotovo =~ /<a href="$/) && ($zbytek =~ /^\./) ) )
    {
      $hotovo = $hotovo.$newstr;
    }
    else
    {
      $hotovo = $hotovo.$class
    }

    $_ = $zbytek;
    /\b$class\b/;
  }
  else
  {
    $hotovo = $hotovo.$_;
    $_ = "";
  }
}

$_ = $hotovo;

      if ( /$newstr\.h/ )
      {
        if( -e "$class.ddl" )
        { $newerstr="<a href=\"$class.ddl.html\">$class.h</a>";
          s/$newstr\.h/$newerstr/g; 
        }
        else
        { $newerstr="<a href=\"$class.h.html\">$class.h</a>";
          s/$newstr\.h/$newerstr/g; 
        }
      }

      if ( /$newstr\.C/ )
      { $newerstr="<a href=\"$class.C.html\">$class.C</a>";
        s/$newstr\.C/$newerstr/g; 
      }

    }
  }

  foreach $enum (keys %fileenum)		# add links for enums
  {
    undef $newstr;
    undef $newerstr;

    ( $locfileenum, $locfileline ) = split( /\/\//, $fileenum{$enum} );

    if( ($file eq $locfileenum) && ($fileline eq $locfileline) )
    { print HLPF "<A NAME=\"$enum\_enumdef\">\n";
    }

    if ( /\b$enum\b/ )
    {
      if( $file eq $locfileenum )
      { $newstr="<a href=\"\#$enum\_enumdef\">$enum</a>";
        s/\b$enum\b/$newstr/g; 
      }
      else
      { $newstr="<a href=\"$locfileenum.html\#$enum\_enumdef\">$enum</a>";
        s/\b$enum\b/$newstr/g; 
      }

      if ( /$newstr\.ddl/ )
      { $newerstr="<a href=\"$enum.ddl.html\">$enum.ddl</a>";
        s/$newstr\.ddl/$newerstr/g; 
      }

      if ( /$newstr\.h/ )
      {
        if( -e "$enum.ddl" )
        { $newerstr="<a href=\"$enum.ddl.html\">$enum.h</a>";
          s/$newstr\.h/$newerstr/g; 
        }
        else
        { $newerstr="<a href=\"$enum.h.html\">$enum.h</a>";
          s/$newstr\.h/$newerstr/g; 
        }
      }

      if ( /$newstr\.C/ )
      { $newerstr="<a href=\"$enum.C.html\">$enum.C</a>";
        s/$newstr\.C/$newerstr/g; 
      }

    }
  }

  print HLPF "$_\n";		# output line to html file
}
}

# ------------------------------------------------------------------------------

