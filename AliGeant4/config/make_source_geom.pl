#!/usr/local/bin/perl
# $Id$
# by I. Hrivnacova, 22.6. 2000
#
# This script creates source directory structured
# according subcategories with links to the flat
# source structure.
  
# main categories (packages)
@DIRLIST = "TGeant4";
@DIRLIST = (@DIRLIST,"AliGeant4");
@NAMELIST = "g4mc";
@NAMELIST = (@NAMELIST,"alice");

# subcategories
@CATLIST = "global";
@CATLIST = (@CATLIST,"geometry");

# link source files and history files
for ( $i = 0 ; $i < $#DIRLIST+1 ; $i++ ) {
  $DIR = @DIRLIST[$i];
  $NAME = @NAMELIST[$i];
  $DIRPATH = $ENV{'AG4_INSTALL'} . "/../" . $DIR;
  $TARGETPATH = $DIRPATH . "_geometry";

  foreach $CAT (@CATLIST) {
    chdir $DIRPATH;

    $CATSTRING = "\"Category: " . $CAT . "\"";
    @FILELIST_H   = `find . -name \"*.h\" -exec grep -l  $CATSTRING  {} \\;`;
    @FILELIST_ICC = `find . -name \"*.icc\" -exec grep -l  $CATSTRING  {} \\;`;
    @FILELIST_CXX = `find . -name \"*.cxx\" -exec grep -l  $CATSTRING  {} \\;`;

    print "Processing category: " . $CAT . "\n"; 

    # .h files
    chdir $TARGETPATH;    
    foreach $FILEPATH (@FILELIST_H) { 
      @TEMP = split('/',$FILEPATH);
      $FILE = @TEMP[@TEMP - 1];
      chop $FILE;
      print "   Linking file " . $FILE . "\n";
      $FILEPATH = $DIRPATH . "/" . $FILE;
      #print "   Linking file " . $FILEPATH . " " . $FILE . "\n";
      `ln -s $FILEPATH $FILE`;
    }

    # .icc files
    foreach $FILEPATH (@FILELIST_ICC) { 
      @TEMP = split('/',$FILEPATH);
      $FILE = @TEMP[@TEMP - 1];
      chop $FILE;
      print "   Linking file " . $FILE . "\n";
      $FILEPATH = $DIRPATH . "/" . $FILE;
      #print "   Linking file " . $FILEPATH . " " . $FILE . "\n";
      `ln -s $FILEPATH $FILE`;
    }

    # .cxx files
    foreach $FILEPATH (@FILELIST_CXX) { 
      @TEMP = split('/',$FILEPATH);
      $FILE = @TEMP[@TEMP - 1];
      chop $FILE;
      print "   Linking file " . $FILE . "\n";
      $FILEPATH = $DIRPATH . "/" . $FILE;
      #print "   Linking file " . $FILEPATH . " " . $FILEDEST . "\n";
      `ln -s $FILEPATH $FILE`;
    }
  }
}

