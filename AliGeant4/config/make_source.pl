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

# subcategories
@CATLIST = "global";
@CATLIST = (@CATLIST,"geometry");
@CATLIST = (@CATLIST,"physics");
@CATLIST = (@CATLIST,"event");
@CATLIST = (@CATLIST,"run");
@CATLIST = (@CATLIST,"visualization");

# create source dir structure
for( $i = 0 ; $i < $#DIRLIST+1 ; $i++ ) {
  $DIR = @DIRLIST[$i];
  $NAME = @NAMELIST[$i];
  $DIRPATH = $ENV{'AG4_INSTALL'} . "/../" . $DIR;
  chdir $DIRPATH;
  if (! grep(/source/, `ls`)) {
    `mkdir source ` ;
    foreach $CAT (@CATLIST) {
       $CATDIRPATH = "source" . "/" . $CAT;
       $INCLUDEPATH = $CATDIRPATH . "/" . "include";
       $SRCPATH = $CATDIRPATH . "/" . "src";
       `mkdir $CATDIRPATH` ;
       `mkdir $INCLUDEPATH` ;
       `mkdir $SRCPATH` ;
    }    
    print $DIR . "/source directory has been created." . "\n";
  }
}  

# link source files
for( $i = 0 ; $i < $#DIRLIST+1 ; $i++ ) {
  $DIR = @DIRLIST[$i];
  $NAME = @NAMELIST[$i];
  $DIRPATH = $ENV{'AG4_INSTALL'} . "/../" . $DIR;
  $RELDIRPATH = "../../..";

  foreach $CAT (@CATLIST) {
    chdir $DIRPATH;

    $CATSTRING = "\"Category: " . $CAT . "\"";
    @FILELIST_H   = `find . -maxdepth 1 -name \"*.h\" -exec grep -l  $CATSTRING  {} \\;`;
    @FILELIST_ICC = `find . -maxdepth 1 -name \"*.icc\" -exec grep -l  $CATSTRING  {} \\;`;
    @FILELIST_CXX = `find . -maxdepth 1 -name \"*.cxx\" -exec grep -l  $CATSTRING  {} \\;`;

    print "Processing category: " . $CAT . "\n"; 
    $CATDIRPATH = "source" . "/" . $CAT;
    $INCLUDEPATH = $CATDIRPATH . "/" . "include";
    $SRCPATH = $CATDIRPATH . "/" . "src";
    $CVSBASE  = $DIRPATH . "/CVS";
 
    chdir $DIRPATH . "/" . $INCLUDEPATH;
    `ln -s $CVSBASE "CVS" `;
    foreach $FILEPATH (@FILELIST_H) { 
      @TEMP = split('/',$FILEPATH);
      $FILE = @TEMP[@TEMP - 1];
      chop $FILE;
      print "   Linking file " . $FILE . "\n";
      $FILEBASE = $RELDIRPATH . "/" . $FILE;
      `ln -s $FILEBASE $FILE`;
    }

    foreach $FILEPATH (@FILELIST_ICC) { 
      @TEMP = split('/',$FILEPATH);
      $FILE = @TEMP[@TEMP - 1];
      chop $FILE;
      print "   Linking file " . $FILE . "\n";
      $FILEBASE = $RELDIRPATH . "/" . $FILE;
      `ln -s $FILEBASE $FILE`;
    }

    chdir $DIRPATH . "/" . $SRCPATH;
    `ln -s $CVSBASE "CVS" `;
    foreach $FILEPATH (@FILELIST_CXX) { 
      @TEMP = split('/',$FILEPATH);
      $FILE = @TEMP[@TEMP - 1];
      chop $FILE;
      print "   Linking file " . $FILE . " in " . $SRCPATH . "\n";
      $FILEBASE = $RELDIRPATH . "/" . $FILE;
      `ln -s $FILEBASE $FILE`;
    }
  }
}  
