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
@CATLIST = (@CATLIST,"physics");
@CATLIST = (@CATLIST,"event");
@CATLIST = (@CATLIST,"run");
@CATLIST = (@CATLIST,"visualization");
@CATLIST = (@CATLIST,"interfaces");

# create source dir structure
for( $i = 0 ; $i < $#DIRLIST+1 ; $i++ ) {
  $DIR = @DIRLIST[$i];
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

# link main history file   
$HISTORYPATH = $ENV{'AG4_INSTALL'} . "/doc/history";
$HISTORYBASE = $HISTORYPATH . "/History";
chdir $ENV{'AG4_INSTALL'};
`ln -s $HISTORYBASE "History" `;
  
# link source files and history files
for ( $i = 0 ; $i < $#DIRLIST+1 ; $i++ ) {
  $DIR = @DIRLIST[$i];
  $NAME = @NAMELIST[$i];
  $DIRPATH = $ENV{'AG4_INSTALL'} . "/../" . $DIR;
  $RELDIRPATH = "../../..";

  # History categories files  
  $HISTORYBASE = $HISTORYPATH . "/" . $NAME . "_History";
  chdir $DIRPATH . "/source";
  `ln -s $HISTORYBASE "History" `;

  foreach $CAT (@CATLIST) {
    chdir $DIRPATH;

    $CATSTRING = "\"Category: " . $CAT . "\"";
    @FILELIST_H   = `find . -name \"*.h\" -exec grep -l  $CATSTRING  {} \\;`;
    @FILELIST_ICC = `find . -name \"*.icc\" -exec grep -l  $CATSTRING  {} \\;`;
    @FILELIST_CXX = `find . -name \"*.cxx\" -exec grep -l  $CATSTRING  {} \\;`;

    print "Processing category: " . $CAT . "\n"; 
    $CATDIRPATH = "source" . "/" . $CAT;
    $INCLUDEPATH = $CATDIRPATH . "/" . "include";
    $SRCPATH = $CATDIRPATH . "/" . "src";
    $CVSBASE  = $DIRPATH . "/CVS";
    $HISTORYBASE = $HISTORYPATH . "/" . $NAME . "_" . $CAT . "_History";
    
    # History subcategories files
    chdir $CATDIRPATH;
    `ln -s $HISTORYBASE "History" `;
 
    # .h files
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

    # .icc files
    foreach $FILEPATH (@FILELIST_ICC) { 
      @TEMP = split('/',$FILEPATH);
      $FILE = @TEMP[@TEMP - 1];
      chop $FILE;
      print "   Linking file " . $FILE . "\n";
      $FILEBASE = $RELDIRPATH . "/" . $FILE;
      `ln -s $FILEBASE $FILE`;
    }

    # .cxx files
    chdir $DIRPATH . "/" . $SRCPATH;
    `ln -s $CVSBASE "CVS" `;
    foreach $FILEPATH (@FILELIST_CXX) { 
      @TEMP = split('/',$FILEPATH);
      $FILE = @TEMP[@TEMP - 1];
      chop $FILE;
      print "   Linking file " . $FILE . "\n";
      $FILEBASE = $RELDIRPATH . "/" . $FILE;
      `ln -s $FILEBASE $FILE`;
    }
  }
}

