#!/usr/local/bin/perl
# $Id$
# Ivana Hrivnacova 15.2.00
#
# This script defines the class categories
# and calls doc_alcategory.pl to generate 
# the html category pages


print "Generating html pages for class categories\n"; 

# main categories (packages)
@DIRLIST = "TGeant4";
@DIRLIST = (@DIRLIST,"AliGeant4");
@NAMELIST = "TGeant4";
@NAMELIST = (@NAMELIST,"AliGeant4");

# subcategories
@CATLIST = "global";
@CATLIST = (@CATLIST,"geometry");
@CATLIST = (@CATLIST,"physics");
@CATLIST = (@CATLIST,"event");
@CATLIST = (@CATLIST,"run");
@CATLIST = (@CATLIST,"visualization");

# categories
for( $i = 0 ; $i < $#DIRLIST+1 ; $i++ ) {
#foreach $DIR (@DIRLIST) {
  $DIR = @DIRLIST[$i];
  $NAME = @NAMELIST[$i];
  $DIRPATH = $ENV{'AG4_INSTALL'} . "/../" . $DIR;
  chdir $DIRPATH;
  foreach $CAT (@CATLIST) {
    #chop $CAT;
    # exclude other subdirectories

    print "Processing category: " . $CAT . " of " . $NAME . "\n"; 
    $CATDIRPATH = $DIRPATH . "/" . $CAT;
    $CATNAME = $NAME . "_" . $CAT;      
 
    # generate the category pages
    system $ENV{'AG4_INSTALL'} . "/config/doc_alcategory.pl " . $CAT . " " . $CATNAME;
  }  
}      
