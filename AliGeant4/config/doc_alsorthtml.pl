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

# AliRoot categories
@EXTDIRLIST = "STEER";
@EXTDIRLIST = (@EXTDIRLIST,"EVGEN");
@EXTDIRLIST = (@EXTDIRLIST,"TGeant3");
@EXTDIRLIST = (@EXTDIRLIST,"ALIFAST");
@EXTDIRLIST = (@EXTDIRLIST,"CASTOR");
@EXTDIRLIST = (@EXTDIRLIST,"FMD");
@EXTDIRLIST = (@EXTDIRLIST,"ITS");
@EXTDIRLIST = (@EXTDIRLIST,"MUON");
@EXTDIRLIST = (@EXTDIRLIST,"PHOS");
@EXTDIRLIST = (@EXTDIRLIST,"PMD");
@EXTDIRLIST = (@EXTDIRLIST,"RICH");
@EXTDIRLIST = (@EXTDIRLIST,"START");
@EXTDIRLIST = (@EXTDIRLIST,"STRUCT");
@EXTDIRLIST = (@EXTDIRLIST,"TOF");
@EXTDIRLIST = (@EXTDIRLIST,"TPC");
@EXTDIRLIST = (@EXTDIRLIST,"TRD");
@EXTDIRLIST = (@EXTDIRLIST,"ZDC");

# categories
for( $i = 0 ; $i < $#DIRLIST+1 ; $i++ ) {
  $DIR = @DIRLIST[$i];
  $NAME = @NAMELIST[$i];
  $DIRPATH = $ENV{'AG4_INSTALL'} . "/../" . $DIR;
  chdir $DIRPATH;
  foreach $CAT (@CATLIST) {

    $CATDIRPATH = $DIRPATH . "/" . $CAT;
    $CATNAME = $NAME . "_" . $CAT;      
 
    # generate the category pages
    system $ENV{'AG4_INSTALL'} . "/config/doc_alcategory.pl " . $CAT . " " . $CATNAME . " TRUE";
  }  
}      

# AliRoot categories
for( $i = 0 ; $i < $#EXTDIRLIST+1 ; $i++ ) {
#foreach $DIR (@EXTDIRLIST) {
  $DIR = @EXTDIRLIST[$i];
  $DIRPATH = $ENV{'ALICE_ROOT'} . "/" . $DIR;
  chdir $DIRPATH;
 
  # generate the category pages
  system $ENV{'AG4_INSTALL'} . "/config/doc_alcategory.pl " . $DIRPATH . " " . $DIR . " FALSE";
}      
