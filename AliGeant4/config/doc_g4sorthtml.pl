#!/usr/local/bin/perl
# $Id$
# Ivana Hrivnacova 15.2.00
#
# This script defines the class categories
# and calls doc_alcategory.pl to generate 
# the html category pages

print "Generating html pages for class categories\n"; 

# main categories list
$G4SOURCE = $ENV{'G4INSTALL'} . "/source";
chdir $G4SOURCE;
@DIRLIST = `ls`;

# categories
foreach $DIR (@DIRLIST) {
  chop $DIR;
  # exclude unwanted subdirectories/files
  if ($DIR ne "GNUmakefile" && $DIR ne "History" && $DIR ne "nohup.out") {  
    # subcategories are considered only in one level
    $DIRPATH = $G4SOURCE . "/" . $DIR;
    chdir $DIRPATH;
    
    # test if subcategories are present
    @INCLUDE = `find . -name "include"`;
    $FIRSTINCLUDE = @INCLUDE[0];
    chop $FIRSTINCLUDE;
    if ($FIRSTINCLUDE eq "./include") {
      @CATLIST = ".";
    } 
    else {      
      @CATLIST = `ls`;    
    }  
   
    # process subcategories
    foreach $CAT (@CATLIST) {
      chop $CAT;
      # exclude other subdirectories
      if ($CAT ne "GNUmakefile" && $CAT ne "History" && $CAT ne "nohup.out") {
        print "Processing category: " . $CAT . " of " . $DIR . "\n"; 
	if ($CAT eq ".") {
          $CATDIRPATH = $DIRPATH;
          $CATNAME = $DIR;      
	}
	else {
          $CATDIRPATH = $DIRPATH . "/" . $CAT;
          $CATNAME = $DIR . "_" . $CAT;      
	}  
 
        # generate the category pages
        system $ENV{'AG4_INSTALL'} . "/config/doc_g4category.pl " . $CATDIRPATH . " " . $CATNAME;
      }	
    }
  }  
}      
