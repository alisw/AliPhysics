#!/usr/local/bin/perl
# $Id$
# Ivana Hrivnacova 12.2.99
#
# HTML documentation is created for all 
# source code files: *.h *.cxx 
# makefiles: Makefile, *.gmk 
# and configuration setup scripts

# create doc directory if it does not exist
$CURDIR = `pwd`;
chdir $ENV{'AG4_INSTALL'};
if (! grep(/doc/, `ls`)) {
  `mkdir doc` ;
  print "Directory doc has been created." . "\n";
};
# move doc/HTML directory to doc/HTML.old
chdir doc;
if (grep(/HTML.old/, `ls`)) {
  print "Cleaning HTML.old" . "\n";
  `rm -fr HTML.old`;
}
if (grep(/HTML/, `ls`)) {
  `mkdir HTML.old`;
  `mv HTML HTML.old`;
  print "Old HTML directory has been saved." . "\n";
}
chdir $ENV{'AG4_INSTALL'};

# create tmpdoc directory is it does not exist
# or clean it
if (! grep(/tmpdoc/, `ls`)) {
  `mkdir tmpdoc` ;
  print "Directory tmpdoc has been created." . "\n";
} else {
  print "Cleaning tmpdoc" . "\n";
  `rm -fr tmpdoc/*`;
}  

# select directory that will be processed
$SOURCEDIR = ". ../TGeant4 ../STEER ../EVGEN ../TGeant3 ../ALIFAST ../ALIROOT ";
$SOURCEDIR = $SOURCEDIR . ". ../CONTAINERS ../THijing ";
$SOURCEDIR = $SOURCEDIR . "../CASTOR ../FMD ../ITS ../MUON ../PHOS ../PMD ";
$SOURCEDIR = $SOURCEDIR . "../RICH ../START ../STRUCT ../TOF ../TPC ../TRD ../ZDC";

# copy everything for documentation to tmpdoc
@FILELIST = `find $SOURCEDIR -name "*.ddl"`;
@FILELIST = (@FILELIST, `find $SOURCEDIR -name "*.h"`);
@FILELIST = (@FILELIST, `find $SOURCEDIR -name "*.cxx"`);
@FILELIST = (@FILELIST, `find $SOURCEDIR -name "*.icc"`);
@FILELIST = (@FILELIST, `find $SOURCEDIR -name "Makefile"`);
@FILELIST = (@FILELIST, `find $SOURCEDIR -name "*.gmk"`);
@FILELIST = (@FILELIST, `find $SOURCEDIR -name "setup*"`);

print "Copying files to tmpdoc" . "\n";
foreach $FILE (@FILELIST) {
  chop $FILE;
  # exclude dictionary classes
  if (!grep(/Dict/,$FILE)) {
    `cp $FILE tmpdoc`;
  }  
  #print "$FILE has been copied to tmpdoc" . "\n";
}

# mv *.cxx to *.C 
# what is recquired by ddl2html.pl
chdir tmpdoc;
print "Renaming files" . "\n";
@CXXLIST = `ls *.cxx`;
foreach $CXXFILE (@CXXLIST) {
  chop $CXXFILE;
  $CFILE = `echo $CXXFILE | sed s/.cxx/.C/g`;
  `mv $CXXFILE $CFILE`;
}

# execute the modified P.Binko's script
system $ENV{'AG4_INSTALL'} . "/config/doc_alddl2html.pl ALICE G4 Project";

# move HTML to doc and remove tmpdoc
$DOCDIR = $ENV{'AG4_INSTALL'} . "/doc";
`mv HTML $DOCDIR`;  
print "Removing tmpdoc" . "\n";
chdir $ENV{'AG4_INSTALL'};
`rm -fr tmpdoc`;

# generate the category pages
system $ENV{'AG4_INSTALL'} . "/config/doc_alsorthtml.pl";

chdir $CURDIR;
