#!/usr/local/bin/perl
# $Id$
# Ivana Hrivnacova 12.2.99
#
# HTML documentation is created for all 
# source code files: *.hh *.cc *.icc 
# generic makefiles: *.gmk 
# and configuration setup scripts

# file extensions
#$INCEXT = ".hh"
#$SRCEXT = ".cc"
#$MKFEXT = ".gmk"

# create doc directory if it does not exist
$CURDIR = `pwd`;
chdir $ENV{'G4INSTALL'};
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
chdir $ENV{'G4INSTALL'};

# create tmpdoc directory is it does not exist
# or clean it
if (! grep(/tmpdoc/, `ls`)) {
  `mkdir tmpdoc` ;
  print "Directory tmpdoc has been created." . "\n";
} else {
  print "Cleaning tmpdoc" . "\n";
  `rm -fr tmpdoc/*`;
}  

# copy everything for documentation to tmpdoc
@FILELIST = `find . -name "*.ddl"`;
@FILELIST = (@FILELIST, `find . -name "*.h"`);
@FILELIST = (@FILELIST, `find . -name "*.hh"`);
@FILELIST = (@FILELIST, `find . -name "*.cc"`);
@FILELIST = (@FILELIST, `find . -name "*.icc"`);
@FILELIST = (@FILELIST, `find config -name "*.gmk"`);
@FILELIST = (@FILELIST, `find config -name "setup*"`);
@FILELIST = (@FILELIST, `find config -name "*boot"`);

print "Copying files to tmpdoc" . "\n";
foreach $FILE (@FILELIST) {
  chop $FILE;
  # exclude dictionary classes
  if (!grep(/Dict/,$FILE)) {
    `cp $FILE tmpdoc`;
  }  
  #print "$FILE has been copied to tmpdoc" . "\n";
}

# mv *.cc to *.C
# what is recquired by ddl2html.pl
chdir tmpdoc;
print "Renaming files" . "\n";
@CXXLIST = `ls *.cc`;
foreach $CXXFILE (@CXXLIST) {
  chop $CXXFILE;
  $CFILE = `echo $CXXFILE | sed s/.cc/.C/g`;
  `mv $CXXFILE $CFILE`;
}

# mv *.hh to *.h 
@HHLIST = `ls *.hh`;
foreach $HHFILE (@HHLIST) {
  chop $HHFILE;
  $HFILE = `echo $HHFILE | sed s/.hh/.h/g`;
  `mv $HHFILE $HFILE`;
}

# execute the modified P.Binko's script
system $ENV{'AG4_INSTALL'} . "/config/doc_g4ddl2html.pl GEANT4 Project";

# move HTML to doc and remove tmpdoc
$DOCDIR = $ENV{'G4INSTALL'} . "/doc";
`mv HTML $DOCDIR`;  
print "Removing tmpdoc" . "\n";
chdir $ENV{'G4INSTALL'};
`rm -fr tmpdoc`;

# generate the category pages
system $ENV{'AG4_INSTALL'} . "/config/doc_g4sorthtml.pl";

chdir $CURDIR;
