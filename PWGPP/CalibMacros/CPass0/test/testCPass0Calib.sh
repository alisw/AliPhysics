#
# This code should be in $ALICE_PHYSICS/../src/PWGPP/CalibMacros/CPass0/testCPass0Calib.sh
# To test:
#   a.) code is not crashing
#   b.) CPU time and VM usage (using AliSysInfo.h)
#   c.) some physic invariants (?)
#
# 0.) Input test sample defintion
#
inputESD="/alice/data/2015/LHC15f/000226495/pass2/15000226495029.919/AliESDs.root"
inputFriend="/alice/data/2015/LHC15f/000226495/pass2/15000226495029.919/AliESDfriends.root"
inputRun=226495
#
# 1.) Dump software info git.log  should be implemented  in some utilities.sh 
#     Current implementation in $ALICE_ROOT/libexec/utilities.sh
#
wdir=`pwd`
( 
    source $ALICE_ROOT/libexec/alilog4bash.sh
    #
    echo ==============================  >$wdir/git.log 
    alilog "ALICE_PHYSICS"  >>$wdir/git.log 
    cd  $ALICE_PHYSICS/../src/PWGPP/CalibMacros/CPass0/;
    alilog "git describe"  >>$wdir/git.log 
    git describe >>$wdir/git.log;  
    alilog "git status" >>$wdir/git.log
    git status --untracked-files=no . >> $wdir/git.log; 
    git diff >$wdir/aliphysics.diff
    #
    echo ==============================  >>$wdir/git.log 
    alilog "ALICE_ROOT"  >>$wdir/git.log 
    alilog "which aliroot: $(which aliroot)"
    cd  $ALICE_ROOT/../src/;
    alilog "git describe"  >>$wdir/git.log 
    git describe >>$wdir/git.log;  
    alilog "git status" >>$wdir/git.log
    git status --untracked-files=no . >> $wdir/git.log; 
    git diff >$wdir/aliroot.diff 
)
#
# 2.)  
#
alien_cp alien://$inputESD .
alien_cp alien://$inputFriend .
cp $ALICE_PHYSICS/../src/PWGPP/CalibMacros/CPass0/*.C .
cp $ALICE_PHYSICS/../src/PWGPP/CalibMacros/CPass0/*.sh .
aliroot -b -q runCalibTrain.C\($inputRun, \"AliESDs.root\", \"raw://\"\)  2>&1 | tee calib.log





