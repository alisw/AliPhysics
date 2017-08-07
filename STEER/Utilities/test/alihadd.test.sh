#!/usr/bin/env bash
# ( source $AliRoot_SRC/STEER/Utilities/test/alihadd.test.sh; configMI; testBadChunks5   2>&1 | tee alihadd.testBadChunks5.log)
# ( source $AliRoot_SRC/STEER/Utilities/test/alihadd.test.sh; configMI; testBatchProduction| tee alihaddbatch.test.log)
configMI(){
    source $ALICE_ROOT/libexec/alihadd.sh
    export prefixList=$NOTES/JIRA/PWGPP-340/testlists/
    batchCommand="/usr/bin/sbatch"
    batchFlags="--get-user-env --mem-per-cpu=6096 --time=8:00:00"
}

testBadChunks5(){
    #
    #   test merging of the files with partially corrupted input
    #       test switch -checkKeys
    #       test switch -s
    #       test prefetch not yet implemented
    #       Input test sample consist of  list selected files with 5
    alilog_info testBadChunks5.Begin
    wdir=`pwd`;
    mkdir -p $wdir/testBadChunks5/
    cd $wdir/testBadChunks5/
    cp $prefixList/testBadChunks5.list input.list
    alihaddsh -k -s 200000000   -checkKeys highPt -prefetch 1 -v 1 -timeOut 30 test.root @input.list  2>&1 |tee testBadChunks5.log
    # check output - valid only for particular input data
    status=0
    cd $wdir/testBadChunks5/
    [ `wc -l <input.list.Good` != 5 ] && alilog_error "alihadd.test.testBadChunks5: Invalid number of good files input.list.Good!=5"
    [ `wc -l <input.list.Bad` != 5 ] && alilog_error "alihadd.test.testBadChunks5: Invalid number of bad files input.list.Bad!=5"
    if [ -z  "`cat $wdir/testBadChunks5/testBadChunks5.log | grep -e "SUCCESS.*alihadd"`" ]; then
        alilog_error "alihadd.test.testBadChunks5: alihadd merging failed";
        return 1;
    fi
    alilog_success alihadd.test.testBadChunks5
    #clean
    #rm *root
    alilog_info testBadChunks5.End
}

testBatchProduction(){
    # submit test jobs to test bigger statistic
    #     makeflow to be used here in future
    alilog_info testBatchProduction.Begin
    wdir=`pwd`
    mkdir -p $wdir/testBatchProduction
    cd $wdir/testBatchProduction
    cp $prefixList/testBatchProduction.list .
    split testBatchProduction.list -d dir -l 25 --additional-suffix .list
    for a in `ls dir*list`; do
        dname=`echo $a |sed s_.list__`
        mkdir -p  $dname
        mv $a $dname/input.list
        cd $dname
        echo '#!/usr/bin/env bash' >command.sh
        echo "(source $ALICE_ROOT/libexec/alihadd.sh; alihaddsh -k -s 2000000000   - highPt -prefetch 1 -v 1 test.root @input.list | tee testBatchProduction.log)" >> command.sh
        chmod a+x ./command.sh
        alilog_info "testBatchProduction.Submit `pwd`  $batchCommand $batchFlags command.sh"
        $batchCommand $batchFlags command.sh
        cd ..
    done
}

testMakeFlow(){
  # not working yet
  alilog_info testMakeFlow.Begin
  wdir=`pwd`
  mkdir -p $wdir/testMakeFlow
  cd $wdir/testMakeFlow
  echo '#!/usr/bin/env bash' >command.sh
  echo "(source $ALICE_ROOT/libexec/alihadd.sh; alihaddsh -k -s 2000000000   -checkKeys highPt -prefetch 1 -v 1 test.root @input.list | tee testBatchProduction.log)" >> command.sh
  chmod a+x ./command.sh
  #
  tail -n 50  $prefixList/testBatchProduction.list > testMakeFlow.list
  split testMakeFlow.list -d dir -l 10 --additional-suffix .list
  echo >  hadd.makeflow
  for a in `ls dir*list`; do
    dname=`echo $a |sed s_.list__`
    mkdir -p  $dname
    mv $a $dname/input.list
    cd $dname
    echo "$dname/testBatchProduction.log: $dname/input.list" >>    ../hadd.makeflow
    echo "(cd $wdir/$dname; touch $dname/testBatchProduction.log)"  >>    ../hadd.makeflow
    cd ..
  done

}

testall(){
    testBadChunks5;
    #[ -n $batchCommand ] && testBatchProduction
}