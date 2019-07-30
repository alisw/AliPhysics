#!/bin/bash
#
# Shell script to compare content of the OCDB entries.
# Usage:
# 1) source functios 
# source $ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh

# ocdbMakeTable() 
#       Usage: bash $inputFile $flag $outputFile
# dumpObject()
#       Usage: bash $inputFile $object_name $dump_type [XML/MI] $outfile
# diffObject
#       Usage: bash $inputFile1 $inputFile2 $object_name $dump_type [XML/MI] $outfile




# 1.) Dump software version including git verson:
#     (source $ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh; ocdbMapInfo > software.list)
#
# 2.) Dump ocdb table what was used in the reconstruction:
#     (source $ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh;   ocdbMakeTable AliESDs_Barrel.root "esd" OCDBrec.list )
#     (source $ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh;   ocdbMakeTable galice.root MC OCDBsim.list)
#
# 3.a) Dump the content of the OCDB entry in human readable format as an XML:
#      (source $ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh;  dumpObject  $ALICE_ROOT/OCDB/TPC/Calib/PadNoise/Run0_999999999_v1_s0.root  "object" "XML"  TPC_Calib_PadNoise_Run0_999999999_v1_s0.xml ); 
#
# 3b.) Dump a summary of the OCDB entry in human readable format as provided by object obj->Print resp. ob->Dump(): 
#      (source $ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh;  dumpObject  $ALICE_ROOT/OCDB/TPC/Calib/PadNoise/Run0_999999999_v1_s0.root  "object" "pocdb0"  TPC_Calib_PadNoise_Run0_999999999_v1_s0.print ); 
#      (source $ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh;  dumpObject  $ALICE_ROOT/OCDB/TPC/Calib/PadNoise/Run0_999999999_v1_s0.root  "object" "docdb"  TPC_Calib_PadNoise_Run0_999999999_v1_s0.dump ); 


# Origin marian.ivanov@cern.ch,  j.wagner@cern.ch
 
if [ "$1" == "-h" ]; then
  echo Usage: 
  echo '(source$ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh <action> <par0> <par1> ... '
  echo "==============================="
  echo Example usage ocdbMakeTable
  echo '(source $ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh; ocdbMakeTable $esdFile   "esd" OCDBrec.list )'
  echo '(source $ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh; ocdbMakeTable $mcGalice  MC $outputDir/OCDBsim.list )'
  #
  echo "==============================="
  echo Example usage ocdbDiffTable
  echo '(source $ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh; ocdbDiffTable  $outputDirMC/OCDBrec.list $outputDir/OCDBrec.list TPC 2)'
fi

AliOCDBtoolkit_INIT=1
[ -z "$ALILOG_HOST" ] && source $ALICE_ROOT/libexec/alilog4bash.sh


dumpOCDBXML(){
   [ -z $1 ] && cat <<HELP_USAGE
   dumpOCDBXML - Dump OCDB file as an xml file
   # * param1: Input OCDB snapshot
   # * param2: Input object
   # * param3: Output XML path
   # Example usage:
   #  dump object TPC*Calib*RecoParam  from OCDB snapshot
   #  dumpOCDBXML /lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/ALIROOT-7077/test1211_PbPb_MCTail_2_2/OCDB/OCDBrec.root TPC*Calib*RecoParam TPC_Calib_RecoParam.xml
HELP_USAGE
   [ -z $1 ] && return
    alilog_info "dumpOCDBXML Begin $1 $2 $3"
    lambdaOCDB $1 $3  $2 "" >>/dev/null
    alilog_info "dumpOCDBXML End  $1 $2 $3"
}

lambdaOCDB(){
    [ -z $1 ] && cat <<HELP_USAGE
# lambdaOCDB - modify OCDB entry with user defined code segment
# modification macro lambdaOCDB.C created
# Usage:
# lambdaOCD  param1 [param2]
# * param1: Input OCDB path
# * param2: Output file name
# * param3: object name
# * param4: Object lambda

#  * Example 1: print content of the Reco param and save object to another file
    lambdaOCDB /cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB/TPC/Calib/RecoParam/Run0_244339_v18_s0.root Run0_244339_v18_s0.xml AliCDBEntry  o.Print\(\"all\"\)\;
#  * Example 2: create lambda code and execute it:
    lambdaCode='AliTPCRecoParam*p=o;p->Dump();'
    lambdaOCDB /cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB/TPC/Calib/RecoParam/Run0_244339_v18_s0.root Run0_244339_v18_s0.xml  AliCDBEntry \$lambdaCode;
#  * Example 4: create lambda code to change ion tail
    lambdaCode='TObjArray*p=(TObjArray*)o;for(i=0;i<4;i++)((AliTPCRecoParam*)p->At(i))->SetCrosstalkCorrection(1.5);'
    lambdaOCDB /cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB/TPC/Calib/RecoParam/Run0_244339_v18_s0.root Run0_244339_v18_s0.xml AliCDBEntry  \$lambdaCode;
HELP_USAGE
[ -z $1 ] && return

cat  << LambdaOCDB   > lambdaOCDB.C
    {
        TFile *fin = TFile::Open("$1");
        AliCDBEntry* entry=(AliCDBEntry*)fin->Get("$3");
        TObject * o=(TObject*)entry->GetObject();
        $4
        TFile *fout= TFile::Open("$2","recreate");
        entry->Write("AliCDBEntry");
        delete fout;
    }
LambdaOCDB
      alilog_info "lambdaOCD  $1 $2 $3 $4"
    root -b -q lambdaOCDB.C
}



ocdbMakeTable(){
#
# create a text file with the OCDB setupt descriptors#
#
# Input: 
#   $1 inputFile name
#   $2 flag type of the input file
#   flags: 
#   log file = "log"
#   AliESDs.root = "ESD"
#   galice.root = "MC"
# Output:
#   $3 output file name
    if [ $# -lt 3 ] ; then
        echo "Usage: $0 \$inputFile \$flag \$outputFile"
        return 1
    fi
    export ALIROOT_FORCE_COREDUMP=1
    local inFile=${1}
    local inFlag=${2}
    local outFile=${3}
    shift 3
    #if [ ! -f ${inFile} ] ; then 
    #    echo ${inFile} not found!
    #    return 1
    #fi
    #if [ -f ${outFile} ] ; then 
    #    >${outFile}
    #fi

    tmpscript=$(mktemp)
    cat > ${tmpscript} <<HEREDOC
        {
            AliOCDBtoolkit::DumpOCDBAsTxt("${inFile}","${inFlag}","${outFile}");
        }
HEREDOC

    aliroot -l -q -b ${tmpscript} 
    sleep 60 && rm ${tmpscript} &
    echo Filtering alien cache part 
    cat ${outFile} | sed -e s_"alien://?User=?DBFolder="_"alien://"_g -e s_"?SE=default?CacheFolder=?OperateDisconnected=1?CacheSize=1073741824?CleanupInterval=0"__g > ${outFile}.tmp
    cp  ${outFile}.tmp ${outFile}
    return 1
}


dumpObject(){
#
#
#  Input:
#    $1 path
#    $2 obj name 
#    $3 type of the dump (Print, Dump, XML or MI recursive dump )
#  Output:
#    $4 output file name   
    export ALIROOT_FORCE_COREDUMP=1
    declare -a dumpOptions=("docdb" "xml" "pocdb" "mi")  # supported options. options are converted into lower case
    if [ $# -lt 4 ] ; then
        echo "Usage: $0 \$inputFile \$object_name \$dump_type [POCDB/DOCDB/XML/MI] \$outfile"
        return 1
    fi
    local inFile=${1}
    local fobject=${2}
    local ftype=${3}
    local outFile=${4}
    shift 4
    # check if type is in list supported dump options
    ftypeLower=`echo $ftype | tr '[:upper:]' '[:lower:]'`
    arraylength=${#dumpOptions[@]};
    isOK=0;    
    for (( i=1; i<${arraylength}+1; i++ )); 
      do if [[ $ftypeLower == *"${dumpOptions[$i-1]}"* ]];  then 
	  isOK=1; 
	  echo match ${dumpOptions[$i-1]};  
      fi;       
    done;
    # exit if not supported option
    if [ $isOK == 0 ]; then
	echo Not supported format $ftype;
	echo Please Use;
	for (( i=1; i<${arraylength}+1; i++ )); do echo $i " / " ${arraylength} " : " ${dumpOptions[$i-1]}; done;
        return 1;
    fi;

    tmpscript=$(mktemp)
    cat > ${tmpscript} <<HEREDOC
        {
            AliOCDBtoolkit::DumpOCDBFile("${inFile}","${outFile}",1,"${ftype}");
        }
HEREDOC

    aliroot -l -q -b ${tmpscript} 
    sleep 60 && rm ${tmpscript} &
    return 1
}

diffObject(){
#
#
#  Input:
#    $1 path0
#    $2 path1
#    $3 obj name 
#    $4 type of the dump (xml or MI recursive dump )
#  Output:
#    $5 output diff file name   
    export ALIROOT_FORCE_COREDUMP=1
    if [ $# -lt 5 ] ; then
        echo "Usage: $0 \$inputFile1 \$inputFile2 \$object_name \$dump_type [XML/MI] \$outfile"
        return 1
    fi
    local inFile1=${1}
    local inFile2=${2}
    local fobject=${3}
    local ftype=${4}
    local outFile=${5}
    shift 5
    local tmp1=$(mktemp)
    local tmp2=$(mktemp)
    if [ ${ftype} = "XML" ] ; then
        isXML=kTRUE
        tmp1="${tmp1}.xml"
        tmp2="${tmp2}.xml"
    elif [ ${ftype} = "MI" ] ; then
        isXML=kFALSE
    else
        echo "option ${ftype} not recognized! Use \"XML\" or \"MI\"!"
        return 1
    fi
    dumpObject ${inFile1} ${fobject} ${ftype} ${tmp1%.xml}
    dumpObject ${inFile2} ${fobject} ${ftype} ${tmp2%.xml}
    diff ${tmp1} ${tmp2} >${outFile}
    rm ${tmp1} ${tmp2} 2>/dev/null
    rm "${tmp1}.xml" "${tmp2}.xml" 2>/dev/null
    return 1
}

dumpOCDBDiffTable(){
#
# Dump  differences between the OCDB tables -comparison based on the hash values of OCDB entries
# Input:
#   $1  - list 1
#   $2  - list 2
# Output:
#   difference is stdout
    export ALIROOT_FORCE_COREDUMP=1
    list1=$1
    list2=$2
    shift 2
    cdbName=$(cut -f1 $list1)
    for i in $cdbName ; do
        line1="$(grep $i $list1)"
        line2="$(grep $i $list2)"
        if [ -z "${line2}" ] ; then
            echo $i not found in $list2!
            continue
        fi
        match1=$(echo $line1 | cut -d' ' -f 5)
        match2=$(echo $line2 | cut -d' ' -f 5)
        if [ "$match1" -ne "$match2" ] ; then
            echo $i doesnt match:
            echo $(echo $line1| awk '{print $2 "/" $3}') \#hash: $match1
            echo $(echo $line2| awk '{print $2 "/" $3}') \#hash: $match2
            echo 
        fi
    done
}

ocdbMapInfo(){
    #
    # Dump basic AliRoot session information
    # Tag values which are considered to be writtent to the OCDB metadata
    # (source $ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh; ocdbMapInfo )
    echo DATE":"`date`
    echo USER":"$USER
    echo HOSTNAME":"$HOSTNAME
    echo WorkingDirectory":"`pwd`
    echo ROOTVERSION":"`root-config --version`
    echo ALICE_ROOT":"$(git -C $ALICE_ROOT/../src/ describe)
    echo ALICE_PHYSICS":"$(git -C $ALICE_PHYSICS/../src/ describe)
}


#
# Example usage+ developer test routines. 
#

example1(){
    ocdbMakeTable "/hera/alice/jwagner/simulations/benchmark/aliroot_tests/ppbench/rec.log" "log" "testout"
}
example2(){
    dumpObject "/hera/alice/jwagner/OCDB/temp/TPC/Calib/RecoParam/Run0_999999999_v0_s0.root" "object" "XML" "testout2XML"
    dumpObject "/hera/alice/jwagner/OCDB/temp/TPC/Calib/RecoParam/Run0_999999999_v0_s0.root" "object" "MI" "testout2MI" 
}
example3(){
    file1="/hera/alice/jwagner/OCDB/temp/TPC/Calib/RecoParam/Run0_999999999_v0_s0.root"
    file2="$ALICE_ROOT/OCDB/TPC/Calib/RecoParam/Run0_999999999_v0_s0.root"
    diffObject ${file1} ${file2} "object" "MI" "testdiffMI"
    diffObject ${file1} ${file2} "object" "XML" "testdiffXML"
}

developerTest(){
    source /hera/alice/jwagner/software/aliroot/loadMyAliroot.sh TPCdev
    example1
    example2
    example3
}


diffConfig(){
    #
    # diff configuaration files ignoring trivial differences between the OCDBprefixes 
    #
    file1=$1
    file2=$2

    cat $file1 | sed s_"alien://folder="_"ocdbprefix"_g | sed s_"alien://Folder="_"ocdbprefix"_g | sed s_"local://\${ALICE\_OCDB}"_"ocdbprefix"_g  >${file1}_ocdbstripped
    cat $file2 | sed s_"alien://folder="_"ocdbprefix"_g | sed s_"alien://Folder="_"ocdbprefix"_g | sed s_"local://\${ALICE\_OCDB}"_"ocdbprefix"_g   >${file2}_ocdbstripped
    diff  ${file1}_ocdbstripped  ${file2}_ocdbstripped  > ${file1}_ocdbstrippeddiff

} 


ocdbDiffTable(){
    #
    # Make a diff between the 2 OCDB tables
    # example usage : (source $ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh; ocdbDiffTable  $outputDirMC/OCDBrec.list $outputDir/OCDBrec.list TPC 2)
    # 
    file1=$1
    file2=$2
    mask=$3
    filterType=$4
    nparam=$#
    if [ $nparam -lt 4 ] ; then
	echo "Error - The total length of all arguments is: $nparam: ";
	echo "Oputput dumped to the stdout"
        echo example usage:   "(source $ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh; ocdbDiffTable $file1 $file2 $mask $filterType; ) > mydiff.txt"
        echo example usage:   "(source $ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh; ocdbDiffTable $file1 $file2 $mask 1 ) > ocdbDiffObject.JIRA"
        echo example usage:   "(source $ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh; ocdbDiffTable $file1 $file2 $mask 2 ) > ocdbDiffTable.JIRA"
        echo example usage:   "(source $ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh; ocdbDiffTable $file1 $file2 $mask 4 ) > ocdbDiffTable.html"
        return 0;
    fi;
    # JIRA diffs for 2 XML OCDB object dumps
    diff  -W 1000 -y  --suppress-common-lines $file1  $file2 | grep $mask> ocdbDiffTable.tmp
    if [ $filterType -eq 1 ] ; then
        #  diff y already properly formatted - add | at the beginning and at the end of the line
	cat ocdbDiffTable.tmp  | sed   -e "s/\([\t ,.]\)\1*/\ /g" -e "s_^_|_g"  -e 's_$_|_g'
    fi;
    # JIRA diffs for the OCDB table dump
    if [ $filterType -eq 2 ] ; then	
        # header
	echo "|Entry|Path0|Path1|Entry0|Entry1|Size0|Size1|Hash0|Hash1|"
        # different content print
	cat ocdbDiffTable.tmp | sed 's_|__g'  | grep -v "<" | grep -v ">" |
	while read CMD; do echo $CMD | gawk '{print "|"$1"|"$2"|"$7"|"$3"|"$8"|"$4"|"$9"|"$5"|"$10"|" }'; done
        # mising content print
	cat ocdbDiffTable.tmp | sed 's_|__g'  | grep -e  ">" |
	while read CMD; do echo $CMD | gawk '{print "|"$2"|-|"$3"|-|"$4"|-|"$5"|-|"$6"|" }'; done
	cat ocdbDiffTable.tmp | sed 's_|__g'  | grep -e  "<" |
	while read CMD; do echo $CMD | gawk '{print "|"$1"|-|"$2"|-|"$3"|-|"$4"|-|"$5"|" }'; done
    fi;
    # HTML  diffs for 2 XML OCDB object dumps NOT YET 
    if [ $filterType -eq 3 ] ; then
        #  diff y already properly formatted - add | at the beginning and at the end of the line
	cat ocdbDiffTable.tmp  | sed   -e "s/\([\t ,.]\)\1*/\ /g" -e "s_^_|_g"  -e 's_$_|_g'
    fi;
    # HTML diffs for the OCDB table dump
    if [ $filterType -eq 4 ] ; then
        # header
	echo '<table style="width:100%">'
	echo "<tr><th>Entry</th><th>Path0</th><th>Path1</th><th>Entry0</th><th>Entry1</th><th>Size0</th><th>Size1</th><th>Hash0</th><th>Hash1</th><th></tr>"
        # different content print
	cat ocdbDiffTable.tmp | sed 's_|__g'  | grep -v "<" | grep -v ">" |
	while read CMD; do echo $CMD | gawk '{print "<tr><td>"$1"</td><td>"$2"</td><td>"$7"</td><td>"$3"</td><td>"$8"</td><td>"$4"</td><td>"$9"</td><td>"$5"</td><td>"$10"</td></tr>" }'; done
        # mising content print
	cat ocdbDiffTable.tmp | sed 's_|__g'  | grep -e  ">" |
	while read CMD; do echo $CMD | gawk '{print "<tr><td>"$2"</td><td>-</td><td>"$3"</td><td>-</td><td>"$4"</td><td>-</td><td>"$5"</td><td>-</td><td>"$6"</td></tr>" }'; done
	cat ocdbDiffTable.tmp | sed 's_|__g'  | grep -e  "<" |
	while read CMD; do echo $CMD | gawk '{print "<tr><td>"$1"</td><td>-</td><td>"$2"</td><td>-</td><td>"$3"</td><td>-</td><td>"$4"</td><td>-</td><td>"$5"</td></tr>" }'; done
	echo '</table>'
   fi;
   rm ocdbDiffTable.tmp
}