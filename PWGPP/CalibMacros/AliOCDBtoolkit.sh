#!/bin/bash
# Shell script to compare content of the OCDB entries.
# Usage:
#    source $ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh
#   ( source $ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh -h)
#         * to get help and list of function
#    enter function whithout parameters to get a help. e.g:
#          ( source $ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh;  dumpObject)
# Origin marian.ivanov@cern.ch,  j.wagner@cern.ch

alias helpCat=cat
[ -x "$(command -v pygmentize)" ] && alias helpCat="pygmentize -O style=borland,linenos=1 -l bash"


if [ "$1" == "-h" ]; then
  helpCat<<HELP_USAGE
  # Usage:
        (source \${ALICE_PHYSICS}/PWGPP/CalibMacros/AliOCDBtoolkit.sh <action> <par0> <par1> ...
  # List of functions:
        dumpOCDBXML       -  NEW help format
        lambdaOCDB        -  NEW help format
        dumpObject        -  NEW help format
        rootFileConvert   -  ROOT->XML and back - does not work for template classes
        ocdbMakeTable     -  NEW help format    - hash does not work  TODO
        cdbMapInfo        -  New help format
        diffObject        -  New help format
        dumpOCDBDiffTable -  TODO
        diffConfig        -  TODO
        ocdbDiffTable     -  TODO

  #    enter function whithout parameters to get a help. e.g:
            ( source \$ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh;  dumpObject)
            ( source \$ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh;  dumpOCDBXML)
            ( source \$ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh;  lambdaOCDB)
            ( source \$ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh;  ocdbMakeTable)
HELP_USAGE
  #===============================
  # Example usage ocdbMakeTable
  #      (source ${ALICE_PHYSICS}/PWGPP/CalibMacros/AliOCDBtoolkit.sh; ocdbMakeTable ${esdFile}   "esd" OCDBrec.list )
  #      (source ${ALICE_PHYSICS}/PWGPP/CalibMacros/AliOCDBtoolkit.sh; ocdbMakeTable ${mcGalice}  MC ${outputDir}/OCDBsim.list )
  #===============================
  #Example usage ocdbDiffTable
  #      (source ${ALICE_PHYSICS}/PWGPP/CalibMacros/AliOCDBtoolf<kit.sh; ocdbDiffTable  ${outputDirMC}/OCDBrec.list ${outputDir}/OCDBrec.list TPC 2)
fi

AliOCDBtoolkit_INIT=1
[ -z "${ALILOG_HOST}" ] && source ${ALICE_ROOT}/libexec/alilog4bash.sh


dumpOCDBXML(){
   [[ $# -ne 3 ]] && helpCat <<HELP_USAGE
   # dumpOCDBXML - Dump OCDB file as an xml file
   # Input parameters
        * param1: Input OCDB snapshot
        * param2: Input object
        * param3: Output XML path
   # Example usage:
        # dump object TPC*Calib*RecoParam  from OCDB snapshot
        dumpOCDBXML /lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/ALIROOT-7077/test1211_PbPb_MCTail_2_2/OCDB/OCDBrec.root TPC*Calib*RecoParam TPC_Calib_RecoParam.xml
HELP_USAGE
    [[ $# -ne 3 ]] && return
    alilog_info "dumpOCDBXML Begin $1 $2 $3"
    lambdaOCDB $1 $3  $2 "" >>/dev/null
    alilog_success "dumpOCDBXML End  $1 $2 $3"
}

rootFileConvert(){
 [[ $# -ne 2 ]] && helpCat <<HELP_USAGE
   # rootFileConvert from ROOT<--> XML
   # NOT WORKING FOR SOME DATAT because of the BUG in TXMLParaser - problem with template classes
                        #  <TMatrixTBase<float> version="5">
   # Input parameters
        * param1: Input  file
        * param2: Output  file
   # Example usage:
        rootFileConvert /cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB/TPC/Calib/RecoParam/Run0_244339_v18_s0.root Run0_244339_v18_s0.xml
        # back conversion does not work because of the problem in TXMLParser for template classes
        rootFileConvert Run0_244339_v18_s0.xml  Run0_244339_v18_s0.root
        #
        rootFileConvert /cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB/TPC/Calib/Parameters/Run0_999999999_v1_s0.root Run0_999999999_v1_s0.xml
        rootFileConvert Run0_999999999_v1_s0.xml Run0_999999999_v1_s1.root
        rootFileConvert Run0_999999999_v1_s1.root Run0_999999999_v1_s1.xml


HELP_USAGE
    [[ $# -ne 2 ]] && return

    tmpscript=$(mktemp)
    inputFile=$1
    outputFile=$2
     tmpscript="test.C"
    #   inputFile=/cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB/TPC/Calib/RecoParam/Run0_244339_v18_s0.root
    #   outputFile=Run0_999999999_v0_s2.xml
    cat > ${tmpscript} <<HEREDOC
        {
            TFile *fin = TFile::Open("${inputFile}");
            TList * keys= fin->GetListOfKeys();
            TFile *fout = TFile::Open("${outputFile}","recreate");
            for (Int_t i=0; i<keys->GetEntries(); i++) (fin->Get(keys->At(i)->GetName()))->Write(keys->At(i)->GetName());
            fout->Close();
        }
HEREDOC
    root.exe -l -b -q ${tmpscript}
}

lambdaOCDB(){
    [[ $# -ne 4 ]]  && helpCat <<HELP_USAGE
    # lambdaOCDB - modify OCDB entry with user defined code segment
    # modification macro lambdaOCDB.C created
    # Usage:
        lambdaOCD  param1 param2 param3 param4
        * param1: Input OCDB path
        * param2: Output file name
        * param3: object name
        * param4: Object lambda

    # Example 1: print content of the Reco param and save object to another file
        lambdaOCDB /cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB/TPC/Calib/RecoParam/Run0_244339_v18_s0.root Run0_244339_v18_s0.xml AliCDBEntry  o.Print\(\"all\"\)\;
    # Example 2: create lambda code and execute it:
        lambdaCode='AliTPCRecoParam*p=o;p->Dump();'
        lambdaOCDB /cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB/TPC/Calib/RecoParam/Run0_244339_v18_s0.root Run0_244339_v18_s0.xml  AliCDBEntry \$lambdaCode;
    # Example 4: create lambda code to change ion tail
        lambdaCode='TObjArray*p=(TObjArray*)o;for(i=0;i<4;i++)((AliTPCRecoParam*)p->At(i))->SetCrosstalkCorrection(1.5);'
        lambdaOCDB /cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB/TPC/Calib/RecoParam/Run0_244339_v18_s0.root Run0_244339_v18_s0.xml AliCDBEntry  \$lambdaCode;
HELP_USAGE
    [[ $# -ne 4 ]] && return
cat  << LambdaOCDB   > lambdaOCDB.C
    {
        TGrid::Connect("alien");
        TFile *fin = TFile::Open("$1");
        AliCDBEntry* entry=(AliCDBEntry*)fin->Get("$3");
        TObject * o=(TObject*)entry->GetObject();
        $4
        TFile *fout= TFile::Open("$2","recreate");
        entry->Write("AliCDBEntry");
        delete fout;
    }
LambdaOCDB
    alilog_info "lambdaOCD BEGIN $1 $2 $3 $4"
    root.exe -l -b -q lambdaOCDB.C
    alilog_success "lambdaOCD  $1 $2 $3"
}



ocdbMakeTable(){
  [[ $# -ne 3 ]]  && helpCat <<HELP_USAGE
    # ocdbMakeTable
    # Usage:
        * param1: Input file
        * param2: Input file type
            *   log file = "log"
            *   AliESDs.root = "ESD"
            *   galice.root = "MC"
        * param3:  output table name
    # Example usage:
        ocdbMakeTable alien:///alice/sim/2020/LHC20c1a/297479/001/AliESDs.root "ESD" LHC20c1a.ESD.table
        ocdbMakeTable alien:///alice/sim/2020/LHC20c1a/297479/001/galice.root  "MC" LHC20c1a.MC.table

HELP_USAGE
    [[ $# -ne 3 ]]  && return
    alilog_info "ocdmMakeTable::Begin $1 $2 $3"
    export ALIROOT_FORCE_COREDUMP=1
    local inFile=${1}
    local inFlag=${2}
    local outFile=${3}
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
    root.exe -l -q -b ${tmpscript}
    sleep 60 && rm ${tmpscript} &
    echo Filtering alien cache part 
    cat ${outFile} | sed -e s_"alien://?User=?DBFolder="_"alien://"_g -e s_"?SE=default?CacheFolder=?OperateDisconnected=1?CacheSize=1073741824?CleanupInterval=0"__g > ${outFile}.tmp
    cp  ${outFile}.tmp ${outFile}
    alilog_success "ocdmMakeTable $1 $2 $3"
    return 1
}


dumpObject(){
   [[ $# -ne 4 ]] && helpCat <<HELP_USAGE
   # dumpOCDBXML - Dump OCDB file as an xml file
   # Input parameters
        * param1: path
        * param2: Input object name
        * param3: type of the dump (POCDB - Print, DOCDB - Dump, "XML",  "MI" recursive dump )
        * param4: output file name
   # Example usage:
        # dump object TPC*Calib*RecoParam  from OCDB snapshot
            dumpObject  alien:///alice/cern.ch/user/p/pwg_pp/JIRA/ALIROOT-8240/TPC/Calib/RecoParam/Run0_999999999_v0_s2.root    "object" "XML"  Run0_999999999_v0_s2.xml
            dumpObject  alien:///alice/cern.ch/user/p/pwg_pp/JIRA/ALIROOT-8240/TPC/Calib/RecoParam/Run0_999999999_v0_s2.root    "object" "MI"  Run0_999999999_v0_s2.mi
            dumpObject  alien:///alice/cern.ch/user/p/pwg_pp/JIRA/ALIROOT-8240/TPC/Calib/RecoParam/Run0_999999999_v0_s2.root    "object" "POCDB"  Run0_999999999_v0_s2.pocdb
            dumpObject  alien:///alice/cern.ch/user/p/pwg_pp/JIRA/ALIROOT-8240/TPC/Calib/RecoParam/Run0_999999999_v0_s2.root    "object" "POCDB"  Run0_999999999_v0_s2.pocdb
            dumpObject alien:///alice/sim/2020/LHC120c1/OCDB/297479/iontail1.0crosstalk0.0/OCDBsim.root "TPC*Calib*RecoParam" "XML" iontail1.0crosstalk0.0.xml
            dumpObject alien:///alice/sim/2020/LHC120c1/OCDB/297479/iontail0.8crosstalk0.2/OCDBsim.root "TPC*Calib*RecoParam" "XML" iontail0.8crosstalk0.2.xml

HELP_USAGE
    [[ $# -ne 4 ]] && return
    alilog_info "dumpObject BEGIN $1 $2 $3 $4"
    export ALIROOT_FORCE_COREDUMP=1
    declare -a dumpOptions=("docdb" "xml" "pocdb" "mi")  # supported options. options are converted into lower case
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

    if [ $fobject == "object" ] ; then
        cat > ${tmpscript} <<HEREDOC
        {
            AliOCDBtoolkit::DumpOCDBFile("${inFile}","${outFile}",1,"${ftype}");
        }
HEREDOC
    else
        cat > ${tmpscript} <<HEREDOC
        {
            AliOCDBtoolkit::DumpOCDBFile("${inFile}","${fobject}","${outFile}",0,"${ftype}");
        }
HEREDOC
    fi

    root.exe -l -q -b ${tmpscript}
    sleep 60 && rm ${tmpscript} &
    alilog_success "dumpObject ${inFile} ${fobject} ${ftype} ${outFile}"
    return 1
}

diffObject(){
   [[ $# -ne 5 ]] && helpCat <<HELP_USAGE
   # diffObject
   # Input parameters
        * param1: path0
        * param2: path1
        * param3: obj name
        * param4:  type of the dump (xml or MI recursive dump )
    #  Output:
        * param5:  output diff file name
    #Example:
        diffObject /cvmfs/alice-ocdb.cern.ch/calibration/data/2018/OCDB/TPC/Calib/PadNoise/Run297459_999999999_v300_s0.root /cvmfs/alice-ocdb.cern.ch/calibration/data/2018/OCDB/TPC/Calib/PadNoise/Run294932_999999999_v257_s0.root AliCDBEntry POCDB PadNoise.diff
HELP_USAGE
    export ALIROOT_FORCE_COREDUMP=1
    [ $# -lt 5 ] && return 1
    alilog_info "diffObject::Beginn $1 $2 $3 $4 $5"
    local inFile1=${1}
    local inFile2=${2}
    local fobject=${3}
    local ftype=${4}
    local outFile=${5}
    #shift 5
    local tmp1=$(mktemp)
    local tmp2=$(mktemp)
    options=( "XML"  "MI" "POCDB" "DOCDB" )
    [[ ${options[@]} != *${ftype}* ]] && alilog_error "diffObject unsuported format $ftype" && return 1
    dumpObject ${inFile1} ${fobject} ${ftype} dumpObject1.txt
    dumpObject ${inFile2} ${fobject} ${ftype} dumpObject2.txt
    diff  dumpObject1.txt dumpObject2.txt >${outFile}
    alilog_success "diffObject $1 $2 $3 $4 $5"
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
[[ $# -ne 1 ]] && helpCat <<HELP_USAGE
   # ocdbMapInfo
        Dump basic AliRoot session information
        Tag values which are considered to be written to the OCDB metadata
HELP_USAGE
    #
    # Dump basic AliRoot session information
    # Tag values which are considered to be written to the OCDB metadata
    # (source $ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh; ocdbMapInfo )
    echo DATE":"`date`
    echo USER":"${USER}
    echo HOSTNAME":"$HOSTNAME
    echo WorkingDirectory":"`pwd`
    echo ROOTVERSION":"`root-config --version`
    echo ALICE_ROOT":"$(git -C ${AliRoot_SRC}/ describe)
    echo ALICE_PHYSICS":"$(git -C ${AliPhysics_SRC}/ describe)
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