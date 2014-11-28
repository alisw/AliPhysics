#!/bin/bash
#
# Shell script to compare content of the OCDB entries.
#
# ocdbMakeTable() 
#       Usage: bash $inputFile $flag $outputFile
# dumpObject()
#       Usage: bash $inputFile $object_name $dump_type [XML/MI] $outfile
# diffObject
#       Usage: bash $inputFile1 $inputFile2 $object_name $dump_type [XML/MI] $outfile
# Example usage:  see example functions below

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
    return 1
}


dumpObject(){
#
#
#  Input:
#    $1 path
#    $2 obj name 
#    $3 type of the dump (XML or MI recursive dump )
#  Output:
#    $4 output file name   
    export ALIROOT_FORCE_COREDUMP=1
    if [ $# -lt 4 ] ; then
        echo "Usage: $0 \$inputFile \$object_name \$dump_type [XML/MI] \$outfile"
        return 1
    fi
    local inFile=${1}
    local fobject=${2}
    local ftype=${3}
    local outFile=${4}
    shift 4
#    if [ ! -f ${inFile} ] ; then 
#        echo ${inFile} not found!
#        return 1
#    fi
    if [ -f ${outFile} ] ; then 
        >${outFile}
    fi
    if [ ${ftype} = "XML" ] ; then
        isXML=kTRUE
    elif [ ${ftype} = "MI" ] ; then
        isXML=kFALSE
    else
        echo "option ${ftype} not recognized! Use \"XML\" or \"MI\"!"
        return 1
    fi
    tmpscript=$(mktemp)
    cat > ${tmpscript} <<HEREDOC
        {
            AliOCDBtoolkit::DumpOCDBFile("${inFile}","${outFile}",1,${isXML});
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
