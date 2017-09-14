#!/bin/bash

#declare macrosDir="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015March13/Macros"
declare -i collsyst=$1
declare -i mesonIndex=$2
declare -a macrosPP=( "RunFeedown_pp_Dzero" "RunFeedown_pp_Dstar" "RunFeedown_pp_Dplus" )
declare -a macrosPPb=( "RunFeedown_pPb_Dzero" "RunFeedown_pPb_Dstar" "RunFeedown_pPb_Dplus" )
declare fpromptfile=$3
declare templatedir=$4
declare inputfiledir=$5
declare inputfileroot=$6
declare -i localcode=$7
declare suffixTemplSystm=$8
declare -i subtrMCclos=$9
declare dirmacroRun="${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros"
declare dirmacroFD="${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros"
if [ ${localcode} = 1 ]; then
    dirmacroRun=${HFCJlocalCodeDir}
fi
if [ ${localcode} = 2 ]; then
    dirmacroFD=${HFCJlocalCodeDir}
fi
if [ ${localcode} = 3 ]; then
    dirmacroRun=${HFCJlocalCodeDir}
    dirmacroFD=${HFCJlocalCodeDir}
fi

if [ $collsyst = 0 ]; then
    echo "DoSubtractFD: subtracting FD for meson $mesonIndex in coll system $collsyst"
    root -b <<EOF &> out.log
Printf("inside root");
.L ${dirmacroRun}/${macrosPP[$mesonIndex]}.C
SetFDmacroDirectory("${dirmacroFD}")
SetFpromptInputFile("$fpromptfile")
SetTemplateDir("${templatedir}")
SetDirectoryInputFiles("${inputfiledir}")
SetInputFileNameRoot("$inputfileroot")
SetFDtemplateSystemString("${suffixTemplSystm}")
//cout<<"file: "<<$file<<endl
//Printf("Analyzing file: %s",${file})
Printf("Coll syst=%d",${collsyst})
Printf("Meson: %d",$mesonIndex)
${macrosPP[$mesonIndex]}()
.q
EOF
    
echo "DoSubtractFD: done"
elif [ $collsyst = 1 ]; then
    echo "DoSubtractFD: subtracting FD for meson $mesonIndex in coll system $collsyst"
    root -b <<EOF &> out.log
Printf("inside root");
.L ${dirmacroRun}/${macrosPPb[$mesonIndex]}.C
SetFDmacroDirectory("${dirmacroFD}")
SetFpromptInputFile("$fpromptfile")
SetTemplateDir("${templatedir}")
SetDirectoryInputFiles("${inputfiledir}")
SetInputFileNameRoot("$inputfileroot")
SetFDtemplateSystemString("${suffixTemplSystm}")
//cout<<"file: "<<$file<<endl
//Printf("Analyzing file: %s",${file})
Printf("Coll syst=%d",${collsyst})
Printf("Meson: %d",$mesonIndex)
${macrosPPb[$mesonIndex]}($collsyst,$subtrMCclos)
.q
EOF

echo "DoSubtractFD: done"
elif [ $collsyst = 2 ]; then
    echo "DoSubtractFD: subtracting FD for meson $mesonIndex in coll system $collsyst (2016)"
    root -b <<EOF &> out.log
Printf("inside root");
.L ${dirmacroRun}/${macrosPPb[$mesonIndex]}.C
SetFDmacroDirectory("${dirmacroFD}")
SetFpromptInputFile("$fpromptfile")
SetTemplateDir("${templatedir}")
SetDirectoryInputFiles("${inputfiledir}")
SetInputFileNameRoot("$inputfileroot")
SetFDtemplateSystemString("${suffixTemplSystm}")
//cout<<"file: "<<$file<<endl
//Printf("Analyzing file: %s",${file})
Printf("Coll syst=%d",${collsyst})
Printf("Meson: %d",$mesonIndex)
${macrosPPb[$mesonIndex]}($collsyst,$subtrMCclos)
.q
EOF
echo "DoSubtractFD: done"  
fi

exit 0