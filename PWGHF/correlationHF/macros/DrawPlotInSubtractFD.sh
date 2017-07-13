#!/bin/bash

#declare macrosDir="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015March13/Macros"
declare file=$1
declare -i mesonIndex=$2
declare -i collsyst=$3
declare -i ptTrig=$4
declare -i ptAssoc=$5
declare -i reflect=$6
declare -i rebinAzi=$7
declare -i localcode=$8
declare dirmacroFD="${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros"

declare -ai ptTrigMin=(3 5 8 16) 
declare -ai ptTrigMax=(5 8 16 24) 

declare -ai ptAssocMin=(3 3 10 20 30 10 20) #have to divide by 10 
declare -ai ptAssocMax=(10 990 990 990 990 20 30)  #have to divide by 10

declare -ai maxrange=(10 15 15)
if [ ${localcode} = 1 ]; then
    dirmacroFD=${HFCJlocalCodeDir}
fi
echo "Analyzing file: $file"
#.L ${macrosDir}/SubtractFDtest.C
root -b <<EOF &> out.log
Printf("inside root");
.L ${dirmacroFD}/SubtractFD.C
//cout<<"file: "<<$file<<endl
//Printf("Analyzing file: %s",${file});
Printf("Meson: %d",$mesonIndex);
Printf("Trig pt: %f to %f",(Double_t)(${ptTrigMin[$ptTrig]}),(Double_t)(${ptTrigMax[$ptTrig]}));
Printf("Min pt assoc: %f",(Double_t)(${ptAssocMin[$ptAssoc]})/10.);
Printf("Coll syst=%d",${collsyst})
OpenOutputFileAndDrawReflect("$file",(Double_t)(${ptTrigMin[$ptTrig]}),(Double_t)(${ptTrigMax[$ptTrig]}),${mesonIndex},(Double_t)(${ptAssocMin[$ptAssoc]})/10.,(Double_t)(${ptAssocMax[$ptAssoc]})/10.,1,${collsyst},${reflect},${maxrange[${collsyst}]},${rebinAzi})
.q
EOF
