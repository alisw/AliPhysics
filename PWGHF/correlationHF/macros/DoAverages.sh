#!/bin/bash

declare -i collsyst=$1
declare -i ptTrig=$2
declare -i ptAssoc=$3
declare -i averOpt=$4
declare -i reflOpt=$5
declare -i baselineOpt=$6
declare dirDzero=$7
declare dirDstar=$8
declare dirDplus=$9
declare -i localCode=${10}
declare -ai ptTrigMin=(3 5 8 16) 
declare -ai ptTrigMax=(5 8 16 24) 

declare -ai ptAssocMin=(3 3 10 20 30 10 20) #have to divide by 10 
declare -ai ptAssocMax=(10 990 990 990 990 20 30)  #have to divide by 10

declare -a ptAssocMinStr=("0.3" "0.3" "1.0" "2.0" "3.0" "1.0" "2.0") 
declare -a ptAssocMaxStr=("1.0" "99.0" "99.0" "99.0" "99.0" "2.0" "3.0") 

declare -ai plotheightpp=(5 8 3)
declare -ai plotheightpPb=(10 14 5 2 1 3 2)
declare -i heightPlot=${plotheightpp[$ptAssoc]}
declare -ai maxrange=(10 15 15)
declare -i year=2010
declare collsyststr="pp"
declare dirmacroAverage="$ALICE_PHYSICS/../src/PWGHF/correlationHF/macros"

if [ ${localCode} = 1 ]; then
    dirmacroAverage=${HFCJlocalCodeDir}
fi

if [ ${collsyst} = 1 ];then
    year=2013
    heightPlot=${plotheightpPb[$ptAssoc]}
    collsyststr="pPb"
fi

if [ ${collsyst} = 2 ];then
    year=2016
    heightPlot=${plotheightpPb[$ptAssoc]}
    collsyststr="pPb2016"
fi

declare avString="Weighted"
if [ ${averOpt} = 1 ];then
    avString="Arithmetic"
fi
root -b <<EOF &> outColl${collsyst}PtD${ptTrig}PtAssoc${ptAssoc}.log
Printf("inside root");
.L ${dirmacroAverage}/MakeAverageDhCorrel.C
//cout<<"file: "<<$file<<endl
//Printf("Analyzing file: %s",${file})
Printf("Trig pt: %f to %f",(Double_t)(${ptTrigMin[$ptTrig]}),(Double_t)(${ptTrigMax[$ptTrig]}));
Printf("Min pt assoc: %f",(Double_t)(${ptAssocMin[$ptAssoc]})/10.);
Printf("Coll syst=%d",${collsyst})
SetFitPlotMacroPath("${dirmacroAverage}")
MakeAverage(${ptTrigMin[$ptTrig]},${ptTrigMax[$ptTrig]},${ptAssocMin[$ptAssoc]}/10.,${ptAssocMax[$ptAssoc]}/10,$collsyst,$year,$reflOpt,$heightPlot,$averOpt,$baselineOpt,"$dirDzero","$dirDplus","$dirDstar","")

.q
EOF


if [ ${collsyst} = 2 ];then
root -b <<EOF &> outStyleColl${collsyst}PtD${ptTrig}PtAssoc${ptAssoc}.log
Printf("inside root");
.L ${dirmacroAverage}/MakeAverageDhCorrel.C
OpenOutputFileAndDraw("${avString}AveragepPbDzeroDstarDplus${ptTrigMin[$ptTrig]}to${ptTrigMax[$ptTrig]}_assoc${ptAssocMinStr[$ptAssoc]}to${ptAssocMaxStr[$ptAssoc]}.root",${ptTrigMin[$ptTrig]},${ptTrigMax[$ptTrig]},"D",${collsyst},${ptAssocMin[$ptAssoc]}/10.,${ptAssocMax[$ptAssoc]}/10,1,"${avString}",$heightPlot)
.q
EOF
else 
root -b <<EOF &> outStyleColl${collsyst}PtD${ptTrig}PtAssoc${ptAssoc}.log
Printf("inside root");
.L ${dirmacroAverage}/MakeAverageDhCorrel.C
OpenOutputFileAndDraw("${avString}Average${collsyststr}DzeroDstarDplus${ptTrigMin[$ptTrig]}to${ptTrigMax[$ptTrig]}_assoc${ptAssocMinStr[$ptAssoc]}to${ptAssocMaxStr[$ptAssoc]}.root",${ptTrigMin[$ptTrig]},${ptTrigMax[$ptTrig]},"D",${collsyst},${ptAssocMin[$ptAssoc]}/10.,${ptAssocMax[$ptAssoc]}/10,1,"${avString}",$heightPlot)
.q
EOF
fi
