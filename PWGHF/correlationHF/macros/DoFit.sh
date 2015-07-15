#!/bin/bash

declare macroPP="${HFCJlocalCodeDir}/FitSystematicsAverage_pp.C"
declare macroPPb="${HFCJlocalCodeDir}/FitSystematicsAverage_pPb.C"
declare -i system=$1
declare -i refl=$2
declare -i avopt=$3
declare inputdir=$4
declare outputdir=$5
declare avstring="Weighted"
declare -i includev2=$6
declare -i plotv2sep=1
declare -i v2had=1
declare -i v2D=1


if [ ${avopt} = 1 ]; then
    avstring="Arithmetic"
fi

if [ ${system} = 0 ]; then
    root -b &>outFitPPlowpthad.log <<EOF
.L ${macroPP}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
 Systematics_pp_lowpthad(${refl})
EOF

root -b &>outFitPPhighpthad.log <<EOF
.L ${macroPP}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
Systematics_pp_highpthad(${refl})
EOF

root -b &>outFitPPintegratedpthad.log <<EOF
.L ${macroPP}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
Systematics_pp_integratedpthad(${refl})
EOF

elif [ ${system} = 1 ]; then
    root -b &>outFitPPblowpthad.log <<EOF
.L ${macroPPb}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
Systematics_pPb_lowpthad(${refl})
SetV2values(${includev2},0.1*((Double_t)v2had),0.1*((Double_t)v2D),${plotv2sep})
EOF

root -b &>outFitPPbhighpthad.log <<EOF
.L ${macroPPb}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
Systematics_pPb_highpthad(${refl})
SetV2values(${includev2},0.1*((Double_t)v2had),0.1*((Double_t)v2D),${plotv2sep})
EOF

root -b &>outFitPPbintegratedpthad.log <<EOF
.L ${macroPPb}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
Systematics_pPb_integratedpthad(${refl})
SetV2values(${includev2},0.1*((Double_t)v2had),0.1*((Double_t)v2D),${plotv2sep})
EOF
fi