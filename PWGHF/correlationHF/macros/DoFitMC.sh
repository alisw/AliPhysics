#!/bin/bash

declare macroPP="${HFCJlocalCodeDir}/FitSystematicsAverage_pp.C"
declare macroPPb="${HFCJlocalCodeDir}/FitSystematicsAverage_pPb.C"
declare -i system=$1
declare -i refl=$2
declare -i avopt=$3
declare inputdir=$4
declare outputdir=$5
declare templNameRoot=$6
declare avstring="Weighted"
declare saveawayside=$7

if [ ${avopt} = 1 ]; then
    avstring="Arithmetic"
fi

if [ ${system} = 0 ]; then

    root -b &>outFitMCPPlowpthad.log <<EOF
.L ${macroPP}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SystematicsMC_pp_lowpthad($refl,"$templNameRoot")
EOF

    root -b &>outFitMCPPhighpthad.log <<EOF
.L ${macroPP}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SystematicsMC_pp_highpthad($refl,"$templNameRoot")
EOF


    root -b &>outFitMCPPintegrpthad.log <<EOF
.L ${macroPP}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SystematicsMC_pp_integrpthad($refl,"$templNameRoot")
EOF

    
elif [ ${system} = 1 ]; then
    root -b &>outFitMCPPblowpthad.log <<EOF
.L ${macroPPb}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SystematicsMC_pPb_lowpthad($refl,"$templNameRoot")
EOF

    root -b &>outFitMCPPbhighpthad.log <<EOF
.L ${macroPPb}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SystematicsMC_pPb_highpthad($refl,"$templNameRoot")
EOF


    root -b &>outFitMCPPbintegrpthad.log <<EOF
.L ${macroPPb}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SystematicsMC_pPb_integrpthad($refl,"$templNameRoot")
EOF
  
fi