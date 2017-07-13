#!/bin/bash

declare macroPP="${HFCJlocalCodeDir}/FitSystematicsAverage_pp.C"
declare macroPPb="${HFCJlocalCodeDir}/FitSystematicsAverage_pPb.C"
declare macroPPb2016="${HFCJlocalCodeDir}/FitSystematicsAverage_pPb2016.C"
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


elif [ ${system} = 2 ]; then
    root -b &>outFitMCPPb03to99had.log <<EOF
.L ${macroPPb2016}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SystematicsMC_pPb_03to99had($refl,"$templNameRoot")
EOF

    root -b &>outFitMCPPb03to1had.log <<EOF
.L ${macroPPb2016}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SystematicsMC_pPb_03to1had($refl,"$templNameRoot")
EOF

    root -b &>outFitMCPPb1to99had.log <<EOF
.L ${macroPPb2016}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SystematicsMC_pPb_1to99had($refl,"$templNameRoot")
EOF

    root -b &>outFitMCPPb2to99had.log <<EOF
.L ${macroPPb2016}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SystematicsMC_pPb_2to99had($refl,"$templNameRoot")
EOF

    root -b &>outFitMCPPb3to99had.log <<EOF
.L ${macroPPb2016}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SystematicsMC_pPb_3to99had($refl,"$templNameRoot")
EOF

    root -b &>outFitMCPPb1to2had.log <<EOF
.L ${macroPPb2016}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SystematicsMC_pPb_1to2had($refl,"$templNameRoot")
EOF

    root -b &>outFitMCPPb2to3had.log <<EOF
.L ${macroPPb2016}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SystematicsMC_pPb_2to3had($refl,"$templNameRoot")
EOF
  
fi
