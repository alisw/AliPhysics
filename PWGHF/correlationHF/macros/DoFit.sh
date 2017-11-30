#!/bin/bash

declare macroPP="${HFCJlocalCodeDir}/FitSystematicsAverage_pp.C"
declare macroPPb="${HFCJlocalCodeDir}/FitSystematicsAverage_pPb.C"
declare macroPPb2016="${HFCJlocalCodeDir}/FitSystematicsAverage_pPb2016.C"
declare -i system=$1
declare -i refl=$2
declare -i avopt=$3
declare inputdir=$4
declare outputdir=$5
declare avstring="Weighted"
declare -i includev2=$6
declare saveawayside=$7
declare -i plotv2sep=1
declare -i v2had=8
declare -i v2D=5
#the following for pPb2016 (percentage points)
declare -i v2had03to1=4
declare -i v2had03to99=6
declare -i v2had1to99=9
declare -i v2had2to99=9
declare -i v2had3to99=10
declare -i v2had1to2=8
declare -i v2had2to3=9
declare -i v2D3to5=5
declare -i v2D5to8=3
declare -i v2D8to16=2
declare -i v2D16to24=2

if [ ${avopt} = 1 ]; then
    avstring="Arithmetic"
fi

if [ ${system} = 0 ]; then
    root -b &>outFitPPlowpthad.log <<EOF
.L ${macroPP}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
 Systematics_pp_lowpthad(${refl})
EOF

root -b &>outFitPPhighpthad.log <<EOF
.L ${macroPP}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
Systematics_pp_highpthad(${refl})
EOF

root -b &>outFitPPintegratedpthad.log <<EOF
.L ${macroPP}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
Systematics_pp_integratedpthad(${refl})
EOF

elif [ ${system} = 1 ]; then
    root -b &>outFitPPblowpthad.log <<EOF
.L ${macroPPb}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SetV2values(${includev2},0.01*((Double_t)${v2had}),0.01*((Double_t)${v2D}),${plotv2sep})
Systematics_pPb_lowpthad(${refl})
EOF

root -b &>outFitPPbhighpthad.log <<EOF
.L ${macroPPb}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SetV2values(${includev2},0.01*((Double_t)${v2had}),0.01*((Double_t)${v2D}),${plotv2sep})
Systematics_pPb_highpthad(${refl})
EOF

root -b &>outFitPPbintegratedpthad.log <<EOF
.L ${macroPPb}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SetV2values(${includev2},0.01*((Double_t)${v2had}),0.01*((Double_t)${v2D}),${plotv2sep})
Systematics_pPb_integratedpthad(${refl})
EOF

elif [ ${system} = 2 ]; then
root -b &>outFitPPb03to99had.log <<EOF
.L ${macroPPb2016}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SetV2values(${includev2},0.01*((Double_t)${v2had03to99}),0.01*((Double_t)${v2D3to5}),0.01*((Double_t)${v2D5to8}),0.01*((Double_t)${v2D8to16}),0.01*((Double_t)${v2D16to24}),${plotv2sep})
Systematics_pPb_03to99had(${refl})
EOF

    root -b &>outFitPPb03to1had.log <<EOF
.L ${macroPPb2016}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SetV2values(${includev2},0.01*((Double_t)${v2had03to1}),0.01*((Double_t)${v2D3to5}),0.01*((Double_t)${v2D5to8}),0.01*((Double_t)${v2D8to16}),0.01*((Double_t)${v2D16to24}),${plotv2sep})
Systematics_pPb_03to1had(${refl})
EOF

root -b &>outFitPPb1to99had.log <<EOF
.L ${macroPPb2016}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SetV2values(${includev2},0.01*((Double_t)${v2had1to99}),0.01*((Double_t)${v2D3to5}),0.01*((Double_t)${v2D5to8}),0.01*((Double_t)${v2D8to16}),0.01*((Double_t)${v2D16to24}),${plotv2sep})
Systematics_pPb_1to99had(${refl})
EOF


root -b &>outFitPPb2to99had.log <<EOF
.L ${macroPPb2016}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SetV2values(${includev2},0.01*((Double_t)${v2had2to99}),0.01*((Double_t)${v2D3to5}),0.01*((Double_t)${v2D5to8}),0.01*((Double_t)${v2D8to16}),0.01*((Double_t)${v2D16to24}),${plotv2sep})
Systematics_pPb_2to99had(${refl})
EOF


root -b &>outFitPPb3to99had.log <<EOF
.L ${macroPPb2016}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SetV2values(${includev2},0.01*((Double_t)${v2had3to99}),0.01*((Double_t)${v2D3to5}),0.01*((Double_t)${v2D5to8}),0.01*((Double_t)${v2D8to16}),0.01*((Double_t)${v2D16to24}),${plotv2sep})
Systematics_pPb_3to99had(${refl})
EOF


root -b &>outFitPPb1to2had.log <<EOF
.L ${macroPPb2016}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SetV2values(${includev2},0.01*((Double_t)${v2had1to2}),0.01*((Double_t)${v2D3to5}),0.01*((Double_t)${v2D5to8}),0.01*((Double_t)${v2D8to16}),0.01*((Double_t)${v2D16to24}),${plotv2sep})
Systematics_pPb_1to2had(${refl})
EOF


root -b &>outFitPPb2to3had.log <<EOF
.L ${macroPPb2016}
SetInputPath("${inputdir}")
SetOutputPath("${outputdir}")
SetAverageString("${avstring}")
SetSaveAwaySidePlots(${saveawayside})
SetV2values(${includev2},0.01*((Double_t)${v2had2to3}),0.01*((Double_t)${v2D3to5}),0.01*((Double_t)${v2D5to8}),0.01*((Double_t)${v2D8to16}),0.01*((Double_t)${v2D16to24}),${plotv2sep})
Systematics_pPb_2to3had(${refl})
EOF
fi
