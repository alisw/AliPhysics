#!/bin/bash
declare pathRef=$1

if [ ! -d ${pathRef} ];then
    echo "REFERENCE DIRECTORY FOR CHECKING CODE CONSISTENCY WILL BE"
    echo "${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/"
    pathRef="${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/"
fi


echo "####################### "
echo " CHECKING CheckDiff.sh " 
echo "####################### "
diff -b ${pathRef}/CheckDiff.sh CheckDiff.sh
echo "####################### "
echo " CHECKING ProducePlotChain.sh " 
echo "####################### "
diff -b ${pathRef}/ProducePlotChain.sh ProducePlotChain.sh
echo "####################### "
echo " CHECKING DoAverages.sh " 
echo "####################### "
diff -b ${pathRef}/DoAverages.sh DoAverages.sh
echo "####################### "
echo " CHECKING DoSubtractFD.sh"
echo "####################### "
diff -b ${pathRef}/DoSubtractFD.sh DoSubtractFD.sh
echo "####################### "
echo " CHECKING DrawPlotInSubtractFD.sh"
echo "####################### "
diff -b ${pathRef}/DrawPlotInSubtractFD.sh DrawPlotInSubtractFD.sh
echo "####################### "
echo " CHECKING DoFit.sh"
echo "####################### "
diff -b ${pathRef}/DoFit.sh DoFit.sh
echo "####################### "
echo " CHECKING DoFitMC.sh"
echo "####################### "
diff -b ${pathRef}/DoFitMC.sh DoFitMC.sh
echo "####################### "
echo " CHECKING RunFeedown_pPb_Dplus.C"
echo "####################### "
diff -b ${pathRef}/RunFeedown_pPb_Dplus.C RunFeedown_pPb_Dplus.C
echo "####################### "
echo " CHECKING  RunFeedown_pPb_Dstar.C"
echo "####################### "
diff -b ${pathRef}/RunFeedown_pPb_Dstar.C RunFeedown_pPb_Dstar.C
echo "####################### "
echo " CHECKING RunFeedown_pPb_Dzero.C"
echo "####################### "
diff -b ${pathRef}/RunFeedown_pPb_Dzero.C RunFeedown_pPb_Dzero.C
echo "####################### "
echo " CHECKING RunFeedown_pp_Dplus.C"
echo "####################### "
diff -b ${pathRef}/RunFeedown_pp_Dplus.C RunFeedown_pp_Dplus.C
echo "####################### "
echo " CHECKING RunFeedown_pp_Dstar.C"
echo "####################### "
diff -b ${pathRef}/RunFeedown_pp_Dstar.C RunFeedown_pp_Dstar.C
echo "####################### "
echo " CHECKING RunFeedown_pp_Dzero.C"
echo "####################### "
diff -b ${pathRef}/RunFeedown_pp_Dzero.C RunFeedown_pp_Dzero.C
echo "####################### "
echo " CHECKING MakeAverageDhCorrel.C"
echo "####################### "
diff -b ${pathRef}/MakeAverageDhCorrel.C MakeAverageDhCorrel.C
echo "####################### "
echo " CHECKING DoPlotInSingleCanvas.C"
echo "####################### "
diff -b ${pathRef}/DoPlotInSingleCanvas.C DoPlotInSingleCanvas.C
echo "####################### "
echo " CHECKING DoPlotCompare1GeVpPb.C"
echo "####################### "
diff -b ${pathRef}/DoPlotCompare1GeVpPb.C DoPlotCompare1GeVpPb.C
echo "####################### "
echo " CHECKING DoPlotCompare1GeVpp.C"
echo "####################### "
diff -b ${pathRef}/DoPlotCompare1GeVpp.C DoPlotCompare1GeVpp.C
echo "####################### "
echo " CHECKING DoPlotComparedot3to1pPb.C"
echo "####################### "
diff -b ${pathRef}/DoPlotComparedot3to1pPb.C DoPlotComparedot3to1pPb.C
echo "####################### "
echo " CHECKING DoPlotComparedot3to1pp.C"
echo "####################### "
diff -b ${pathRef}/DoPlotComparedot3to1pp.C DoPlotComparedot3to1pp.C
echo "####################### "
echo " CHECKING DoComparison_pPbVsMC.C"
echo "####################### "
diff -b ${pathRef}/DoComparison_pPbVsMC.C DoComparison_pPbVsMC.C
echo "####################### "
echo " CHECKING DoComparison_ppVsMC.C"
echo "####################### "
diff -b ${pathRef}/DoComparison_ppVsMC.C DoComparison_ppVsMC.C
echo "####################### "
echo " CHECKING DoComparison_ppVspPb.C"
echo "####################### "
diff -b ${pathRef}/DoComparison_ppVspPb.C DoComparison_ppVspPb.C
echo "####################### "
echo " CHECKING FitSystematicsAverage_pp.C"
echo "####################### "
diff -b ${pathRef}/FitSystematicsAverage_pp.C FitSystematicsAverage_pp.C
echo "####################### "
echo " CHECKING  FitSystematicsAverage_pPb.C"
echo "####################### "
diff -b ${pathRef}/FitSystematicsAverage_pPb.C FitSystematicsAverage_pPb.C
echo "####################### "
echo " CHECKING FitPlots.C"
echo "####################### "
diff -b ${pathRef}/FitPlots.C FitPlots.C
echo "####################### "
echo " CHECKING SubtractFD.C"
echo "####################### "
diff -b ${pathRef}/SubtractFD.C SubtractFD.C
echo "####################### "
echo " CHECKING DoNiceSpecieComparisonPlot.C"
echo "####################### "
diff -b ${pathRef}/DoNiceSpecieComparisonPlot.C DoNiceSpecieComparisonPlot.C
echo "####################### "
echo " CHECKING CompareFitResults.C"
echo "####################### "
diff -b ${pathRef}/CompareFitResults.C CompareFitResults.C
echo "####################### "
echo " CHECKING DoNiceFitPlots.C"
echo "####################### "
diff -b ${pathRef}/DoNiceFitPlots.C DoNiceFitPlots.C
echo "####################### "
echo " CHECKING DoComparison_ppVsMCallPanelsNew.C"
echo "####################### "
diff -b ${pathRef}/DoComparison_ppVsMCallPanelsNew.C DoComparison_ppVsMCallPanelsNew.C
echo "####################### "
echo " CHECKING DoComparison_ppVspPballPanels.C"
echo "####################### "
diff -b ${pathRef}/DoComparison_ppVspPballPanels.C DoComparison_ppVspPballPanels.C
echo "####################### "
echo " CHECKING DoPlotInSingleCanvasNoSpaces.C"
echo "####################### "
diff -b ${pathRef}/DoPlotInSingleCanvasNoSpaces.C DoPlotInSingleCanvasNoSpaces.C
echo "####################### "
echo " CHECKING CompareFitResultspPb2016.C"
echo "####################### "
diff -b ${pathRef}/CompareFitResultspPb2016.C CompareFitResultspPb2016.C
echo "####################### "
echo " CHECKING DoComparison_ppVspPballPanels2016.C"
echo "####################### "
diff -b ${pathRef}/DoComparison_ppVspPballPanels2016.C DoComparison_ppVspPballPanels2016.C
echo "####################### "
echo " CHECKING DoComparison_pPb2016VsMCallPanelsNew.C"
echo "####################### "
diff -b ${pathRef}/DoComparison_pPb2016VsMCallPanelsNew.C DoComparison_pPb2016VsMCallPanelsNew.C
echo "####################### "
echo " CHECKING FitSystematicsAverage_pPb2016.C"
echo "####################### "
diff -b ${pathRef}/FitSystematicsAverage_pPb2016.C FitSystematicsAverage_pPb2016.C
echo "####################### "
echo " CHECKING DoComparison_pPbVsMC2016.C"
echo "####################### "
diff -b ${pathRef}/DoComparison_pPbVsMC2016.C DoComparison_pPbVsMC2016.C
echo "####################### "
echo " CHECKING Restyle_pPb2016_Prel_Plots.C"
echo "####################### "
diff -b ${pathRef}/Restyle_pPb2016_Prel_Plots.C Restyle_pPb2016_Prel_Plots.C

