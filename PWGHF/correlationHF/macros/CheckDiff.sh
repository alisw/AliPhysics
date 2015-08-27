#!/bin/bash

echo "####################### "
echo " CHECKING DoAverages.sh " 
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoAverages.sh DoAverages.sh
echo "####################### "
echo " CHECKING DoSubtractFD.sh"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoSubtractFD.sh DoSubtractFD.sh
echo "####################### "
echo " CHECKING DrawPlotInSubtractFD.sh"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DrawPlotInSubtractFD.sh DrawPlotInSubtractFD.sh
echo "####################### "
echo " CHECKING DoFit.sh"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoFit.sh DoFit.sh
echo "####################### "
echo " CHECKING RunFeedown_pPb_Dplus.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/RunFeedown_pPb_Dplus.C RunFeedown_pPb_Dplus.C
echo "####################### "
echo " CHECKING  RunFeedown_pPb_Dstar.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/RunFeedown_pPb_Dstar.C RunFeedown_pPb_Dstar.C
echo "####################### "
echo " CHECKING RunFeedown_pPb_Dzero.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/RunFeedown_pPb_Dzero.C RunFeedown_pPb_Dzero.C
echo "####################### "
echo " CHECKING RunFeedown_pp_Dplus.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/RunFeedown_pp_Dplus.C RunFeedown_pp_Dplus.C
echo "####################### "
echo " CHECKING RunFeedown_pp_Dstar.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/RunFeedown_pp_Dstar.C RunFeedown_pp_Dstar.C
echo "####################### "
echo " CHECKING RunFeedown_pp_Dzero.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/RunFeedown_pp_Dzero.C RunFeedown_pp_Dzero.C
echo "####################### "
echo " CHECKING MakeAverageDhCorrel.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/MakeAverageDhCorrel.C MakeAverageDhCorrel.C
echo "####################### "
echo " CHECKING DoPlotInSingleCanvas.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoPlotInSingleCanvas.C DoPlotInSingleCanvas.C
echo "####################### "
echo " CHECKING DoPlotCompare1GeVpPb.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoPlotCompare1GeVpPb.C DoPlotCompare1GeVpPb.C
echo "####################### "
echo " CHECKING DoPlotCompare1GeVpp.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoPlotCompare1GeVpp.C DoPlotCompare1GeVpp.C
echo "####################### "
echo " CHECKING DoPlotComparedot3to1pPb.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoPlotComparedot3to1pPb.C DoPlotComparedot3to1pPb.C
echo "####################### "
echo " CHECKING DoPlotComparedot3to1pp.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoPlotComparedot3to1pp.C DoPlotComparedot3to1pp.C
echo "####################### "
echo " CHECKING DoComparison_pPbVsMC.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoComparison_pPbVsMC.C DoComparison_pPbVsMC.C
echo "####################### "
echo " CHECKING DoComparison_ppVsMC.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoComparison_ppVsMC.C DoComparison_ppVsMC.C
echo "####################### "
echo " CHECKING DoComparison_ppVspPb.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoComparison_ppVspPb.C DoComparison_ppVspPb.C
echo "####################### "
echo " CHECKING FitSystematicsAverage_pp.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/FitSystematicsAverage_pp.C FitSystematicsAverage_pp.C
echo "####################### "
echo " CHECKING  FitSystematicsAverage_pPb.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/FitSystematicsAverage_pPb.C FitSystematicsAverage_pPb.C
echo "####################### "
echo " CHECKING FitPlots.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/FitPlots.C FitPlots.C
echo "####################### "
echo " CHECKING SubtractFD.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/SubtractFD.C SubtractFD.C
echo "####################### "
echo " CHECKING DoNiceSpecieComparisonPlot.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoNiceSpecieComparisonPlot.C DoNiceSpecieComparisonPlot.C
echo "####################### "
echo " CHECKING CompareFitResults.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/CompareFitResults.C CompareFitResults.C
echo "####################### "
echo " CHECKING DoNiceFitPlots.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoNiceFitPlots.C DoNiceFitPlots.C
echo "####################### "
echo " CHECKING DoComparison_ppVsMCallPanelsNew.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoComparison_ppVsMCallPanelsNew.C DoComparison_ppVsMCallPanelsNew.C
echo "####################### "
echo " CHECKING DoPlotInSingleCanvasNoSpaces.C"
echo "####################### "
diff ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoPlotInSingleCanvasNoSpaces.C DoPlotInSingleCanvasNoSpaces.C


