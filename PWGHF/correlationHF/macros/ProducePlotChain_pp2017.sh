#!/bin/bash
#######################################################################################################
###
###    Script to running in sequence 
###         - FD subtraction of D-h correlations in pp and p-Pb
###         - reflection of plots in 0-pi range
###         - average of different D meson results
###         - visualization: 
###              - comparison of different meson results 
###              - comparison of pp and p-Pb results
###              - comparison of data with MC simulations
###          
###   Note that, with standard settings (but there is a dedicate flag to remove this feature) the 
###   macros and scripts are used are copied from PWGHF/correlationHF/macros directory in a working directory (HFCJlocalCodeDir)
###   which is exported and used throught the code: this is not done only for convenience (you can modify the macro locally and run again the script 
###   with the cp option =0) but mainly for storing the version of the code used to produce the results.
###
###   A part from the options (see comments below), only the paths of the D meson files before FD subtraction
###   and of the MC templates should be set (will be replaced soon with wget calls to retrieve the files from the TWiki
###
###   Andrea Rossi, andrea.rossi@cern.ch  
### 
#######################################################################################################

######## THE FOLLOWING LINES ARE THE ONLY VARIABLES TO BE CHANGED ACCORDING TO YOUR LOCAL PATHS #####

declare baseStartingDir="$PWD/ScriptOutput/"
export HFCJlocalCodeDir="$PWD"
# "/Users/administrator/soft/alisoft/aliphysics/master/PWGHF/correlationHF/macros"
declare templateDirPP="/home/colamaria/Scrivania/Codici_Ausiliari_Dh/Dhadron_Final_Output_pp2017_PAPER/Inputs/Templates_pp"
declare templateDirPPb="/home/colamaria/Scrivania/Codici_Ausiliari_Dh/Dhadron_Final_Output_pp2017_PAPER/Inputs/Templates_pPb"
declare -a templateDirSystemSuffix=( "none" "none" "none" "none" ) #### THIS IS KEPT JUST FOR BACKWARD COMPATIBILITY WITH OLD TEMPLATES! NO NEED TO TOUCH IT UNLESS YOU WANT TO USE OLD TEMPLATES
declare -a templateDir=( "$templateDirPP" "$templateDirPPb" "$templateDirPPb"  "$templateDirPP")

###these for loading beta values from MC, in case of fixed beta approach!
declare betaTemplDirPP="/home/colamaria/Scrivania/Codici_Ausiliari_Dh/Dhadron_Final_Output_pp2017_PAPER/Inputs/Templates_pp/FitResults/Trends_pp/PYTHIA8"
declare betaTemplDirPPb="/home/colamaria/Scrivania/Codici_Ausiliari_Dh/Dhadron_Final_Output_pp2017_PAPER/Inputs/Templates_pPb/FitResults/Trends_pPb/PYTHIA8"
declare -a betaTemplDir=( "$betaTemplDirPP" "$betaTemplDirPPb" "$betaTemplDirPPb"  "$betaTemplDirPP")

### the following is needed for hte comparison to MC (as well as MC fitting)
declare -a Nmccase=( 6 6 6 8 ) #warning! If you change this, the number of fitted MC models will change (and also their order!) I suggest not changing it but only the array below...
declare -a mccasePP=( 0 0 1 1 1 1 1 1 ) # according to CompareFitResults array: Perugia0, Perugia2010, Perugia2011, PYTHIA8, HERWIG, POHWEG+Perugia2011, POWHEG+Perugia2011 with EPS09, EPOS 3
declare -a mccasePPb=( 1 1 1 1 0 0 1 0 )
declare -a isreflectedMC=( 0 0 0 0 0 0 0 1 ) # used only to determine the fit range and the transverse region range, it does not however change the results. Only EPOS is already reflected
declare -a templRootNamepp=( "CorrelationPlotsPerugia0PtAveragefromC" "CorrelationPlotsPerugia2010PtAveragefromC" "CorrelationPlotsPerugia2011PtAveragefromC" "CorrelationPlotsPYTHIA8PtAveragefromC" "CorrelationPlotsHERWIGPtAveragefromC" "CorrelationPlotsPOWHEGPtAveragefromC" "CorrelationPlotsPOW_LOPtAveragefromC" "CorrelationPlotsEPOS3PtAveragefromC")
declare -a templRootNamepPb=( "CorrelationPlotsPerugia0wBoostPtAveragefromC" "CorrelationPlotsPerugia2010wBoostPtAveragefromC" "CorrelationPlotsPerugia2011wBoostPtAveragefromC" "CorrelationPlotsPYTHIA8wBoostPtAveragefromC" "CorrelationPlotsHERWIGPtAveragefromC" "CorrelationPlotsPOWHEGPtAveragefromC" "CorrelationPlotsPOW_LOPtAveragefromC" "CorrelationPlotsEPOS3PtAveragefromC")

declare puritytemplateDirpp="/home/colamaria/Scrivania/Codici_Ausiliari_Dh/Dhadron_Final_Output_pp2017_PAPER/Inputs/PurityTemplates_pp"
declare puritytemplateDirpPb="/home/colamaria/Scrivania/Codici_Ausiliari_Dh/Dhadron_Final_Output_pp2017_PAPER/Inputs/PurityTemplates_pPb"
declare -a puritytemplateDir=( "$puritytemplateDirpp" "$puritytemplateDirpPb" "$puritytemplateDirpPb" "$puritytemplateDirpp")

########## THE FOLLOWING DIRECTORIES SHOULD CONTAIN THE RESULTS BEFORE FD SUBTRACTION #####
declare dirppDzeroNotFDsubt="/home/colamaria/Scrivania/Codici_Ausiliari_Dh/Dhadron_Final_Output_pp2017_PAPER/Inputs/Dzero_pp"
declare dirpPbDzeroNotFDsubt="/home/colamaria/Scrivania/Codici_Ausiliari_Dh/Dhadron_Final_Output_pp2017_PAPER/Inputs/Dzero_pPb"
declare -a dirDzeroNotFDsubt=( "$dirppDzeroNotFDsubt" "$dirpPbDzeroNotFDsubt" "$dirpPbDzeroNotFDsubt" "$dirppDzeroNotFDsubt" )
declare -a fpromptfileDzero=( "HFPtSpectrum_Nb.root" "HFPtSpectrum_Nb.root" "HFPtSpectrum_Nb.root"  "HFPtSpectrum_Nb_Dzero.root")
declare -a filerootDzero=( "AzimCorrDistr_Dzero_Canvas_PtIntBins" "AzimCorrDistr_Dzero_Canvas_PtIntBins" "AzimCorrDistr_Dzero_Canvas_PtIntBins" "AzimCorrDistr_Dzero_Canvas_PtIntBins" )

declare dirppDstarNotFDsubt="/home/colamaria/Scrivania/Codici_Ausiliari_Dh/Dhadron_Final_Output_pp2017_PAPER/Inputs/Dstar_pp"
declare dirpPbDstarNotFDsubt="/home/colamaria/Scrivania/Codici_Ausiliari_Dh/Dhadron_Final_Output_pp2017_PAPER/Inputs/Dstar_pPb"
declare -a dirDstarNotFDsubt=( "$dirppDstarNotFDsubt" "$dirpPbDstarNotFDsubt" "$dirpPbDstarNotFDsubt" "$dirppDstarNotFDsubt" )
declare -a fpromptfileDstar=( "HFPtSpectrum_Nb.root" "HFPtSpectrum_Nb.root" "HFPtSpectrum_Nb.root"  "HFPtSpectrum_Nb_Dstar.root")
declare -a filerootDstar=( "AzimCorrDistr_Dstar_Canvas_PtIntBins" "AzimCorrDistr_Dstar_Canvas_PtIntBins" "AzimCorrDistr_Dstar_Canvas_PtIntBins" "AzimCorrDistr_Dstar_Canvas_PtIntBins" )

declare dirppDplusNotFDsubt="/home/colamaria/Scrivania/Codici_Ausiliari_Dh/Dhadron_Final_Output_pp2017_PAPER/Inputs/Dplus_pp"
declare dirpPbDplusNotFDsubt="/home/colamaria/Scrivania/Codici_Ausiliari_Dh/Dhadron_Final_Output_pp2017_PAPER/Inputs/Dplus_pPb"
declare -a dirDplusNotFDsubt=( "$dirppDplusNotFDsubt" "$dirpPbDplusNotFDsubt" "$dirpPbDplusNotFDsubt" "$dirppDplusNotFDsubt" )
declare -a fpromptfileDplus=( "HFPtSpectrum_Nb.root" "HFPtSpectrum_Nb.root" "HFPtSpectrum_Nb.root"  "HFPtSpectrum_Nb_Dplus.root")
declare -a filerootDplus=( "AzimCorrDistr_Dplus_Canvas_PtIntBins" "AzimCorrDistr_Dplus_Canvas_PtIntBins" "AzimCorrDistr_Dplus_Canvas_PtIntBins" "AzimCorrDistr_Dplus_Canvas_PtIntBins" )

###### THE FOLLOWING DIRECTORIES WILL BE USED ONLY IN CASE THE FD IS NOT DONE WITH THIS SCRIPT FOR A GIVEN MESON
#####  THAT IS, IF useScriptFDpaths and doFeedDownMes are both set to 0 for that meson

declare dirppDzero="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2014Nov9paperProp/InputPlotsDzero/pp/Singlev2Envelope"
declare dirpPbDzero="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2014Nov9paperProp/InputPlotsDzero/pPb/pPb/TotalEnvelope"

declare dirppDplus="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2014Nov9paperProp/InputPlotsDplus/pp"
declare dirpPbDplus="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2014Nov9paperProp/InputPlotsDplus/pPb"

declare dirppDstar="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2014Nov9paperProp/InputPlotsDstar/pp/CorrectedPlots"
declare dirpPbDstar="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2014Nov9paperProp/InputPlotsDstar/pPb/CorrectedPlots"

declare -a dirppInput=( "$dirppDzero" "$dirppDplus" "$dirppDstar" )
declare -a dirpPbInput=( "$dirpPbDzero" "$dirpPbDplus" "$dirpPbDstar" )

#################################################################################################

#################################################################################################
##### THE PARAMETERS BELOW CAN BE CHANGED TO SPECIFY WHETHER YOU WANT TO REFLECT OR NOT, REBIN OR NOT

## default rebin = -1 (no rebin)
declare rebin=-1  
declare reflect=1   #0 : do not reflect, 1 reflect
declare averageOpt=0  #0 = weighted, 1=arithmetic (BE AWARE THAT SETTINGS FOR FINAL STYLE PLOTS MIGHT NOT BE OK FOR ARITHMETIC AVERAGE AND THEY MIGHT USE PATHS EXPECTED FOR WEIGHTING AV)
declare -a includev2=( 0 0 1 0 ) 
##### Subtraction of b-origin modulation (MC closure) - Set 0 in pp2010, pPb2013; set 1 in pPb2016. IF 1, REMEMBER TO SET THE UNCERTAINTY UPWARDS!
declare subtrMCclos=1  #0: don't do anything for MC closure structures (use old approach), #1 apply correction for MC closure b->D structures
declare purityOpt=30 #0: simple rescale (flat), 10: multiply by prim/inclusive histogram, 20: multiply by template fit on prim/inclusive histogram, 30: multiply by moving mean of prim/inclusive histogram
declare doNewPedestalVars=1 #USE ALWAYS 1 UNLESS FOR REPRODUCING 2010+2013 PAPER RESULTS; #0: use paper pedestal variations; #1: use new variations (less transverse regions, added baseline +/- isigma_stat)

###############################################################################
############ YOU CAN CHOOSE TO DO ONLY SOME STEPS           ###################
############  IN CASE SOME WERE ALREADY DONE WITH THIS VERY SAME SCRIPT #######
declare -i cpCode=0 # THIS WILL MAKE THE COMMITTED MACRO TO BE COPIED AND USED IN THE HFCJlocalCodeDir DIRECTORY, WHICH IS EXPORTED. IF YOU WANT TO MODIFY CODE YOU CAN RUN WITH THIS SET TO 1 THE FIRST TIME AND THEN SET IT TO 0. 
declare useoldFilenames=0 #USE 0 IF YOU TAKE USE COMMON EXTRACTION FRAMEWORK - USE 1 ONLY TO REPRODUCE pp2010/pPb2013 PAPER RESULTS! 
declare -ai useScriptFDpaths=( 1 1 1 ) #DO NOT CHANGE IT, UNLESS YOU KNOW!!  THIS IS USEFUL IN CASE YOU DO NOT WANT TO RECOMPUTE THE FD BUT USE THE PATHS SET BY THE SCRIPT FOR THE FILES COMING FROM THE FD SUBTRACTION
declare -i doFeedDownGlob=0
declare -ia doFeedDownMes=( 1 1 1 ) ## Dzero, Dstar, Dplus values
declare doInitAndReflStep=0 ## NOTE THAT THIS STEP IS NECESSARY ALSO IN CASE THE PLOTS DO NOT HAVE TO BE REFLECTED. THE "reflect" PARAMETER ABOVE IS WHAT DETERMINED WHETHER OR NOT THE PLOTS WILL BE REFLECTED INSIDE THIS STEP. THIS PARAMETER IS JUST A SWITCH: YOU CAN SET IT TO 0 IN CASE YOU ALREADY DID IT AND YOU DO NOT WANT TO REPEAT IT
declare doAverage=0
declare dofit=0
declare doDrawFitFigure=0
declare dofitMC=1
declare dofitawayside=1
#Paper 2010&2013
declare doNicePlot=0
declare doCompareMesons=0 #
declare doFitResultComparisonPPpPb=0 #For paper 2010 & 2013
declare doFitResultComparisonPPtoMC=0 #For paper 2010 & 2013
declare doFitResultComparisonPPbtoMC=0 #For paper 2010 & 2013
declare doFitResultComparisonPPtoMCawayside=0 #For paper 2010 & 2013
declare doFitResultComparisonPPbtoMCawayside=0 #For paper 2010 & 2013
declare doFitResultComparisonPPtoPPbtoMCPP=0 #For paper 2010 & 2013
declare doCompareWithMC=0 #For paper 2010 & 2013
declare doComparepppPb=0 #For paper 2010 & 2013
declare doOldPlots=0 #For paper 2010 & 2013
#pPb 2016 preliminaries
declare doCompareWithMC2016=0 #For pPb 2016
declare doFitResultComparisonPPbtoMC2016=0 #For paper 2016
declare doFitResultComparisonPPbtoMCawayside2016=0 #For paper 2016
declare doComparepppPb2016=0 #For pPb 2016
declare doFitResultComparisonPPpPb2016=0 #For pPb 2016
declare doProduceSQMplots=0 #For pPb 2016
#pp 2017 5 TeV preliminaries
declare doNiceStyleSinglePanelAverages=1 #For pp 2017 @ 5 TeV prels
declare doCompareWithOtherSystems=1 #For pp 2017 @ 5 TeV prels
declare doCompareWithMCPP2017=1 #For pp 2017 @ 5 TeV prels
declare doFitResultComparisonPPtoMC2017=1 #For pp 2017 @ 5 TeV prels
declare doFitResultComparisonPPtoMC2017awayside=1 #For pp 2017 @ 5 TeV prels

#MODIFY DIR CONTAINING OTHER SYSTEMS OUTPUT FOR SYSTEMS COMPARISON IN PP 2017 MODE
declare othersystdir="$PWD/OtherSystemsForCfr"

### MODIFY THE PARAMETERS BELOW ONLY IF YOU WANT TO RUN ONLY OVER FEW KINE CASES, OTHERWISE DO NOT TOUCH THEM
declare -i firstcollsyst=3   #default 0
declare -i lastcollsyst=3   #default 1 
declare -i firstmeson=0     #default 0
declare -i lastmeson=2      #default 2

declare -i firstassocbin=0  #default 0
declare -i lastassocbin=6   #default 2
declare -i firsttrigbin=0   #default 0, note that 0=2-3, 1=3-5, 2=5-8, 3=8-16, 4=16-24
declare -i lasttrigbin=4    #default 2

declare -i centralitybin=0  #default 0, use always 0 for cent-integrated cases! then, 1 is 0-20%, 2 is 20-60%, 3 is 60-100%

#################################################################################################
## BELOW THIS LINE NOTHING SHOULD BE CHANGED
#################################################################################################

declare rebinDir="StdRebin"
declare -a meson=("Dzero" "Dstar" "Dplus")
declare baseFile="FDsub"
declare endFilepp="_v2D0.00_v2had0.00.root"

declare -a pttrig=("2to3" "3to5" "5to8" "8to16" "16to24")
declare -a ptassoc=("0.3to1.0" "0.3to99.0" "1.0to99.0" "2.0to99.0" "3.0to99.0" "1.0to2.0" "2.0to3.0")
declare -a collsystdir=("pp" "pPb" "pPb2016" "pp2017")


declare -i collsyst=${firstcollsyst}
declare -i imeson=${firstmeson}
declare -i itrigbin=${firsttrigbin}
declare -i iassocbin=${firstassocbin}

##############################################################################################
##### CREATE A DIRECTORY WHERE PAPER FIGURES WILL BE LINKED, JUST FOR CONVENIENCE #####
###############################################################################################

cd ${baseStartingDir}

ln -s ${ALICE_PHYSICS} ${ALICE_PHYSICS}/../../git #needed for AliBuild
ln -s ${ALICE_PHYSICS} ${ALICE_PHYSICS}/../src

############# GO TO MACRO PATH AND CP CODE ############
if [ ${cpCode} = 1 ]; then
    mkdir -p ${HFCJlocalCodeDir}
    cd ${HFCJlocalCodeDir}
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoAverages.sh .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoSubtractFD.sh .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DrawPlotInSubtractFD.sh .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoFit.sh .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/RunFeedown_pPb_Dplus.C .	
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/RunFeedown_pPb_Dstar.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/RunFeedown_pPb_Dzero.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/RunFeedown_pp_Dplus.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/RunFeedown_pp_Dstar.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/RunFeedown_pp_Dzero.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/MakeAverageDhCorrel.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoPlotInSingleCanvas.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoPlotInSingleCanvasNoSpaces.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoPlotCompare1GeVpPb.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoPlotCompare1GeVpp.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoPlotComparedot3to1pPb.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoPlotComparedot3to1pp.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoComparison_pPbVsMC.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoComparison_ppVsMC.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoComparison_ppVspPb.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/FitSystematicsAverage_pp.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/FitSystematicsAverage_pp2017.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/FitSystematicsAverage_pPb.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/FitSystematicsAverage_pPb2016.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/FitPlots.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/SubtractFD.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoNiceSpecieComparisonPlot.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/CompareFitResults.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/CompareFitResults_DiffppEnergies.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/CompareFitResultspp2017.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/CompareFitResultspPb2016.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/CompareFitResults_ppVspPb_1.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/CompareFitResults_ppVspPb_2.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoNiceFitPlots.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoComparison_ppVspPb.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoComparison_ppVspPballPanels.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoComparison_ppVspPballPanels2016.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoComparison_ppVsMCallPanelsNew.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoComparison_Distributions_DiffppEnergies.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoComparison_Distributions_ppVspPb_1.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoComparison_Distributions_ppVspPb_2.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoComparison_pp2017VsMCallPanels.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoComparison_pPb2016VsMCallPanelsNew.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoComparison_ppVsMC.C . 
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoComparison_pPbVsMC.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoComparison_pPbVsMC2016.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DrawSinglePanel_pp2017_CorrelationDistributions_1.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DrawSinglePanel_pp2017_CorrelationDistributions_2.C .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/DoFitMC.sh .
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/CheckDiff.sh .
else 
    cd ${HFCJlocalCodeDir}
    cp ${ALICE_PHYSICS}/../../git/PWGHF/correlationHF/macros/CheckDiff.sh .
    ${HFCJlocalCodeDir}/CheckDiff.sh &>outCodeDiff.log
fi

rm ${ALICE_PHYSICS}/../../git #was needed for AliBuild

#############  CREATE BASIC DIRECTORY TREE #############
declare baseDir=${baseStartingDir}
cd ${baseDir}
if [ ${reflect} = 1 ]; then
    baseDir=${baseDir}/ReflectedPlots
elif [ ${reflect} = 0 ]; then
    baseDir=${baseDir}/2PiRange
fi

if [ ${rebin} -ne -1 ]; then
    rebinDir=Rebin${rebin}
fi
baseDir=${baseDir}/${rebinDir}
mkdir -p ${baseDir}
cd ${baseDir}
echo "BASE DIR WILL BE ${baseDir}"
mkdir AllPlots
cd AllPlots
mkdir Averages
cd Averages
mkdir FitResults

declare baseDirFD=${baseStartingDir}/FDsubtr
if [ ${doFeedDownGlob} = 1 ]; then
    mkdir -p ${baseDirFD}
fi


collsyst=${firstcollsyst}
imeson=${firstmeson}
itrigbin=${firsttrigbin}
iassocbin=${firstassocbin}

while [ ${collsyst} -le ${lastcollsyst} ]; do
    cd ${baseDir}
    mkdir ${collsystdir[${collsyst}]}
    imeson=${firstmeson}
    while [ ${imeson} -le ${lastmeson} ]; do
	cd ${baseDir}/${collsystdir[${collsyst}]}
	mkdir ${meson[${imeson}]}
	if [ ${doFeedDownMes[${imeson}]} = 1 ];then
	    cd ${meson[${imeson}]}
	    mkdir -p ${baseDirFD}
	    ## FORCE USAGE OF SCRIPT FD PATHS EVEN IF THIS IS DONE FOR A SINGLE COLLISION SYSTEM... TAKE CARE!!!
	    useScriptFDpaths[${imeson}]=1
	fi
	if [ ${useScriptFDpaths[${imeson}]} = 1 ];then
	    if [ ${collsyst} = 0 ]; then
		echo "changing dir with input pp files for meson ${imeson} to that expected from this script after FD subtr"
		dirppInput[${imeson}]="${baseDirFD}/${collsystdir[${collsyst}]}/${meson[imeson]}/Final_Plots_pp/Singlev2Envelope"
		echo "New dir will be: ${dirppInput[${imeson}]}"
	    elif [ ${collsyst} = 1 ]; then
		echo "changing dir with input pPb files for meson ${imeson} to that expected from this script after FD subtr"
		dirpPbInput[${imeson}]="${baseDirFD}/${collsystdir[${collsyst}]}/${meson[imeson]}/Final_Plots_pPb/TotalEnvelope"
		#dirpPbInput[${imeson}]="${baseDirFD}/${collsystdir[${collsyst}]}/${meson[imeson]}/Final_Plots_pPb/Singlev2Envelope"
		echo "New dir will be: ${dirpPbInput[${imeson}]}"
        elif [ ${collsyst} = 2 ]; then
        echo "changing dir with input pPb2016 files for meson ${imeson} to that expected from this script after FD subtr"
        dirpPbInput[${imeson}]="${baseDirFD}/${collsystdir[${collsyst}]}/${meson[imeson]}/Final_Plots_pPb/TotalEnvelope"
        #dirpPbInput[${imeson}]="${baseDirFD}/${collsystdir[${collsyst}]}/${meson[imeson]}/Final_Plots_pPb/Singlev2Envelope"
        echo "New dir will be: ${dirpPbInput[${imeson}]}"
        elif [ ${collsyst} = 3 ]; then
        echo "changing dir with input pp files for meson ${imeson} to that expected from this script after FD subtr"
        dirppInput[${imeson}]="${baseDirFD}/${collsystdir[${collsyst}]}/${meson[imeson]}/Final_Plots_pp/Singlev2Envelope"
        echo "New dir will be: ${dirppInput[${imeson}]}"        
	    fi
	fi
	imeson=${imeson}+1
    done
    collsyst=${collsyst}+1
done

collsyst=${firstcollsyst}
imeson=${firstmeson}
itrigbin=${firsttrigbin}
iassocbin=${firstassocbin}



####### Feed Down subtraction ####
### IF YOU WANT TO USE LOCAL CODE, CHECK THE DoSubtractFD.sh script (you have 3 possibilities)
if [ $doFeedDownGlob = 1 ]; then
    while [ ${collsyst} -le ${lastcollsyst} ]; do
	imeson=${firstmeson}
	while [ ${imeson} -le ${lastmeson} ]; do
	    if [ ${doFeedDownMes[${imeson}]} = 0 ]; then
		imeson=${imeson}+1
		continue
	    fi
	    echo "ProducePlotChain: subtracting FD for meson $imeson in coll system $collsyst"
	    mkdir -p ${baseDirFD}/${collsystdir[${collsyst}]}/${meson[${imeson}]}
	    cd ${baseDirFD}/${collsystdir[${collsyst}]}/${meson[${imeson}]}
	    if [ ${imeson} = 0 ]; then
		$HFCJlocalCodeDir/DoSubtractFD.sh ${collsyst} ${imeson} ${dirDzeroNotFDsubt[${collsyst}]}/${fpromptfileDzero[${collsyst}]} ${templateDir[${collsyst}]} ${dirDzeroNotFDsubt[${collsyst}]} ${filerootDzero[${collsyst}]} 3 ${templateDirSystemSuffix[${collsyst}]} ${subtrMCclos} ${useoldFilenames} ${centralitybin} ${purityOpt} ${puritytemplateDir[${collsyst}]}
	    elif [ ${imeson} = 1 ]; then
		$HFCJlocalCodeDir/DoSubtractFD.sh ${collsyst} ${imeson} ${dirDstarNotFDsubt[${collsyst}]}/${fpromptfileDstar[${collsyst}]} ${templateDir[${collsyst}]} ${dirDstarNotFDsubt[${collsyst}]} ${filerootDstar[${collsyst}]} 3 ${templateDirSystemSuffix[${collsyst}]} ${subtrMCclos} ${useoldFilenames} ${centralitybin} ${purityOpt} ${puritytemplateDir[${collsyst}]}
	    elif [ ${imeson} = 2 ]; then
		$HFCJlocalCodeDir/DoSubtractFD.sh ${collsyst} ${imeson} ${dirDplusNotFDsubt[${collsyst}]}/${fpromptfileDplus[${collsyst}]} ${templateDir[${collsyst}]} ${dirDplusNotFDsubt[${collsyst}]} ${filerootDplus[${collsyst}]} 3 ${templateDirSystemSuffix[${collsyst}]} ${subtrMCclos} ${useoldFilenames} ${centralitybin} ${purityOpt} ${puritytemplateDir[${collsyst}]}
	    fi
	    imeson=${imeson}+1
	done
	collsyst=${collsyst}+1
    done
    echo "Feed down done"
fi


######### NOW DO REFLECTION ##########
if [ $doInitAndReflStep = 1 ]; then
    collsyst=${firstcollsyst}
    while [ ${collsyst} -le ${lastcollsyst} ]; do
	imeson=${firstmeson}
	while [ ${imeson} -le ${lastmeson} ]; do
	    itrigbin=${firsttrigbin}
	    while [ ${itrigbin} -le ${lasttrigbin} ]; do
		cd ${baseDir}/${collsystdir[${collsyst}]}/${meson[${imeson}]}
		echo "Curent path: ${baseDir}/${collsystdir[${collsyst}]}/${meson[${imeson}]}"
		iassocbin=${firstassocbin}
		while [ ${iassocbin} -le ${lastassocbin} ]; do
		    echo "DOING REFLECTION Case: system ${collsystdir[collsyst]}, meson ${meson[$imeson]}, pt trig ${pttrig[$itrigbin]}, pt assoc ${ptassoc[$iassocbin]}"
		if [ ${collsyst} = 0 ]; then 
			if [ ${imeson} = 0 ]; then
			#echo "ciao"
			    $HFCJlocalCodeDir/DrawPlotInSubtractFD.sh ${dirppInput[${imeson}]}/${collsystdir[${collsyst}]}_${baseFile}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assoc${ptassoc[${iassocbin}]}${endFilepp} ${imeson} ${collsyst} ${itrigbin} ${iassocbin} ${reflect} ${rebin} 1 ${centralitybin}
			elif [ ${imeson} = 1 ]; then
			#echo "ciao"
			    $HFCJlocalCodeDir/DrawPlotInSubtractFD.sh ${dirppInput[${imeson}]}/${collsystdir[${collsyst}]}_${baseFile}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assoc${ptassoc[${iassocbin}]}${endFilepp} ${imeson} ${collsyst} ${itrigbin} ${iassocbin} ${reflect} ${rebin} 1 ${centralitybin}
			elif [ ${imeson} = 2 ]; then
			#echo "ciao"
			    $HFCJlocalCodeDir/DrawPlotInSubtractFD.sh ${dirppInput[${imeson}]}/${collsystdir[${collsyst}]}_${baseFile}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assoc${ptassoc[${iassocbin}]}${endFilepp} ${imeson} ${collsyst} ${itrigbin} ${iassocbin} ${reflect} ${rebin} 1 ${centralitybin}
			fi
		elif [ ${collsyst} = 1 ]; then 
			if [ ${imeson} = 0 ]; then
			#echo "ciao"
			    $HFCJlocalCodeDir/DrawPlotInSubtractFD.sh ${dirpPbInput[${imeson}]}/${collsystdir[${collsyst}]}_${baseFile}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assoc${ptassoc[${iassocbin}]}.root ${imeson} ${collsyst} ${itrigbin} ${iassocbin} ${reflect} ${rebin} 1 ${centralitybin}
			elif [ ${imeson} = 1 ]; then
			#echo "ciao"
			    $HFCJlocalCodeDir/DrawPlotInSubtractFD.sh ${dirpPbInput[${imeson}]}/${collsystdir[${collsyst}]}_${baseFile}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assoc${ptassoc[${iassocbin}]}.root ${imeson} ${collsyst} ${itrigbin} ${iassocbin} ${reflect} ${rebin} 1  ${centralitybin}
			elif [ ${imeson} = 2 ]; then
			#echo "ciao"
			    $HFCJlocalCodeDir/DrawPlotInSubtractFD.sh ${dirpPbInput[${imeson}]}/${collsystdir[${collsyst}]}_${baseFile}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assoc${ptassoc[${iassocbin}]}.root ${imeson} ${collsyst} ${itrigbin} ${iassocbin} ${reflect} ${rebin} 1 ${centralitybin}
			fi
        elif [ ${collsyst} = 2 ]; then 
            if [ ${imeson} = 0 ]; then
            #echo "ciao"
                $HFCJlocalCodeDir/DrawPlotInSubtractFD.sh ${dirpPbInput[${imeson}]}/${collsystdir[1]}_${baseFile}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assoc${ptassoc[${iassocbin}]}.root ${imeson} ${collsyst} ${itrigbin} ${iassocbin} ${reflect} ${rebin} 1 ${centralitybin}
            elif [ ${imeson} = 1 ]; then
            #echo "ciao"
                $HFCJlocalCodeDir/DrawPlotInSubtractFD.sh ${dirpPbInput[${imeson}]}/${collsystdir[1]}_${baseFile}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assoc${ptassoc[${iassocbin}]}.root ${imeson} ${collsyst} ${itrigbin} ${iassocbin} ${reflect} ${rebin} 1 ${centralitybin}
            elif [ ${imeson} = 2 ]; then
            #echo "ciao"
                $HFCJlocalCodeDir/DrawPlotInSubtractFD.sh ${dirpPbInput[${imeson}]}/${collsystdir[1]}_${baseFile}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assoc${ptassoc[${iassocbin}]}.root ${imeson} ${collsyst} ${itrigbin} ${iassocbin} ${reflect} ${rebin} 1 ${centralitybin}
            fi
        elif [ ${collsyst} = 3 ]; then 
            if [ ${imeson} = 0 ]; then
            #echo "ciao"
                $HFCJlocalCodeDir/DrawPlotInSubtractFD.sh ${dirppInput[${imeson}]}/${collsystdir[0]}_${baseFile}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assoc${ptassoc[${iassocbin}]}${endFilepp} ${imeson} ${collsyst} ${itrigbin} ${iassocbin} ${reflect} ${rebin} 1 ${centralitybin}
            elif [ ${imeson} = 1 ]; then
            #echo "ciao"
                $HFCJlocalCodeDir/DrawPlotInSubtractFD.sh ${dirppInput[${imeson}]}/${collsystdir[0]}_${baseFile}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assoc${ptassoc[${iassocbin}]}${endFilepp} ${imeson} ${collsyst} ${itrigbin} ${iassocbin} ${reflect} ${rebin} 1 ${centralitybin}
            elif [ ${imeson} = 2 ]; then
            #echo "ciao"
                $HFCJlocalCodeDir/DrawPlotInSubtractFD.sh ${dirppInput[${imeson}]}/${collsystdir[0]}_${baseFile}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assoc${ptassoc[${iassocbin}]}${endFilepp} ${imeson} ${collsyst} ${itrigbin} ${iassocbin} ${reflect} ${rebin} 1 ${centralitybin}
            fi
		    fi
        #now make symbolic link of newly created plots
        if [ ${collsyst} = 2 ]; then     
		    if [ -e CanvaAndVariedHisto${collsystdir[1]}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assocPt${ptassoc[${iassocbin}]}.root ]; then
			ln -s $PWD/CanvaAndVariedHisto${collsystdir[1]}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assocPt${ptassoc[${iassocbin}]}.root ${baseDir}/AllPlots/
			ln -s $PWD/CanvaAndVariedHisto${collsystdir[1]}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assocPt${ptassoc[${iassocbin}]}.eps ${baseDir}/AllPlots/
		    ln -s $PWD/CanvaAndVariedHisto${collsystdir[1]}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assocPt${ptassoc[${iassocbin}]}.png ${baseDir}/AllPlots/
		    fi
        elif [ ${collsyst} = 3 ]; then     
            if [ -e CanvaAndVariedHisto${collsystdir[0]}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assocPt${ptassoc[${iassocbin}]}.root ]; then
            ln -s $PWD/CanvaAndVariedHisto${collsystdir[0]}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assocPt${ptassoc[${iassocbin}]}.root ${baseDir}/AllPlots/
            ln -s $PWD/CanvaAndVariedHisto${collsystdir[0]}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assocPt${ptassoc[${iassocbin}]}.eps ${baseDir}/AllPlots/
            ln -s $PWD/CanvaAndVariedHisto${collsystdir[0]}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assocPt${ptassoc[${iassocbin}]}.png ${baseDir}/AllPlots/
            fi
        else 
            if [ -e CanvaAndVariedHisto${collsystdir[$collsyst]}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assocPt${ptassoc[${iassocbin}]}.root ]; then
            ln -s $PWD/CanvaAndVariedHisto${collsystdir[$collsyst]}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assocPt${ptassoc[${iassocbin}]}.root ${baseDir}/AllPlots/
            ln -s $PWD/CanvaAndVariedHisto${collsystdir[$collsyst]}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assocPt${ptassoc[${iassocbin}]}.eps ${baseDir}/AllPlots/
            ln -s $PWD/CanvaAndVariedHisto${collsystdir[$collsyst]}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assocPt${ptassoc[${iassocbin}]}.png ${baseDir}/AllPlots/
            fi
        fi

		iassocbin=${iassocbin}+1
		done
	    itrigbin=${itrigbin}+1
	    done
	    imeson=${imeson}+1
	done
	collsyst=${collsyst}+1
    done
fi
 

collsyst=${firstcollsyst}
imeson=${firstmeson}
itrigbin=${firsttrigbin}
iassocbin=${firstassocbin}

######### NOW DO AVERAGES #########
cd ${baseDir}/AllPlots/
if [ ${doAverage} = 1 ]; then
    cd Averages
    while [ ${collsyst} -le ${lastcollsyst} ]; do
	itrigbin=${firsttrigbin}
	while [ ${itrigbin} -le ${lasttrigbin} ]; do
	    iassocbin=${firstassocbin}
	    while [ ${iassocbin} -le ${lastassocbin} ]; do
		echo "Averaging for ${collsystdir[${collsyst}]} collisions,  pt(D) bin ${itrigbin}, pt(assoc) bin  ${iassocbin} "
		$HFCJlocalCodeDir/DoAverages.sh ${collsyst} ${itrigbin} ${iassocbin} ${averageOpt} ${reflect} 5 ${baseDir}/AllPlots/ ${baseDir}/AllPlots/ ${baseDir}/AllPlots/ 1 ${centralitybin}
		iassocbin=${iassocbin}+1
		done
	    itrigbin=${itrigbin}+1
	    done	    
	collsyst=${collsyst}+1
    done	       
fi

collsyst=${firstcollsyst}
imeson=${firstmeson}
itrigbin=${firsttrigbin}
iassocbin=${firstassocbin}


###### NOW PRODUCE CANVASES WITH NICE STYLE AND THOSE COMPARING MESONS ###########
#### NOTE THAT THE POSSIBILITY TO SKIP CASES (e.g. pp collisions or p-Pb collisions) it is not yet implemented

if [ ${doNicePlot} = 1 ]; then
    cd ${baseDir}/AllPlots
    mkdir NiceStylePlots
    cd NiceStylePlots
    root -b <<EOF &> PlotNiceStyle.log
.L ${HFCJlocalCodeDir}/DoPlotInSingleCanvas.C
SetInputDirectory("${baseDir}/AllPlots")
DoAllPlotInSingleCanvasStandardPaths(${averageOpt})
.q
EOF

    cd ${baseDir}/AllPlots/NiceStylePlots
    root -b <<EOF &> PlotNiceStyleNoSpace.log
.L ${HFCJlocalCodeDir}/DoPlotInSingleCanvasNoSpaces.C
DoPlotInSingleCanvasNoSpaces()
.q
EOF
    
    cd ${baseDir}/AllPlots/NiceStylePlots/Output_Plots/WeightedAverageDzeroDstarDplus/
    root -b <<EOF &> PlotNiceStyleNoSpaceMergePPandpPb.log
.L ${HFCJlocalCodeDir}/DoPlotInSingleCanvasNoSpaces.C
MergePPandPPbInSingleCanvas("${baseDir}/AllPlots/NiceStylePlots/Output_Plots/WeightedAverageDzeroDstarDplus/CanvasNoSpaces_WeightedAverageDzeroDstarDplus_pp.root","${baseDir}/AllPlots/NiceStylePlots/Output_Plots/WeightedAverageDzeroDstarDplus/CanvasNoSpaces_WeightedAverageDzeroDstarDplus_pPb.root")
EOF
fi

if [ ${doCompareMesons} = 1 ];then
    cd ${baseDir}/AllPlots
    mkdir CompareMesons
    cd CompareMesons
    root -b <<EOF &> CompareMesonsPtAssoc1to99pPb.log
.L ${HFCJlocalCodeDir}/DoPlotCompare1GeVpPb.C
SetInputDirectory("${baseDir}/AllPlots")
DoPlotCompare1GeVpPb()
.q
EOF
   
 root -b <<EOF &> CompareMesonsPtAssoc1to99pp.log
.L ${HFCJlocalCodeDir}/DoPlotCompare1GeVpp.C
SetInputDirectory("${baseDir}/AllPlots")
DoPlotCompareAbove1pp()
.q
EOF
 
    root -b <<EOF &> CompareMesonsPtAssoc03to1pp.log
.L ${HFCJlocalCodeDir}/DoPlotComparedot3to1pp.C
SetInputDirectory("${baseDir}/AllPlots")
DoPlotComparedot3to1pp()
.q
EOF
    
    root -b <<EOF &> CompareMesonsPtAssoc03to1pPb.log
.L ${HFCJlocalCodeDir}/DoPlotComparedot3to1pPb.C
SetInputDirectory("${baseDir}/AllPlots")
DoPlotComparedot3to1pPb()
.q
EOF
    

    root -b <<EOF &> CompareMesonsNiceStyle.log
.L ${HFCJlocalCodeDir}/DoNiceSpecieComparisonPlot.C
SetInputDirectory("${baseDir}/AllPlots/CompareMesons/Output_SngCav_Comparison")
DoNiceSpecieComparisonPlot("${pttrig[1]}","${ptassoc[2]}","${collsystdir[0]}","${pttrig[2]}","${ptassoc[2]}","${collsystdir[1]}")
.q
EOF
fi

######## NOW FIT DISTRIBUTIONS ############
if [ ${dofitMC} = 1 ]; then    
    echo "Produce Plot Chain: fit MC distributions"
    while [ ${collsyst} -le ${lastcollsyst} ]; do
	mkdir -p ${templateDir[${collsyst}]}/FitResults/
	cd ${templateDir[${collsyst}]}/FitResults
	for (( mccase=0; mccase<${Nmccase[${collsyst}]}; mccase++ ))
	do 
	    if [ ${collsyst} = 0  ]; then
		$HFCJlocalCodeDir/DoFitMC.sh ${collsyst} ${isreflectedMC[mccase]} ${averageOpt} ${templateDir[${collsyst}]}  ${templateDir[${collsyst}]}/FitResults/ ${collsystdir[${collsyst}]}${templRootNamepp[$mccase]} ${dofitawayside} ${betaTemplDir[${collsyst}]}
	    elif [ ${collsyst} = 1  ]; then
		$HFCJlocalCodeDir/DoFitMC.sh ${collsyst} ${isreflectedMC[mccase]} ${averageOpt} ${templateDir[${collsyst}]}  ${templateDir[${collsyst}]}/FitResults/ ${collsystdir[${collsyst}]}${templRootNamepPb[$mccase]} ${dofitawayside} ${betaTemplDir[${collsyst}]}
        elif [ ${collsyst} = 2  ]; then
        $HFCJlocalCodeDir/DoFitMC.sh ${collsyst} ${isreflectedMC[mccase]} ${averageOpt} ${templateDir[${collsyst}]}  ${templateDir[${collsyst}]}/FitResults/ ${collsystdir[1]}${templRootNamepPb[$mccase]} ${dofitawayside} ${betaTemplDir[${collsyst}]}
	    elif [ ${collsyst} = 3  ]; then
        $HFCJlocalCodeDir/DoFitMC.sh ${collsyst} ${isreflectedMC[mccase]} ${averageOpt} ${templateDir[${collsyst}]}  ${templateDir[${collsyst}]}/FitResults/ ${collsystdir[0]}${templRootNamepp[$mccase]} ${dofitawayside} ${betaTemplDir[${collsyst}]}
        fi
	done
	collsyst=${collsyst}+1
    done
fi
collsyst=${firstcollsyst}

if [ ${dofit} = 1 ]; then
    echo "Produce Plot Chain: fit data distributions"
    cd ${baseDir}/AllPlots/Averages/FitResults/    
    while [ ${collsyst} -le ${lastcollsyst} ]; do
	$HFCJlocalCodeDir/DoFit.sh ${collsyst} ${reflect} ${averageOpt} ${baseDir}/AllPlots/Averages/  ${baseDir}/AllPlots/Averages/FitResults/ ${includev2[${collsyst}]} ${dofitawayside} ${doNewPedestalVars} ${betaTemplDir[${collsyst}]}
	collsyst=${collsyst}+1
    done
fi
collsyst=${firstcollsyst}

if [ ${doDrawFitFigure} = 1 ]; then
  cd ${baseDir}/AllPlots/Averages/FitResults
  mkdir NiceStylePlots
  cd NiceStylePlots
  root -b <<EOF &> NiceStylePlots.log
  .L ${HFCJlocalCodeDir}/DoNiceFitPlots.C
  SetInputDirectory("${baseDir}/AllPlots/Averages/FitResults")
  DoNiceFitPlots()
  .q
EOF
fi

if [ ${doFitResultComparisonPPpPb} = 1 ];then
    mkdir -p ${baseDir}/AllPlots/Averages/FitResults/ComparisonPPtoPPb
    cd ${baseDir}/AllPlots/Averages/FitResults/ComparisonPPtoPPb
    root -b <<EOF &> CompareFitResultsUniqueCanvasPPtoPPb.log
.L ${HFCJlocalCodeDir}/CompareFitResults.C
SetDirectoryFitResultPP("${baseDir}/AllPlots/Averages/FitResults/")
SetDirectoryFitResultPPb("${baseDir}/AllPlots/Averages/FitResults/")
CompareFitResultsPPtoPPbUniqueCanvas();
EOF

    root -b <<EOF &> CompareFitResultsPPtoPPb.log
.L ${HFCJlocalCodeDir}/CompareFitResults.C
SetDirectoryFitResultPP("${baseDir}/AllPlots/Averages/FitResults/")
SetDirectoryFitResultPPb("${baseDir}/AllPlots/Averages/FitResults/")
CompareFitResultsPPtoPPb()
EOF
fi


if [ ${doFitResultComparisonPPtoMC} = 1 ];then
    collsyst=0
    mkdir -p ${baseDir}/AllPlots/Averages/FitResults/ComparisonPPtoMC
    cd ${baseDir}/AllPlots/Averages/FitResults/ComparisonPPtoMC
    root -b <<EOF &> CompareFitResultsPPtoMC.log
.L ${HFCJlocalCodeDir}/CompareFitResults.C
SetDirectoryFitResultPP("${baseDir}/AllPlots/Averages/FitResults/")
SetDirectoryFitResultsMCpp("${templateDir[${collsyst}]}/FitResults/")
CompareFitResultsPPDataToMC()
EOF

    root -b <<EOF &> CompareFitResultsPPtoMCuniqueCanvas.log
.L ${HFCJlocalCodeDir}/CompareFitResults.C
SetDirectoryFitResultPP("${baseDir}/AllPlots/Averages/FitResults/")
SetDirectoryFitResultsMCpp("${templateDir[${collsyst}]}/FitResults/")
IncludeModel(0,${mccasePP[0]})
IncludeModel(1,${mccasePP[1]})
IncludeModel(2,${mccasePP[2]})
IncludeModel(3,${mccasePP[3]})
IncludeModel(4,${mccasePP[4]})
IncludeModel(5,${mccasePP[5]})
IncludeModel(6,${mccasePP[6]})
IncludeModel(7,${mccasePP[7]})
SetDrawSystMC(kFALSE)
CompareFitResultsPPtoMCUniqueCanvas()
EOF
fi

if [ ${doFitResultComparisonPPtoPPbtoMCPP} = 1 ]; then
 root -b <<EOF &> CompareFitResultsPPtoPPbtoPPMCuniqueCanvas.log
.L ${HFCJlocalCodeDir}/CompareFitResults.C
SetDirectoryFitResultPP("${baseDir}/AllPlots/Averages/FitResults/")
SetDirectoryFitResultsMCpp("${templateDir[${collsyst}]}/FitResults/")
SetDirectoryFitResultPPb("${baseDir}/AllPlots/Averages/FitResults/")
IncludeModel(0,${mccasePP[0]}||${mccasePPb[0]})
IncludeModel(1,${mccasePP[1]}||${mccasePPb[1]})
IncludeModel(2,${mccasePP[2]}||${mccasePPb[2]})
IncludeModel(3,${mccasePP[3]}||${mccasePPb[3]})
IncludeModel(4,${mccasePP[4]}||${mccasePPb[4]})
IncludeModel(5,${mccasePP[5]}||${mccasePPb[5]})
IncludeModel(6,${mccasePP[6]}||${mccasePPb[6]})
IncludeModel(7,${mccasePP[7]}||${mccasePPb[7]})
CompareFitResultsPPtoPpbAndMCUniqueCanvas();
EOF

collsyst=${firstcollsyst}
fi

if [ ${doFitResultComparisonPPbtoMC} = 1 ];then
    collsyst=1
    mkdir -p ${baseDir}/AllPlots/Averages/FitResults/ComparisonPPbtoMC
    cd ${baseDir}/AllPlots/Averages/FitResults/ComparisonPPbtoMC
    root -b <<EOF &> CompareFitResultsPPbtoMC.log
.L ${HFCJlocalCodeDir}/CompareFitResults.C
SetDirectoryFitResultPPb("${baseDir}/AllPlots/Averages/FitResults/")
SetDirectoryFitResultsMCpPb("${templateDir[${collsyst}]}/FitResults/")
IncludePowheg(kFALSE)
CompareFitResultsPPbDataToMC()
EOF

    root -b <<EOF &> CompareFitResultsPPbtoMCuniqueCanvas.log
.L ${HFCJlocalCodeDir}/CompareFitResults.C
SetDirectoryFitResultPPb("${baseDir}/AllPlots/Averages/FitResults/")
SetDirectoryFitResultsMCpPb("${templateDir[${collsyst}]}/FitResults/")
IncludeModel(0,${mccasePPb[0]})
IncludeModel(1,${mccasePPb[1]})
IncludeModel(2,${mccasePPb[2]})
IncludeModel(3,${mccasePPb[3]})
IncludeModel(4,${mccasePPb[4]})
IncludeModel(5,${mccasePPb[5]})
IncludeModel(6,${mccasePPb[6]})
IncludeModel(7,${mccasePPb[7]})
SetDrawSystMC(kFALSE)
SetMinPtDisplayData((Double_t)${minptdisplaypPb})
SetMinPtDisplayMC((Double_t)${minptdisplaypPb})
CompareFitResultsPPbtoMCUniqueCanvas()
EOF

collsyst=${firstcollsyst}
fi


##### do away side if requested

if [ ${doFitResultComparisonPPtoMCawayside} = 1 ];then
    collsyst=0
    mkdir -p ${baseDir}/AllPlots/Averages/FitResults/ComparisonPPtoMC
    cd ${baseDir}/AllPlots/Averages/FitResults/ComparisonPPtoMC

    root -b <<EOF &> CompareFitResultsPPtoMCuniqueCanvasAwaySide.log
.L ${HFCJlocalCodeDir}/CompareFitResults.C
SetDirectoryFitResultPP("${baseDir}/AllPlots/Averages/FitResults/")
SetDirectoryFitResultsMCpp("${templateDir[${collsyst}]}/FitResults/")
IncludeModel(0,${mccasePP[0]})
IncludeModel(1,${mccasePP[1]})
IncludeModel(2,${mccasePP[2]})
IncludeModel(3,${mccasePP[3]})
IncludeModel(4,${mccasePP[4]})
IncludeModel(5,${mccasePP[5]})
IncludeModel(6,${mccasePP[6]})
IncludeModel(7,${mccasePP[7]})
SetDrawSystMC(kFALSE)
CompareFitResultsPPtoMCUniqueCanvasAwaySide()
EOF

collsyst=${firstcollsyst}
fi

if [ ${doFitResultComparisonPPbtoMCawayside} = 1 ];then
    collsyst=1
    mkdir -p ${baseDir}/AllPlots/Averages/FitResults/ComparisonPPbtoMC
    cd ${baseDir}/AllPlots/Averages/FitResults/ComparisonPPbtoMC

    root -b <<EOF &> CompareFitResultsPPbtoMCuniqueCanvasAS.log
.L ${HFCJlocalCodeDir}/CompareFitResults.C
SetDirectoryFitResultPPb("${baseDir}/AllPlots/Averages/FitResults/")
SetDirectoryFitResultsMCpPb("${templateDir[${collsyst}]}/FitResults/")
IncludeModel(0,${mccasePPb[0]})
IncludeModel(1,${mccasePPb[1]})
IncludeModel(2,${mccasePPb[2]})
IncludeModel(3,${mccasePPb[3]})
IncludeModel(4,${mccasePPb[4]})
IncludeModel(5,${mccasePPb[5]})
IncludeModel(6,${mccasePPb[6]})
IncludeModel(7,${mccasePPb[7]})
SetDrawSystMC(kFALSE)
SetMinPtDisplayData((Double_t)${minptdisplaypPb})
SetMinPtDisplayMC((Double_t)${minptdisplaypPb})
CompareFitResultsPPbtoMCUniqueCanvasAwaySide()
EOF

collsyst=${firstcollsyst}
fi

######## NOW COMPARE DATA AND MC ################# 
if [ $doCompareWithMC = 1 ]; then
    cd ${baseDir}/AllPlots/Averages
    mkdir ComparisonToModels
    cd ComparisonToModels
    ### START FROM PP ###
    if [ ${doOldPlots} = 1 ]; then
	echo "Running DoComparison_ppVsMC_PaperProposal.C for ptassoc 0.3to1"
	root -b <<EOF &> CompareToModelsPP03to1.log
.L ${HFCJlocalCodeDir}/DoComparison_ppVsMC.C
Printf("Macro loaded")
SetInputDataDirectory("${baseDir}/AllPlots/Averages")
SetInputTemplateDirectory("${templateDir[0]}")
SetFitPlotMacroPath("${HFCJlocalCodeDir}")
SetReflectTemplate($reflect)
SetIsDataReflected($reflect)
SetSkip3to5pPb(kFALSE)
SetFDtemplateSystemString("${templateDirSystemSuffix[${collsyst}]}");
SetBaselineDirectory("${baseDir}/AllPlots/Averages/FitResults")
Printf("READY TO GO")
SetAverageMode($averageOpt)
DoComparison_ppVsMCTEST("0.3to1.0")
EOF
	echo "Running DoComparison_ppVsMC_PaperProposal.C for ptassoc 0.3to99"
	root -b <<EOF &> CompareToModelsPP03to99.log
.L ${HFCJlocalCodeDir}/DoComparison_ppVsMC.C
Printf("Macro loaded")
SetInputDataDirectory("${baseDir}/AllPlots/Averages")
SetInputTemplateDirectory("${templateDir[0]}")
SetFitPlotMacroPath("${HFCJlocalCodeDir}")
SetReflectTemplate($reflect)
SetIsDataReflected($reflect)
SetSkip3to5pPb(kFALSE)
SetFDtemplateSystemString("${templateDirSystemSuffix[${collsyst}]}");
SetBaselineDirectory("${baseDir}/AllPlots/Averages/FitResults")
Printf("READY TO GO")
SetAverageMode($averageOpt)
DoComparison_ppVsMCTEST("0.3to99.0")
EOF
	echo "Running DoComparison_ppVsMC_PaperProposal.C for ptassoc 1.0to99.0"
	root -b <<EOF &> CompareToModelsPP1to99.log
.L ${HFCJlocalCodeDir}/DoComparison_ppVsMC.C
Printf("Macro loaded")
SetInputDataDirectory("${baseDir}/AllPlots/Averages")
SetInputTemplateDirectory("${templateDir[0]}")
SetFitPlotMacroPath("${HFCJlocalCodeDir}")
SetReflectTemplate($reflect)
SetIsDataReflected($reflect)
SetSkip3to5pPb(kFALSE)
SetBaselineDirectory("${baseDir}/AllPlots/Averages/FitResults")
Printf("READY TO GO")
SetAverageMode($averageOpt)
DoComparison_ppVsMCTEST("1.0to99.0")
EOF
	
    ### NOW PPb ###
	echo "Running DoComparison_pPbVsMC_PaperProposal.C for ptassoc 0.3to1.0"
	root -b <<EOF &> CompareToModelsPPB03to1.log
.L ${HFCJlocalCodeDir}/DoComparison_pPbVsMC.C
Printf("Macro loaded")
SetInputDataDirectory("${baseDir}/AllPlots/Averages")
SetSkip3to5pPb(kFALSE)
SetInputTemplateDirectory("${templateDir[1]}")
SetFitPlotMacroPath("${HFCJlocalCodeDir}")
SetReflectTemplate($reflect)
SetIsDataReflected($reflect)
SetBaselineDirectory("${baseDir}/AllPlots/Averages/FitResults")
SetFDtemplateSystemString("${templateDirSystemSuffix[${collsyst}]}");
Printf("READY TO GO")
SetAverageMode($averageOpt)
DoComparison_pPbVsMCTEST("0.3to1.0")
EOF
	
	echo "Running DoComparison_pPbVsMC_PaperProposal.C for ptassoc 0.3to99.0"
	root -b <<EOF &> CompareToModelsPPB03to99.log
.L ${HFCJlocalCodeDir}/DoComparison_pPbVsMC.C
Printf("Macro loaded")
SetInputDataDirectory("${baseDir}/AllPlots/Averages")
SetSkip3to5pPb(kFALSE)
SetInputTemplateDirectory("${templateDir[1]}")
SetReflectTemplate($reflect)
SetIsDataReflected($reflect)
SetBaselineDirectory("${baseDir}/AllPlots/Averages/FitResults")
SetFitPlotMacroPath("${HFCJlocalCodeDir}")
SetFDtemplateSystemString("${templateDirSystemSuffix[${collsyst}]}");
Printf("READY TO GO")
SetAverageMode($averageOpt)
DoComparison_pPbVsMCTEST("0.3to99.0")
EOF
	
	echo "Running DoComparison_pPbVsMC_PaperProposal.C for ptassoc 1.0to99.0"
	root -b <<EOF &> CompareToModelsPPB1to99.log
.L ${HFCJlocalCodeDir}/DoComparison_pPbVsMC.C
Printf("Macro loaded")
SetInputDataDirectory("${baseDir}/AllPlots/Averages")
SetSkip3to5pPb(kFALSE)
SetInputTemplateDirectory("${templateDir[1]}")
SetReflectTemplate($reflect)
SetIsDataReflected($reflect)
SetFDtemplateSystemString("${templateDirSystemSuffix[${collsyst}]}");
SetBaselineDirectory("${baseDir}/AllPlots/Averages/FitResults")
SetFitPlotMacroPath("${HFCJlocalCodeDir}")
Printf("READY TO GO")
SetAverageMode($averageOpt)
DoComparison_pPbVsMCTEST("1.0to99.0")
EOF
    fi    
    
### now make figur with all panels in same picture for pp
    echo "Running DoComparison_ppVsMCallPanelsNew.C"
    root -b <<EOF &> ComparePPtoMCuniqueCanvas.log
.L ${HFCJlocalCodeDir}/DoComparison_ppVsMCallPanelsNew.C
SetInputDataDirectory("${baseDir}/AllPlots/Averages")
SetInputTemplateDirectory("${templateDir[0]}")
SetFitPlotMacroPath("${HFCJlocalCodeDir}")
SetReflectTemplate($reflect)
SetIsDataReflected($reflect)
SetBaselineDirectory("${baseDir}/AllPlots/Averages/FitResults")
SetAverageMode($averageOpt)
SetSplitMClegendInTwoPanels(kTRUE)
SetIncludeAllMCmodels(kTRUE)
SetIncludeEPOS(kTRUE)
DoComparison_ppVsMCallPanels()
EOF 
fi

#Comparison of pp 2010 and pPb 2013
if [ $doComparepppPb = 1 ]; then
    cd ${baseDir}/AllPlots/Averages
    mkdir ComparisonPPtoPPB
    cd ComparisonPPtoPPB
### NB THE FOLLOWING ARE SET TO KFALSE ON PURPOSE (THEY PRODUCE OBSOLETE PLOTS)
    if [ $doComparepppPb = 0 ]; then
	echo "Running DoComparison_ppVspPb_PaperProposal.C for ptassoc 1to99.0"
	root -b <<EOF &> ComparePPtoPPB1to99.log
.L ${HFCJlocalCodeDir}/DoComparison_ppVspPb.C
Printf("Macro loaded")
SetInputDataDirectory("${baseDir}/AllPlots/Averages")
SetBaselineDirectory("${baseDir}/AllPlots/Averages/FitResults")
SetSkip3to5pPb(kTRUE)
SetIsReflected($reflect)
Printf("READY TO GO")
SetAverageMode($averageOpt)
DoComparison_ppVspPbTEST("1.0to99.0")
.q
EOF
    
    echo "Running DoComparison_ppVspPb_PaperProposal.C for ptassoc0.3to1.0"
    root -b <<EOF &> ComparePPtoPPB03to1.log
.L ${HFCJlocalCodeDir}/DoComparison_ppVspPb.C
Printf("Macro loaded")
SetInputDataDirectory("${baseDir}/AllPlots/Averages")
SetBaselineDirectory("${baseDir}/AllPlots/Averages/FitResults")
SetSkip3to5pPb(kTRUE)
SetIsReflected($reflect)
Printf("READY TO GO")
SetAverageMode($averageOpt)
DoComparison_ppVspPbTEST("0.3to1.0")
.q
EOF
    echo "Running DoComparison_ppVspPb_PaperProposal.C for ptassoc0.3to99.0"
    root -b <<EOF &> ComparePPtoPPB0.3to99.log
.L ${HFCJlocalCodeDir}/DoComparison_ppVspPb.C
Printf("Macro loaded")
SetInputDataDirectory("${baseDir}/AllPlots/Averages")
SetBaselineDirectory("${baseDir}/AllPlots/Averages/FitResults")
SetSkip3to5pPb(kTRUE)
SetIsReflected($reflect)
Printf("READY TO GO")
SetAverageMode($averageOpt)
DoComparison_ppVspPbTEST("0.3to99.0")
.q
EOF
    echo "Produce Unique CANVAS"
    root -b <<EOF &> ProduceUniqueCanvas.log
.L ${HFCJlocalCodeDir}/DoComparison_ppVspPballPanels.C
Printf("Macro loaded")
SetInputDataDirectory("${baseDir}/AllPlots/Averages")
SetBaselineDirectory("${baseDir}/AllPlots/Averages/FitResults")
SetSkip3to5pPb(kTRUE)
SetIsReflected($reflect)
Printf("READY TO GO")
SetAverageMode($averageOpt)
DoComparison_ppVspPballPanels()
.q
EOF
fi

#########################
### 2016 Comparisons! ###
#########################

######## NOW COMPARE DATA AND MC FOR pPb 2016 ################# 
if [ $doCompareWithMC2016 = 1 ]; then
    cd ${baseDir}/AllPlots/Averages
    mkdir ComparisonToModels
    cd ComparisonToModels
    
### make figure with all panels
    echo "Running DoComparison_pPb2016VsMCallPanelsNew.C"
    root -b <<EOF &> ComparePb2016toMCuniqueCanvas.log
.L ${HFCJlocalCodeDir}/DoComparison_pPb2016VsMCallPanelsNew.C
SetInputDataDirectory("${baseDir}/AllPlots/Averages")
SetInputTemplateDirectory("${templateDir[2]}")
SetFitPlotMacroPath("${HFCJlocalCodeDir}")
SetReflectTemplate($reflect)
SetIsDataReflected($reflect)
SetBaselineDirectory("${baseDir}/AllPlots/Averages/FitResults")
SetAverageMode($averageOpt)
SetSplitMClegendInTwoPanels(kTRUE)
//SetIncludeAllMCmodels(kFALSE)
//SetIncludeEPOS(kFALSE)
DoComparison_pPb2016VsMCallPanels()
DoComparison_pPb2016VsMCSinglePanel()
EOF
fi

### compare observables from fits now
if [ ${doFitResultComparisonPPbtoMC2016} = 1 ];then
    collsyst=1
    mkdir -p ${baseDir}/AllPlots/Averages/FitResults/ComparisonPPbtoMC
    cd ${baseDir}/AllPlots/Averages/FitResults/ComparisonPPbtoMC
    echo "Running CompareFitResultspPb2016.C"
    root -b <<EOF &> CompareFitResultsPPbtoMC.log
.L ${HFCJlocalCodeDir}/CompareFitResultspPb2016.C
SetDirectoryFitResultPPb("${baseDir}/AllPlots/Averages/FitResults/")
SetDirectoryFitResultsMCpPb("${templateDir[${collsyst}]}/FitResults/")
IncludeModel(0,${mccasePPb[0]})
IncludeModel(1,${mccasePPb[1]})
IncludeModel(2,${mccasePPb[2]})
IncludeModel(3,${mccasePPb[3]})
IncludeModel(4,${mccasePPb[4]})
IncludeModel(5,${mccasePPb[5]})
IncludeModel(6,${mccasePPb[6]}) 
IncludeModel(7,${mccasePPb[7]})
CompareFitResultsPPbDataToMC()
CompareFitResultsPPbtoMCUniqueCanvas()
CompareFitResultsPPbtoMCUniqueCanvas_2()
EOF

    root -b <<EOF &> CompareFitResultsPPbtoMCuniqueCanvas.log
.L ${HFCJlocalCodeDir}/CompareFitResultspPb2016.C
SetDirectoryFitResultPPb("${baseDir}/AllPlots/Averages/FitResults/")
SetDirectoryFitResultsMCpPb("${templateDir[${collsyst}]}/FitResults/")
IncludeModel(0,${mccasePPb[0]})
IncludeModel(1,${mccasePPb[1]})
IncludeModel(2,${mccasePPb[2]})
IncludeModel(3,${mccasePPb[3]})
IncludeModel(4,${mccasePPb[4]})
IncludeModel(5,${mccasePPb[5]})
IncludeModel(6,${mccasePPb[6]}) 
IncludeModel(7,${mccasePPb[7]})
SetDrawSystMC(kTRUE)
SetMinPtDisplayData((Double_t)${minptdisplaypPb})
SetMinPtDisplayMC((Double_t)${minptdisplaypPb})
CompareFitResultsPPbtoMCUniqueCanvasAwaySide()
CompareFitResultsPPbtoMCUniqueCanvasAwaySide_2()
EOF

fi

##### do away side if requested

if [ ${doFitResultComparisonPPbtoMCawayside2016} = 1 ];then
    collsyst=1
    mkdir -p ${baseDir}/AllPlots/Averages/FitResults/ComparisonPPbtoMC
    cd ${baseDir}/AllPlots/Averages/FitResults/ComparisonPPbtoMC

    root -b <<EOF &> CompareFitResultsPPbtoMCuniqueCanvasAS.log
.L ${HFCJlocalCodeDir}/CompareFitResultspPb2016.C
SetDirectoryFitResultPPb("${baseDir}/AllPlots/Averages/FitResults/")
SetDirectoryFitResultsMCpPb("${templateDir[${collsyst}]}/FitResults/")
IncludeModel(0,${mccasePPb[0]})
IncludeModel(1,${mccasePPb[1]})
IncludeModel(2,${mccasePPb[2]})
IncludeModel(3,${mccasePPb[3]})
IncludeModel(4,${mccasePPb[4]})
IncludeModel(5,${mccasePPb[5]})
IncludeModel(6,${mccasePPb[6]})
IncludeModel(7,${mccasePPb[7]})
SetDrawSystMC(kFALSE)
CompareFitResultsPPbtoMCUniqueCanvasAwaySide()
EOF
fi

##### COMPARE DATA pPb 2016 vs pp 2010
if [ $doComparepppPb2016 = 1 ]; then
    cd ${baseDir}/AllPlots/Averages
    mkdir ComparisonPPtoPPB
    cd ComparisonPPtoPPB
    echo "Produce Unique CANVAS"
    root -b <<EOF &> ProduceUniqueCanvas.log
.L ${HFCJlocalCodeDir}/DoComparison_ppVspPballPanels2016.C
Printf("Macro loaded")
SetInputDataDirectory("${baseDir}/AllPlots/Averages")
SetBaselineDirectory("${baseDir}/AllPlots/Averages/FitResults")
SetSkip3to5pPb(kFALSE)
SetIsReflected($reflect)
Printf("READY TO GO")
SetAverageMode($averageOpt)
DoComparison_ppVspPballPanels2016()
.q
EOF
fi

#Comparison of pp 2010 and pPb 2016 For Fit Observables (NS)
if [ ${doFitResultComparisonPPpPb2016} = 1 ];then
    mkdir -p ${baseDir}/AllPlots/Averages/FitResults/ComparisonPPtoPPb
    cd ${baseDir}/AllPlots/Averages/FitResults/ComparisonPPtoPPb
    root -b <<EOF &> CompareFitResultsUniqueCanvasPPtoPPb.log
.L ${HFCJlocalCodeDir}/CompareFitResultspPb2016.C
SetDirectoryFitResultPP("${baseDir}/AllPlots/Averages/FitResults/")
SetDirectoryFitResultPPb("${baseDir}/AllPlots/Averages/FitResults/")
SetRunOn2016(kTRUE);
CompareFitResultsPPtoPPbUniqueCanvas();
EOF

    root -b <<EOF &> CompareFitResultsPPtoPPb.log
.L ${HFCJlocalCodeDir}/CompareFitResultspPb2016.C
SetDirectoryFitResultPP("${baseDir}/AllPlots/Averages/FitResults/")
SetDirectoryFitResultPPb("${baseDir}/AllPlots/Averages/FitResults/")
CompareFitResultsPPtoPPb()
EOF
fi

##### Produce SQM plots
if [ $doProduceSQMplots = 1 ]; then
    cd ${HFCJlocalCodeDir}
    echo "Produce SQM plots"
    root -b <<EOF &> SQM_Plots.log
.L ${HFCJlocalCodeDir}/Restyle_pPb2016_Prel_Plots.C
Restyle_pPb2016_Prel_Plots()
.q
EOF
fi

################################
#### pp 2017 specific plots ####
################################

## Nice style, single panel average plots
if [ $doNiceStyleSinglePanelAverages = 1 ]; then
  cd ${baseDir}/AllPlots/Averages
  mkdir MultiPanelAverages
  cd MultiPanelAverages

  echo "Doing multi-panel averages"
    root -b <<EOF &> DrawSinglePanel_pp2017_CorrelationDistributions_1.log
.L ${HFCJlocalCodeDir}/DrawSinglePanel_pp2017_CorrelationDistributions_1.C
SetBasicDir("${baseStartingDir}/..");
SetInputDataDirectory("${othersystdir}/pp_5TeV")
SetBaselineDirectory("${othersystdir}/pp_5TeV/FitResults/Trends_pp")
DrawSinglePanel_CorrelationDistributions_1()
EOF
    root -b <<EOF &> DrawSinglePanel_pp2017_CorrelationDistributions_2.log
.L ${HFCJlocalCodeDir}/DrawSinglePanel_pp2017_CorrelationDistributions_2.C
SetBasicDir("${baseStartingDir}/..");
SetInputDataDirectory("${othersystdir}/pp_5TeV")
SetBaselineDirectory("${othersystdir}/pp_5TeV/FitResults/Trends_pp")
DrawSinglePanel_CorrelationDistributions_2()
EOF
fi

## Comparison with other systems (first pPb 5 TeV, then pp at 7, 13 TeV)
if [ $doCompareWithOtherSystems = 1 ]; then
  
  if [ ! -d ${othersystdir}/pp_5TeV ]; then
    ln -s ${baseDir}/AllPlots/Averages/ ${othersystdir}
    mv ${othersystdir}/Averages ${othersystdir}/pp_5TeV
  fi
  cd ${baseDir}/AllPlots/Averages
  mkdir ComparisonTopPb
  cd ComparisonTopPb

  ##pp vs pPb
  echo "Running DoComparison_Distributions_ppVspPb_1.C and _2.C"
    root -b <<EOF &> DoComparison_Distributions_ppVspPb_1.log
.L ${HFCJlocalCodeDir}/DoComparison_Distributions_ppVspPb_1.C
SetBasicDir("${baseStartingDir}/..");
SetInputDataDirectory("${othersystdir}/pp_5TeV","${othersystdir}/pPb_5TeV")
SetBaselineDirectory("${othersystdir}/pp_5TeV/FitResults/Trends_pp","${othersystdir}/pPb_5TeV/FitResults/Trends_pPb")
DoComparison_Distributions()
EOF
    root -b <<EOF &> DoComparison_Distributions_ppVspPb_2.log
.L ${HFCJlocalCodeDir}/DoComparison_Distributions_ppVspPb_2.C
SetBasicDir("${baseStartingDir}/..");
SetInputDataDirectory("${othersystdir}/pp_5TeV","${othersystdir}/pPb_5TeV")
SetBaselineDirectory("${othersystdir}/pp_5TeV/FitResults/Trends_pp","${othersystdir}/pPb_5TeV/FitResults/Trends_pPb")
DoComparison_Distributions()
EOF

#ADDED...!!!
  echo "Running CompareFitResults_ppVspPb vs AssocPt track.C"
    root -b <<EOF &> CompareFitResults_ppVspPb_vsPtAss.log
.L ${HFCJlocalCodeDir}/CompareFitResults_ppVspPb_vsPtAss.C++
SetBasicDir("${baseStartingDir}/..");
SetDirectoryFitResult("${othersystdir}/pp_5TeV/FitResults/Trends_pp","${othersystdir}/pPb_5TeV/FitResults/Trends_pPb")
CompareFitResults_vsPtAss_UniqueCanvas()
CompareFitResults_vsPtAss_UniqueCanvas_AwaySide()
EOF
#...ADDED!!!

  echo "Running CompareFitResults_ppVspPb_1.C and _2.C"
    root -b <<EOF &> CompareFitResults_ppVspPb_1.log
.L ${HFCJlocalCodeDir}/CompareFitResults_ppVspPb_1.C++
SetBasicDir("${baseStartingDir}/..");
SetDirectoryFitResult("${othersystdir}/pp_5TeV/FitResults/Trends_pp","${othersystdir}/pPb_5TeV/FitResults/Trends_pPb")
CompareFitResults_UniqueCanvas()
CompareFitResults_UniqueCanvas_AwaySide()
EOF
    root -b <<EOF &> CompareFitResults_ppVspPb_2.log
.L ${HFCJlocalCodeDir}/CompareFitResults_ppVspPb_2.C++
SetBasicDir("${baseStartingDir}/..");
SetDirectoryFitResult("${othersystdir}/pp_5TeV/FitResults/Trends_pp","${othersystdir}/pPb_5TeV/FitResults/Trends_pPb")
CompareFitResults_UniqueCanvas()
CompareFitResults_UniqueCanvas_AwaySide()
EOF

  ##pp vs other pp
  cd ${baseDir}/AllPlots/Averages
  mkdir ComparisonToOtherpp
  cd ComparisonToOtherpp
  
  echo "Running DoComparison_Distributions_DiffppEnergies.C"
    root -b <<EOF &> DoComparison_Distributions_DiffppEnergies.log
.L ${HFCJlocalCodeDir}/DoComparison_Distributions_DiffppEnergies.C
SetBasicDir("${baseStartingDir}/..");
SetInputDataDirectory("${othersystdir}/pp_5TeV","${othersystdir}/pp_7TeV","${othersystdir}/pp_13TeV")
SetBaselineDirectory("${othersystdir}/pp_5TeV/FitResults/Trends_pp","${othersystdir}/pp_7TeV/FitResults/Trends_pp","${othersystdir}/pp_13TeV/FitResults/Trends_pp")
DoComparison_Distributions()
EOF

  echo "Running CompareFitResults_DiffppEnergies.C"
    root -b <<EOF &> CompareFitResults_DiffppEnergies.log
.L ${HFCJlocalCodeDir}/CompareFitResults_DiffppEnergies.C
SetBasicDir("${baseStartingDir}/..");
SetDirectoryFitResult("${othersystdir}/pp_5TeV/FitResults/Trends_pp","${othersystdir}/pp_7TeV/FitResults/Trends_pp","${othersystdir}/pp_13TeV/FitResults/Trends_pp")
CompareFitResults_UniqueCanvas()
CompareFitResults_UniqueCanvas_AwaySide()
EOF

fi

######## NOW COMPARE DATA AND MC FOR pp 2017 ################# 
if [ $doCompareWithMCPP2017 = 1 ]; then
    cd ${baseDir}/AllPlots/Averages
    mkdir ComparisonToModels
    cd ComparisonToModels
    
### make figure with all panels
    echo "Running DoComparison_pp2017VsMCallPanels.C"
    root -b <<EOF &> ComparePP2017toMCuniqueCanvas.log
.L ${HFCJlocalCodeDir}/DoComparison_pp2017VsMCallPanels.C
SetInputDataDirectory("${baseDir}/AllPlots/Averages")
SetInputTemplateDirectory("${templateDir[${collsyst}]}")
SetFitPlotMacroPath("${HFCJlocalCodeDir}")
SetReflectTemplate($reflect)
SetIsDataReflected($reflect)
SetBaselineDirectory("${baseDir}/AllPlots/Averages/FitResults")
SetAverageMode($averageOpt)
SetSplitMClegendInTwoPanels(kTRUE)
SetIncludePerugia2011()
SetIncludePYTHIA8()
SetIncludeHERWIG()
SetIncludePOWHEG()
SetIncludePOWHEG_LO()
SetIncludeEPOS()
DoComparison_pp2017VsMCallPanels()
DoComparison_pp2017VsMCSinglePanel()
EOF
fi


### compare observables from fits now
if [ ${doFitResultComparisonPPtoMC2017} = 1 ];then
    collsyst=3
    mkdir -p ${baseDir}/AllPlots/Averages/ComparisonToModels
    cd ${baseDir}/AllPlots/Averages/ComparisonToModels
    echo "Running CompareFitResultspp2017.C"

    root -b <<EOF &> CompareFitResultspp2017toMC_uniquecanvas.log
.L ${HFCJlocalCodeDir}/CompareFitResultspp2017.C
SetDirectoryFitResultPP("${baseDir}/AllPlots/Averages/FitResults/")
SetDirectoryFitResultsMCPP("${templateDir[${collsyst}]}/FitResults/")
IncludeModel(0,${mccasePP[0]})
IncludeModel(1,${mccasePP[1]})
IncludeModel(2,${mccasePP[2]})
IncludeModel(3,${mccasePP[3]})
IncludeModel(4,${mccasePP[4]})
IncludeModel(5,${mccasePP[5]})
IncludeModel(6,${mccasePP[6]}) 
IncludeModel(7,${mccasePP[7]})
CompareFitResultsPPDataToMC()
CompareFitResultsPPtoMCUniqueCanvas()
CompareFitResultsPPtoMCUniqueCanvas_2()
CompareFitResults_Ratios_NS_1()
CompareFitResults_Ratios_NS_2()
EOF

#ADDED...!!!
  echo "Running CompareFitResultspp2017toMC vs AssocPt track.C"
    root -b <<EOF &> CompareFitResultspp2017_vsPtAss.log
.L ${HFCJlocalCodeDir}/CompareFitResultspp2017_vsPtAss.C
SetDirectoryFitResultPP("${baseDir}/AllPlots/Averages/FitResults/")
SetDirectoryFitResultsMCPP("${templateDir[${collsyst}]}/FitResults/")
IncludeModel(0,${mccasePP[0]})
IncludeModel(1,${mccasePP[1]})
IncludeModel(2,${mccasePP[2]})
IncludeModel(3,${mccasePP[3]})
IncludeModel(4,${mccasePP[4]})
IncludeModel(5,${mccasePP[5]})
IncludeModel(6,${mccasePP[6]}) 
IncludeModel(7,${mccasePP[7]})
SetDrawSystMC(kTRUE)
CompareFitResultsPPtoMCUniqueCanvas_vsPtAss()
CompareFitResults_Ratios_NS_vsPtAss()
EOF
#...ADDED!!!

fi

##### do away side if requested

if [ ${doFitResultComparisonPPtoMC2017awayside} = 1 ];then
    collsyst=3
    mkdir -p ${baseDir}/AllPlots/Averages/ComparisonToModels
    cd ${baseDir}/AllPlots/Averages/ComparisonToModels
    echo "Running CompareFitResultspp2017.C for away side"

    root -b <<EOF &> CompareFitResultspp2017toMC_AS_uniquecanvas.log
.L ${HFCJlocalCodeDir}/CompareFitResultspp2017.C
SetDirectoryFitResultPP("${baseDir}/AllPlots/Averages/FitResults/")
SetDirectoryFitResultsMCPP("${templateDir[${collsyst}]}/FitResults/")
IncludeModel(0,${mccasePP[0]})
IncludeModel(1,${mccasePP[1]})
IncludeModel(2,${mccasePP[2]})
IncludeModel(3,${mccasePP[3]})
IncludeModel(4,${mccasePP[4]})
IncludeModel(5,${mccasePP[5]})
IncludeModel(6,${mccasePP[6]}) 
IncludeModel(7,${mccasePP[7]})
SetDrawSystMC(kTRUE)
CompareFitResultsPPtoMCUniqueCanvasAwaySide()
CompareFitResultsPPtoMCUniqueCanvasAwaySide_2()
CompareFitResults_Ratios_AS_1()
CompareFitResults_Ratios_AS_2()
EOF

#ADDED...!!!
  echo "Running CompareFitResultspp2017toMC vs AssocPt track.C"
    root -b <<EOF &> CompareFitResultspp2017_AS_vsPtAss.log
.L ${HFCJlocalCodeDir}/CompareFitResultspp2017_vsPtAss.C
SetDirectoryFitResultPP("${baseDir}/AllPlots/Averages/FitResults/")
SetDirectoryFitResultsMCPP("${templateDir[${collsyst}]}/FitResults/")
IncludeModel(0,${mccasePP[0]})
IncludeModel(1,${mccasePP[1]})
IncludeModel(2,${mccasePP[2]})
IncludeModel(3,${mccasePP[3]})
IncludeModel(4,${mccasePP[4]})
IncludeModel(5,${mccasePP[5]})
IncludeModel(6,${mccasePP[6]}) 
IncludeModel(7,${mccasePP[7]})
SetDrawSystMC(kTRUE)
CompareFitResultsPPtoMCUniqueCanvasAwaySide_vsPtAss()
CompareFitResults_Ratios_AS_vsPtAss()
EOF
#...ADDED!!!
fi

rm ${ALICE_PHYSICS}/../src #was needed by §AliBuild


    
