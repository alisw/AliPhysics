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


declare baseStartingDir="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015May21testFabioInputs/TEST"
#declare macrosDir="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015March13/Macros"
export HFCJlocalCodeDir="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015May21testFabioInputs/TEST"
# "/Users/administrator/soft/alisoft/aliphysics/master/src/PWGHF/correlationHF/macros"
declare templateDirPP="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015May21testFabioInputs/Templates_pp"
declare templateDirPPb="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015May21testFabioInputs/Templates_pPb"
declare -a templateDirSystemSuffix=( "none" "none" ) #### THIS IS KEPT JUST FOR BACKWARD COMPATIBILITY WITH OLD TEMPLATES! NO NEED TO TOUCH IT UNLESS YOU WANT TO USE OLD TEMPLATES
declare -a templateDir=( "$templateDirPP" "$templateDirPPb" )


########## THE FOLLOWING DIRECTORIES SHOULD CONTAIN THE RESULTS BEFORE FD SUBTRACTION #####
declare dirppDzeroNotFDsubt="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015May21testFabioInputs/Dzero_pp"
declare dirpPbDzeroNotFDsubt="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015May21testFabioInputs/Dzero_pPb"
declare -a dirDzeroNotFDsubt=( "$dirppDzeroNotFDsubt" "$dirpPbDzeroNotFDsubt" )
declare -a fpromptfileDzero=( "HFPtSpectrum_pp.root" "HFPtSpectrum_DrawFpromptVsRaaElossHypoCombined.root" )
declare -a filerootDzero=( "1D_Signal_WithEMCorr_Normal_Charg_OriginSuper_Integrated_Bins" "1D_Signal_WithEMCorr_Normal_Charg_OriginSuper_Integrated_Bins" )

declare dirppDstarNotFDsubt="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015May21testFabioInputs/Dstar_pp"
declare dirpPbDstarNotFDsubt="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015May21testFabioInputs/Dstar_pPb"
declare -a dirDstarNotFDsubt=( "$dirppDstarNotFDsubt" "$dirpPbDstarNotFDsubt" )
declare -a fpromptfileDstar=( "outputkfcB6_23mb.root" "fPromptWithBeautyRpA.root" )
declare -a filerootDstar=( "FinalDphiCorrelationsCanvas_" "FinalDphiCorrelationsCanvas_" )

declare dirppDplusNotFDsubt="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015May21testFabioInputs/Dplus_pp"
declare dirpPbDplusNotFDsubt="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015May21testFabioInputs/Dplus_pPb"
declare -a dirDplusNotFDsubt=( "$dirppDplusNotFDsubt" "$dirpPbDplusNotFDsubt" )
declare -a fpromptfileDplus=( "HFPtSpectrum_ppDplus_kfc_kpp7.root" "DrawFpromptVsRaaElossHypoCombined.root" )
declare -a filerootDplus=( "1D_pp_DplusHCorr_" "1D_pPb_DplusHCorr_" )


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
declare averageOpt=0  #0 = weighted, 1=arithmetic

###############################################################################
############ YOU CAN CHOOSE TO DO ONLY SOME STEPS           ###################
############  IN CASE SOME WERE ALREADY DONE WITH THIS VERY SAME SCRIPT #######
declare -i cpCode=1 # THIS WILL MAKE THE COMMITTED MACRO TO BE COPIED AND USED IN THE HFCJlocalCodeDir DIRECTORY, WHICH IS EXPORTED. IF YOU WANT TO MODIFY CODE YOU CAN RUN WITH THIS SET TO 1 THE FIRST TIME AND THEN SET IT TO 0. 
declare -ai useScriptFDpaths=( 1 1 1 ) #THIS IS USEFUL IN CASE YOU DO NOT WANT TO RECOMPUTE THE FD BUT USE THE PATHS SET BY THE SCRIPT FOR THE FILES COMING FROM THE FD SUBTRACTION
declare -i doFeedDownGlob=1
declare -ia doFeedDownMes=( 1 1 1 ) ## Dzero, Dstar, Dplus values
declare doInitAndReflStep=1  ## NOTE THAT THIS STEP IS NECESSARY ALSO IN CASE THE PLOTS DO NOT HAVE TO BE REFLECTED. THE "reflect" PARAMETER ABOVE IS WHAT DETERMINED WHETHER OR NOT THE PLOTS WILL BE REFLECTED INSIDE THIS STEP. THIS PARAMETER IS JUST A SWITCH: YOU CAN SET IT TO 0 IN CASE YOU ALREADY DID IT AND YOU DO NOT WANT TO REPEAT IT
declare doAverage=1
declare doNicePlot=1
declare doCompareMesons=1
declare dofit=1
declare doCompareWithMC=1
declare doComparepppPb=1

### MODIFY THE PARAMETERS BELOW ONLY IF YOU WANT TO RUN ONLY OVER FEW KINE CASES, OTHERWISE DO NOT TOUCH THEM
declare -i firstcollsyst=0   #default 0
declare -i lastcollsyst=1    #default 1 
declare -i firstmeson=0     #default 0
declare -i lastmeson=2      #default 2
declare -i firstassocbin=0  #default 0
declare -i lastassocbin=2   #default 2

declare -i firsttrigbin=0   #default 0
declare -i lasttrigbin=2    #default 2



#################################################################################################
## BELOW THIS LINE NOTHING SHOULD BE CHANGED
#################################################################################################

declare rebinDir="StdRebin"
declare -a meson=("Dzero" "Dstar" "Dplus")
declare baseFile="FDsub"
declare endFilepp="_v2D0.00_v2had0.00.root"

declare -a pttrig=("3to5" "5to8" "8to16")
declare -a ptassoc=("0.3to1.0" "0.3to99.0" "1.0to99.0")
declare -a collsystdir=("pp" "pPb")


declare -i collsyst=${firstcollsyst}
declare -i imeson=${firstmeson}
declare -i itrigbin=${firsttrigbin}
declare -i iassocbin=${firstassocbin}



############# GO TO MACRO PATH AND CP CODE ############
if [ ${cpCode} = 1 ]; then
    mkdir -p ${HFCJlocalCodeDir}
    cd ${HFCJlocalCodeDir}
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoAverages.sh .
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoSubtractFD.sh .
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DrawPlotInSubtractFD.sh .
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoFit.sh .
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/RunFeedown_pPb_Dplus.C .	
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/RunFeedown_pPb_Dstar.C .
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/RunFeedown_pPb_Dzero.C .
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/RunFeedown_pp_Dplus.C .
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/RunFeedown_pp_Dstar.C .
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/RunFeedown_pp_Dzero.C .
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/MakeAverageDhCorrel.C .
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoPlotInSingleCanvas.C .
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoPlotCompare1GeVpPb.C .
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoPlotCompare1GeVpp.C .
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoPlotComparedot3to1pPb.C .
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoPlotComparedot3to1pp.C .
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoComparison_pPbVsMC.C .
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoComparison_ppVsMC.C .
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoComparison_ppVspPb.C .
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/FitSystematicsAverage_pp.C .
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/FitSystematicsAverage_pPb.C .
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/FitPlots.C .
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/SubtractFD.C .
    cp ${ALICE_PHYSICS}/../src/PWGHF/correlationHF/macros/DoNiceSpecieComparisonPlot.C .
fi

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
		echo "New dir will be: ${dirpPbInput[${imeson}]}"
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
		$HFCJlocalCodeDir/DoSubtractFD.sh ${collsyst} ${imeson} ${dirDzeroNotFDsubt[${collsyst}]}/${fpromptfileDzero[${collsyst}]} ${templateDir[${collsyst}]} ${dirDzeroNotFDsubt[${collsyst}]} ${filerootDzero[${collsyst}]} 3 ${templateDirSystemSuffix[${collsyst}]}
	    elif [ ${imeson} = 1 ]; then
		$HFCJlocalCodeDir/DoSubtractFD.sh ${collsyst} ${imeson} ${dirDstarNotFDsubt[${collsyst}]}/${fpromptfileDstar[${collsyst}]} ${templateDir[${collsyst}]} ${dirDstarNotFDsubt[${collsyst}]} ${filerootDstar[${collsyst}]} 3 ${templateDirSystemSuffix[${collsyst}]}
	    elif [ ${imeson} = 2 ]; then
		$HFCJlocalCodeDir/DoSubtractFD.sh ${collsyst} ${imeson} ${dirDplusNotFDsubt[${collsyst}]}/${fpromptfileDplus[${collsyst}]} ${templateDir[${collsyst}]} ${dirDplusNotFDsubt[${collsyst}]} ${filerootDplus[${collsyst}]} 3 ${templateDirSystemSuffix[${collsyst}]}
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
			    $HFCJlocalCodeDir/DrawPlotInSubtractFD.sh ${dirppInput[${imeson}]}/${collsystdir[${collsyst}]}_${baseFile}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assoc${ptassoc[${iassocbin}]}${endFilepp} ${imeson} ${collsyst} ${itrigbin} ${iassocbin} ${reflect} ${rebin} 1
			    
			elif [ ${imeson} = 1 ]; then
			#echo "ciao"
			    $HFCJlocalCodeDir/DrawPlotInSubtractFD.sh ${dirppInput[${imeson}]}/${collsystdir[${collsyst}]}_${baseFile}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assoc${ptassoc[${iassocbin}]}${endFilepp} ${imeson} ${collsyst} ${itrigbin} ${iassocbin} ${reflect} ${rebin} 1
			elif [ ${imeson} = 2 ]; then
			#echo "ciao"
			    $HFCJlocalCodeDir/DrawPlotInSubtractFD.sh ${dirppInput[${imeson}]}/${collsystdir[${collsyst}]}_${baseFile}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assoc${ptassoc[${iassocbin}]}${endFilepp} ${imeson} ${collsyst} ${itrigbin} ${iassocbin} ${reflect} ${rebin} 1
			fi
		else 
			if [ ${imeson} = 0 ]; then
			#echo "ciao"
			    $HFCJlocalCodeDir/DrawPlotInSubtractFD.sh ${dirpPbInput[${imeson}]}/${collsystdir[${collsyst}]}_${baseFile}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assoc${ptassoc[${iassocbin}]}.root ${imeson} ${collsyst} ${itrigbin} ${iassocbin} ${reflect} ${rebin} 1
			elif [ ${imeson} = 1 ]; then
			#echo "ciao"
			    $HFCJlocalCodeDir/DrawPlotInSubtractFD.sh ${dirpPbInput[${imeson}]}/${collsystdir[${collsyst}]}_${baseFile}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assoc${ptassoc[${iassocbin}]}.root ${imeson} ${collsyst} ${itrigbin} ${iassocbin} ${reflect} ${rebin} 1
			elif [ ${imeson} = 2 ]; then
			#echo "ciao"
			    $HFCJlocalCodeDir/DrawPlotInSubtractFD.sh ${dirpPbInput[${imeson}]}/${collsystdir[${collsyst}]}_${baseFile}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assoc${ptassoc[${iassocbin}]}.root ${imeson} ${collsyst} ${itrigbin} ${iassocbin} ${reflect} ${rebin} 1
			fi
		    fi
		    if [ -e CanvaAndVariedHisto${collsystdir[$collsyst]}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assocPt${ptassoc[${iassocbin}]}.root ]; then
			ln -s $PWD/CanvaAndVariedHisto${collsystdir[$collsyst]}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assocPt${ptassoc[${iassocbin}]}.root ${baseDir}/AllPlots/
			ln -s $PWD/CanvaAndVariedHisto${collsystdir[$collsyst]}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assocPt${ptassoc[${iassocbin}]}.eps ${baseDir}/AllPlots/
		    ln -s $PWD/CanvaAndVariedHisto${collsystdir[$collsyst]}${meson[$imeson]}Pt${pttrig[${itrigbin}]}assocPt${ptassoc[${iassocbin}]}.png ${baseDir}/AllPlots/
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
		$HFCJlocalCodeDir/DoAverages.sh ${collsyst} ${itrigbin} ${iassocbin} ${averageOpt} ${reflect} 5 ${baseDir}/AllPlots/ ${baseDir}/AllPlots/ ${baseDir}/AllPlots/ 0
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
if [ ${dofit} = 1 ]; then
    cd ${baseDir}/AllPlots/Averages/FitResults    
    while [ ${collsyst} -le ${lastcollsyst} ]; do
	$HFCJlocalCodeDir/DoFit.sh ${collsyst} ${reflect} ${averageOpt} ${baseDir}/AllPlots/Averages/  ${baseDir}/AllPlots/Averages/FitResults
	collsyst=${collsyst}+1
    done
fi
collsyst=${firstcollsyst}

######## NOW COMPARE DATA AND MC ################# 
if [ $doCompareWithMC = 1 ]; then
    cd ${baseDir}/AllPlots/Averages
    mkdir ComparisonToModels
    cd ComparisonToModels
    ### START FROM PP ###
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
SetSkip3to5pPb(kTRUE)
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
SetSkip3to5pPb(kTRUE)
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
SetSkip3to5pPb(kTRUE)
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

if [ $doComparepppPb = 1 ]; then
    cd ${baseDir}/AllPlots/Averages
    mkdir ComparisonPPtoPPB
    cd ComparisonPPtoPPB
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

    

fi
