/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///
/// \file doxymodules.h
/// \brief The definitions of EMCAL performance/calibration analysis class categories for Doxygen
///

/**
 * \defgroup EMCALPerformance EMCALPerformance
 * \brief Analysis package for EMCal performance, calibrations tasks
 */

/**
 * \defgroup EMCALPerfAddTaskMacros EMCALPerfAddTaskMacros
 * \ingroup EMCALPerformance EMCALPerformance
 * \brief Analysis task configuration macros for EMCal performance package
 */

/**
 * \defgroup EMCALOfflineMacros EMCALOfflineMacros
 * \ingroup EMCALPerformance EMCALPerformance
 * \brief Macros used for data analysis of train outputs for EMCal calibration
 *
 * This is a collection of analysis macros and general code used for EMCal QA and Calibration.
 * They are grouped in the following folder structure.
 *
 * + QAMacros contains macros for a run-by-run quality assurance of the EMcal performance
 * + TiCMacros contains macros for time calibration of the EMCal cell signals
 * + TeCMacros contains macros for temperature of the EMCal cell signals
 * + ECPi0Macros contains macros for energy calibration of the EMCal
 * + ECSingleChMacros contains macros for single channel energy calibration of the EMCal
 * + BCMacros contains macros for masking bad cells in the EMCal
 * 
 *
 * QA code
 * ---------------------
 * + runEMCALQA.pl
 * + PlotEMCALQATrendingTree.C
 * + MakeQAPdf.C
 * + CreateEMCALRunQA.C
 * + CopyQAFile.C
 *
 * Time Calibration code
 * ---------------------
 * + macro1.C
 * + macro2.C
 * + macro3.C
 *
 *
 * Bad Channel Analysis code
 * ---------------------
 * + BadChannelAna.cxx
 * + BadChannelAna.h
 * + runCopyFromLegoTrainBC.sh 
 * + runRByRAnalysisBC.sh      
 * + helperMacrosBC.C          
 * + runAnalysisBC.C           
 * + helperMacrosOADBBC.C      
 * + helperMacrosRunByRunBC.C  
 * 
 *
 * Energy Calibration (pi0) code
 * ---------------------
 * + CheckNoisePeakVariationWithTimeCut.C
 * + ComparePi0CalibResults.C
 * + ConvertOCDBtoText.C
 * + HVrecalculation.C
 * + HVrecalculation_EMCALthirds.C
 * + MergeHVscanFilesVariousScans2.C
 * + MergeThirdSMfilesIntoOne.C
 * + MultiplyPi0CalibrationFactors_TextToHisto.C
 * + MultiplyPi0CalibrationFactors_TextToHisto_Final.C
 * + Pi0CalibInvMassAnalysis3.C
 * + PlotLEDruns.C
 * + PrepareHV_SM1819.C
 * + SetCDB.C
 * + UnderstandDCALthirdShuffle.C
 * + UpdateEMCAL_OADB_Recalib.C
 *
 *
 * Energy Calibration (Single Channel) code
 * ---------------------
 * + paint_emcal.py
 * + calib_emcal.py
 *
 *
 * Temperature Calibration code
 * ---------------------
 * + plot.C
 * + plotProf.C
 * + scanAll.C
 * + scanTemp.C
 *
 */
