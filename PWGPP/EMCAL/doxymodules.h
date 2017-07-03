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
 * + QA contains macros for a run-by-run quality assurance of the EMcal performance
 * + TiC contains macros for time calibration of the EMCal cell signals
 * + TeC contains macros for temperature of the EMCal cell signals
 * + EC contains macros for energy calibration of the EMCal
 * + BC contains macros for masking bad cells in the EMCal
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
 * 
 * + macro1.C
 * + macro2.C
 * + macro3.C
 *
 * Bad Channel Analysis code
 * ---------------------
 *
 * + BadChannelAna.cxx
 * + BadChannelAna.h
 * + runCopyFromLegoTrainBC.sh 
 * + runRByRAnalysisBC.sh      
 * + helperMacrosBC.C          
 * + runAnalysisBC.C           
 * + helperMacrosOADBBC.C      
 * + helperMacrosRunByRunBC.C  
 * 
 * Energy Calibration code
 * ---------------------
 *
 * + macro1.C
 * + macro2.C
 * + macro3.C
 *
 * Temperature Calibration code
 * ---------------------
 *
 * + macro1.C
 * + macro2.C
 * + macro3.C
 *
 */
