/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///
/// \file doxymodules.h
/// \brief The definitions of CaloTrackCorr analysis class categories for Doxygen
///

/**
 * \defgroup CaloTrackCorrelationsAnalysis CaloTrackCorr Analysis Tasks
 * \ingroup CaloTrackCorrelations
 * \brief Analysis package for calorimeter PID and correlation with tracks, main analysis classes
 */

/**
 * \defgroup CaloTrackCorrMacros CaloTrackCorr Analysis Configuration Macros
 * \ingroup CaloTrackCorrelations
 * \brief Analysis task configuration macros for CaloTrackCorr package
 */

/**
 * \defgroup CaloTrackCorrMacrosQA CaloTrackCorr QA Configuration and Postprocessing Macros
 * \ingroup CaloTrackCorrelations
 * \brief Analysis task configuration and postprocessing macros for CaloTrackCorr package devoted to QA
 */

/**
 * \defgroup CaloTrackCorrMacrosQAPtHard Postprocessing files for Pt hard binned productions
 * \ingroup CaloTrackCorrMacrosQA
 * \brief Postprocessing files for EMCal analysis QA wagon for Pt hard binned productions
 *
 * The macros and scritps recover the merged output per pT hard bin and run, merges per run, 
 * extracts the wagon histograms and scales each pT hard bin by its cross section and produce the QA plots
 *
 * The files are:
 *  + DownloadExtractScaleMergePtHardAnalysisFiles.sh : execute the chain of download, mergin and scaling
 *  + mergePtHardRunFiles.sh : merge outpout histogram files in same pT hard bin
 *  + ScaleExtractPtHardBinHistograms.C : extract EMCal QA histograms to separate in 2 files, one scaled with cross section
 *  + DrawPtHardBins.C : read the scaled and unscaled merged files and do QA plots
 *
 */

