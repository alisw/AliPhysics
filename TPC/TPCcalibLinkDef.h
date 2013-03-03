#ifdef __CINT__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
 
// add calib base here


#pragma link C++ class AliTPCAnalysisTaskcalib+;   // Task for processing TPC calibartion, contains arrays of calib base



#pragma link C++ class AliTPCcalibTracks+;         // Histogram cluster to track residuals in order to get
                                                    // clustershape and cluster error paramterization
                                                    // results go to AliTPCClusterParam
                                                    // --- clean up
#pragma link C++ class AliTPCcalibTracksCuts+;     // specify input (cuts)for AliTPCcalibTracks and AliTPCcalibTracksGain
#pragma link C++ class AliTPCcalibTracksGain+;     // Get gain calibration per pad
                                                    // --- clean up
#pragma link C++ class AliTPCFitPad+;              // Helper class (fitting) for the above 
                                                    // duplicted also (most probably) AliTPCCalPad --- to be cleaned up
#pragma link C++ class AliTPCCalPadRegion+;        // Helper class (fitting) for the above 
                                                   // duplicted also (most probably) AliTPCCalPad --- to be cleaned up

#pragma link C++ class AliTPCcalibCalib+;          // re-applying calib on cluster level - refitting of the tracks
#pragma link C++ class AliTPCcalibAlign+;          // (sector) alignment calibration 
                                                   //  Histogram cluster to track residuals in order to 
                                                   //  create later on distortion maps
                                                   // --- update documentation (current docu is obsolete)
#pragma link C++ class AliTPCcalibAlignment+;      // obsolete --- remove
#pragma link C++ class AliTPCcalibV0+;             // Used Filterin of V0s - to be removed later
#pragma link C++ class AliTPCcalibCosmic+;         // Used Filterin of Cosmics - to be removed later

#pragma link C++ class AliTPCCalibKr+;             // Krypton Calibration
#pragma link C++ class AliTPCCalibKrTask+;         // Analysis task for  Krypton Calibration 
#pragma link C++ class AliTPCcalibLaser+;          // Histograms residuals of clusters to ideal laser tracks

#pragma link C++ class AliTPCcalibTime+;           // Calibrates v_Drift + Histograms residuals of track extrapolations to outer detectors 
                                                   // output is used in OCDB object TPC/Calib/TimeDrift
#pragma link C++ class AliTPCcalibTimeGain+;       // gas gain calibartion and multiplicity correction, normalizing the MIP - per padregion and sector
                                                   // output is used in OCDB object TPC/Calib/TimeGain

#pragma link C++ class AliTPCMisAligner+;          // remove --- after checking with Raffaele Grosso
                                                   // documentation needed
#pragma link C++ class AliTPCcalibTrigger+;        // obsolete --- remove

#pragma link C++ class AliTPCPreprocessorOffline+; // proccess output of calibration and create OCDB entry
#pragma link C++ class AliTPCcalibGainMult+;       // not used for the moment / contains 'new' dE/dx algorithm

#pragma link C++ class AliTPCCorrectionFit;        // Fitting Methods for space-point calibration classes
#pragma link C++ class AliTPCkalmanAlign+;         // combines relative alignment with global alignmet
                                                   // --- move functionality to AliTPCCorrectionFit

#pragma link C++ class AliTPCcalibSummary;         // tree creation of calibration parameters

#endif





