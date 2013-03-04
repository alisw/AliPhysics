#ifdef __CINT__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// Data Structure for Simulated Data 
#pragma link C++ class AliSegmentID+;                  // Base class for AliDigts and AliClusters (before also used for TRD)
#pragma link C++ class AliSegmentArray-;               // Keeps all AliDigits (all rows) (before also used for TRD)
#pragma link C++ class AliDigits+;                     // Digits per row -> (pad,timebin) (before also used for TRD)

#pragma link C++ class AliH2F+;                        // Additional functionality to 2D Histogram (used in Draw padResponse func)
                                                       // --- remove it, check miminal code needed for drawing

#pragma link C++ class AliTPCLoader+;                  // TPC Loader - derived from AliLoader ... 

#pragma link C++ class AliTPCRF1D-;                    // 1D Response Function (used for Time Response Function)
#pragma link C++ class AliTPCPRF2D-;                   // 2D Pad Response Function

#pragma link C++ class AliDetectorParam+;              // Base class for AliTPCParam (before also used for TRD)
#pragma link C++ class AliTPCParam+;                   // Parameterize the Geometry, Diffusion, ResponseFunction, Default HV, ...
                                                       // Base class for AliTPCParamSR
#pragma link C++ class AliTPCParamSR-;                 // SR = Straight Rows
                                                       // --- In principle only 1 class of (AliDetectorParam, AliTPCParam, 
                                                       //     AliTPCParamSR) is needed - can be merged, but breaks OCDB
 
#pragma link C++ class AliSimDigits+;                  // Derived from AliDigits - MC track labels in addition
                                                       //  --- Maybe combine AliDigits and AliSimDigits to new AliTPCdigits
#pragma link C++ class AliDigitsArray+;                // Derived from AliSegmentArray - Adds only AliDetecorParam
                                                       // -> Keeps AliDigits (all rows)
                                                       // --- Is this ptr still nedded? use singleton 
#pragma link C++ class AliTPCDigitsArray+;             // --- merge with AliDigitsArray -> new name AliTPCDigitsArray

#pragma link C++ class AliTPCROC+;                     // Geometry for 1 ROC (ReadOutChamber) - hardcoded
                                                       // --- (possible) duplication of AliTPCParam 
#pragma link C++ class AliTPCmapper+;                  // Hardware address mapping 
                                                       // --- investigate if it can be merged with AliTPCROC and AliTPCParam
#pragma link C++ class AliTPCCalROC-;                  // Calibration for 1 ROC - contains 1 entry per pad in a ROC
#pragma link C++ class AliTPCCalPad+;                  // Calibration for 1 TPC - contains AliTPCCalROC -> all pads
#pragma link C++ class AliTPCcalibDB+;                 // Main class to access OCDB and derived info - caches those
#pragma link C++ class AliTPCcalibDButil+;             // Helper methods for AliTPCcalibDB
#pragma link C++ class AliTPCSensorTemp+;              // Temperature info from DCS of 1 sensor- Derived from AliDCSSensor
#pragma link C++ class AliTPCSensorTempArray+;         // Array of AliTPCSensorTemp sensors

#pragma link C++ class AliTPCAltroMapping+;            // Maps Altro and physical HW address - used by AliTPCmapper
#pragma link C++ class AliTPCRawStreamV3+;             // TPC interface to Altro RAW decoder - using RCU format 3

// Calibration classes using RAW data (used in DA and HLT)
#pragma link C++ class AliTPCCalibRawBase+;            // Base class for the next 4 classes
#pragma link C++ class AliTPCCalibPedestal+;           // Calib : Calculate pedestal and noise from BlackEvents
                                                       //         pedestal value are wrong by 0.5 ADC counts
                                                       // --- To be checked with Christian and made consistent
#pragma link C++ class AliTPCCalibPulser+;             // Calib : Analysis of Pulser events (timing, amplitude 
                                                       //         and dead channel detection)
#pragma link C++ class AliTPCCalibCE+;	               // Calib : Drift velocity from CentralElectrode
                                                       //         + new method using laser tracks - not used so far in reco
#pragma link C++ class AliTPCCalibRaw+;	               // Monitors/analyzes Altro header info (eg Altro phase)

#pragma link C++ class AliTPCPreprocessor+;            // The Preprocessor
#pragma link C++ class AliTPCPreprocessorOnline+;      // Combines calibration info on pad level from 
                                                       // different sources, creates tree
                                                       // --- rename with meaningfull name"AliTPCCalibTreeCreator"
                                                       //     check overlap with calibviewer,  move MakeTree here
#pragma link C++ class AliTPCCalibViewer+;                     // Uses the trees from above -- move to Util
#pragma link C++ class AliTPCCalibViewerGUI+;                  // GUI for AliTPCCalibViewer -- move to Util
#pragma link C++ class AliTPCCalibViewerGUItime+;              // Timing trend for OCDB entries 
                                                               // --- rename like "CalibGUITimeDependent"  -- move to Util
#pragma link C++ class AliTPCCalibViewerGUItimeAddAliasFrame+; // Helper class for above  -- move to util
#pragma link C++ class AliTPCConfigDA+;                // Config Parser for DA configuration (= file in DAQ DB)
#pragma link C++ class AliTPCConfigParser+;            // same as above
                                                       // --- Jens to check the differences - one might obsolete

#pragma link C++ class AliTPCCalibVdrift+;             // Describes v_D(E,B,T,p) used in CalibDB
#pragma link C++ class AliTPCTempMap+;                 // Calculate temperature map and the gradient of the temperature

#pragma link C++ class AliTPCExB+;                     // Base class of AliTPCExBExact, AliTPCExBFirst
#pragma link C++ class AliTPCExBFirst+ ;               // Used only in MC - in CalibDB 
#pragma link C++ class AliTPCExBExact+;                // Benchmark for EXBFirst, to confirm its output
                                                       // --- after moving to RAW OCDB -> move all 3 to attic

#pragma link C++ class AliTransform+;                  // Base class for AliTPCTransform 
                                                       // --- combine with AliTPCTransform
#pragma link C++ class AliTPCTransform+;               // Full transformation of the space points 

#pragma link C++ class AliTPCdataQA-;                  // Used in DA, DQM, QA
#pragma link C++ class AliTPCQAChecker+;               // Offline QA framework

#pragma link C++ class AliTPCPointCorrection+;         // remove calling code form Transfrom : do quadrant alignment
                                                       // --- move to attic
#pragma link C++ class AliTPCLaserTrack+;              // Stores mirror positions and track angles of the Laser

#pragma link C++ class AliXRDPROOFtoolkit+;            // Toolkit for file tree based checks

#pragma link C++ class AliTPCCorrection+;              // Base class for space point distortions
#pragma link C++ class AliTPCInverseCorrection+;       // Inverse of AliTPCCorrection
#pragma link C++ class AliTPCComposedCorrection+;      // List of AliTPCCorrection - is itself a AliTPCCorrection
#pragma link C++ class AliTPCCorrectionDrift+;         // Parameterization of v_D 

// different space point corrections - implementations of AliTPCCorrection
#pragma link C++ class AliTPCExBBShape+;               // ExB due to B field inhomogenity
#pragma link C++ class AliTPCExBTwist+;                // ExB due to misalignment of TPC in respect of B field axis
#pragma link C++ class AliTPCGGVoltError+;             // ExB due to misalignment/bad voltage matching of Gating Grid 
#pragma link C++ class AliTPCFCVoltError3D+;           // ExB due to misalignment of components of the Field Cage
                                                       //  resitor rods, FC strip clamps
#pragma link C++ class AliTPCROCVoltError3D+;          // ExB due to ROC misalignment
#pragma link C++ class AliTPCBoundaryVoltError+;       // ExB due to voltage errors on main boundaries of the TPC
#pragma link C++ class AliTPCCalibGlobalMisalignment+; // Alignment expressed in terms of AliTPCCorrection 
                                                       //  + OROC qudarant (OROC has 4 separate pad planes) alignment 
#pragma link C++ class AliTPCSpaceCharge+;             // Distortions due to space charge in the TPC - rotational symetric
#pragma link C++ class AliTPCSpaceCharge3D+;           // Distortions due to space charge in the TPC - 3D calculation

#pragma link C++ class AliTPCExBEffective+;            // Cover ExB effect of non-explained physical model - not used
                                                       // --- still used in CalibMacros --- move to attic if removed there
#pragma link C++ class AliTPCExBEffectiveSector+;      // sectorwise above
                                                       // --- still used in CalibMacros --- move to attic if removed there

#endif
