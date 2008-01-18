// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 * 
 **************************************************************************/
#pragma link off all functions;
#pragma link off all globals;
#pragma link off all classes;


//================================
// base/
//================================

// AliEveEventManager
#pragma link C++ class  AliEveEventManager+;
#pragma link C++ global gEvent;

#pragma link C++ class AliEveKineTools+;

#pragma link C++ class AliEveVSDCreator+;

// Fit
#pragma link C++ class AliEveTrackFitter+;
#pragma link C++ class AliEveTrackFitterEditor+;

// AliEveJetPlane
#pragma link C++ class AliEveJetPlane+;
#pragma link C++ class AliEveJetPlaneGL+;
#pragma link C++ class AliEveJetPlaneEditor+;
#pragma link C++ class AliEveJetPlaneEditor::StaticDataWindow+;

// Removed. Messy code, tons of violations and incompatible with TEve
// classes. Author Ludovic Gaudichet left ALICE.
// Should be thoroughly revised.
//
// AliEveCascade
//#pragma link C++ class AliEveCascade+;
//#pragma link C++ class CascadeList+;
//#pragma link C++ class CascadeListEditor+;
//
// AliEveV0
//#pragma link C++ class AliEveV0+;
//#pragma link C++ class V0List+;
//#pragma link C++ class V0ListEditor+;


//================================
// detectors/
//================================

// ITS
#pragma link C++ class AliEveITSModuleSelection+;
#pragma link C++ class AliEveITSDigitsInfo+;
#pragma link C++ class AliEveITSModule+;
#pragma link C++ class AliEveDigitScaleInfo+;
#pragma link C++ class AliEveITSScaledModule+;
#pragma link C++ class AliEveITSScaledModuleEditor+;
#pragma link C++ class AliEveITSModuleStepper+;
#pragma link C++ class AliEveITSModuleStepperEditor;

// MUON
#pragma link C++ class AliEveMUONData+;
#pragma link C++ class AliEveMUONChamber+;
#pragma link C++ class AliEveMUONChamberData+;
#pragma link C++ class AliEveMUONChamberEditor+;
#pragma link C++ class AliEveMUONChamberGL+;
#pragma link C++ class AliEveMUONTrack+;

// PMD
#pragma link C++ class AliEvePMDModule+;
#pragma link C++ class AliEvePMDModuleEditor+;

// T0
#pragma link C++ class AliEveT0Module+;

// TPC
#pragma link C++ class AliEveTPCData+;

#pragma link C++ class AliEveTPCSectorData+;
#pragma link C++ class AliEveTPCSectorData::PadData;
#pragma link C++ class AliEveTPCSectorData::PadIterator;
#pragma link C++ class AliEveTPCSectorData::RowIterator;
#pragma link C++ class AliEveTPCSectorData::SegmentInfo;

#pragma link C++ class AliEveTPCSectorData::PadRowHack;

#pragma link C++ class AliEveTPCSectorViz+;
#pragma link C++ class AliEveTPCSectorVizEditor+;
#pragma link C++ class AliEveTPCSector2D+;
#pragma link C++ class AliEveTPCSector2DEditor+;
#pragma link C++ class AliEveTPCSector2DGL+;
#pragma link C++ class AliEveTPCSector3D+;
#pragma link C++ class AliEveTPCSector3DEditor+;
#pragma link C++ class AliEveTPCSector3DGL+;

#pragma link C++ class AliEveTPCLoader+;
#pragma link C++ class AliEveTPCLoaderEditor+;

// TRD
#pragma link C++ class AliEveTRDLoaderManager+;
#pragma link C++ class AliEveTRDLoaderManagerEditor+;
#pragma link C++ class AliEveTRDLoader+;
#pragma link C++ class AliEveTRDLoaderEditor+;
#pragma link C++ class AliEveTRDLoaderSim+;
#pragma link C++ class AliEveTRDLoaderSimEditor+;
#pragma link C++ class AliEveTRDLoaderRaw+;
//#pragma link C++ class TRDLoaderRawEditor+;
#pragma link C++ class AliEveTRDModule+;
#pragma link C++ class AliEveTRDChamber+;
#pragma link C++ class AliEveTRDNode+;
#pragma link C++ class AliEveTRDModuleEditor+;
#pragma link C++ class AliEveTRDDigits+;
#pragma link C++ class AliEveTRDDigitsEditor+;
#pragma link C++ class AliEveTRDHits+;
#pragma link C++ class AliEveTRDHitsEditor+;
#pragma link C++ class AliEveTRDClusters+;

// TOF
#pragma link C++ class AliEveTOFDigitsInfo+;
#pragma link C++ class AliEveTOFSector+;
#pragma link C++ class AliEveTOFStrip+;

#pragma link C++ class AliEveTOFDigitsInfoEditor+;
#pragma link C++ class AliEveTOFSectorEditor+;
#pragma link C++ class AliEveTOFStripEditor+;


//================================
// HLT/
//================================

#pragma link C++ class AliEveHOMERManager+;
#pragma link C++ class AliEveHOMERManagerEditor+;
#pragma link C++ class AliEveHOMERSource+;
#pragma link C++ class AliEveHOMERSourceList+;
