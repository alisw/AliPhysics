#pragma link off all functions;
#pragma link off all globals;
#pragma link off all classes;

#pragma link C++ namespace Alieve;

//================================
// base/
//================================

// Event
#pragma link C++ class  Alieve::Event+;
#pragma link C++ global Alieve::gEvent;

#pragma link C++ class Alieve::KineTools+;

#pragma link C++ class Alieve::VSDCreator+;


//================================
// compound classes
//================================

// Cascade
#pragma link C++ class Alieve::Cascade+;
#pragma link C++ class Alieve::CascadeList+;
#pragma link C++ class Alieve::CascadeListEditor+;

// V0
#pragma link C++ class Alieve::V0+;
#pragma link C++ class Alieve::V0List+;
#pragma link C++ class Alieve::V0ListEditor+;


//================================
// detectors/
//================================

// ITS
#pragma link C++ class Alieve::ITSModuleSelection+;
#pragma link C++ class Alieve::ITSDigitsInfo+;
#pragma link C++ class Alieve::ITSModule+;
#pragma link C++ class Alieve::DigitScaleInfo+;
#pragma link C++ class Alieve::ITSScaledModule+;
#pragma link C++ class Alieve::ITSScaledModuleEditor+;
#pragma link C++ class Alieve::ITSModuleStepper+;
#pragma link C++ class Alieve::ITSModuleStepperEditor;

// MUON
#pragma link C++ class Alieve::MUONData+;
#pragma link C++ class Alieve::MUONChamber+;
#pragma link C++ class Alieve::MUONChamberData+;
#pragma link C++ class Alieve::MUONChamberEditor+;
#pragma link C++ class Alieve::MUONChamberGL+;
#pragma link C++ class Alieve::MUONTrack+;

// PMD
#pragma link C++ class Alieve::PMDModule+;
#pragma link C++ class Alieve::PMDModuleEditor+;

// T0
#pragma link C++ class Alieve::T0Module+;

// TPC
#pragma link C++ class Alieve::TPCData+;

#pragma link C++ class Alieve::TPCSectorData+;
#pragma link C++ class Alieve::TPCSectorData::PadData;
#pragma link C++ class Alieve::TPCSectorData::PadIterator;
#pragma link C++ class Alieve::TPCSectorData::RowIterator;
#pragma link C++ class Alieve::TPCSectorData::SegmentInfo;

#pragma link C++ class Alieve::TPCSectorData::PadRowHack;

#pragma link C++ class Alieve::TPCSectorViz+;
#pragma link C++ class Alieve::TPCSectorVizEditor+;
#pragma link C++ class Alieve::TPCSector2D+;
#pragma link C++ class Alieve::TPCSector2DEditor+;
#pragma link C++ class Alieve::TPCSector2DGL+;
#pragma link C++ class Alieve::TPCSector3D+;
#pragma link C++ class Alieve::TPCSector3DEditor+;
#pragma link C++ class Alieve::TPCSector3DGL+;

#pragma link C++ class Alieve::TPCLoader+;
#pragma link C++ class Alieve::TPCLoaderEditor+;

// TRD
#pragma link C++ class Alieve::TRDLoaderManager+;
#pragma link C++ class Alieve::TRDLoaderManagerEditor+;
#pragma link C++ class Alieve::TRDLoader+;
#pragma link C++ class Alieve::TRDLoaderEditor+;
#pragma link C++ class Alieve::TRDLoaderSim+;
#pragma link C++ class Alieve::TRDLoaderSimEditor+;
#pragma link C++ class Alieve::TRDLoaderRaw+;
//#pragma link C++ class Alieve::TRDLoaderRawEditor+;
#pragma link C++ class Alieve::TRDModule+;
#pragma link C++ class Alieve::TRDChamber+;
#pragma link C++ class Alieve::TRDNode+;
#pragma link C++ class Alieve::TRDModuleEditor+;
#pragma link C++ class Alieve::TRDDigits+;
#pragma link C++ class Alieve::TRDDigitsEditor+;
#pragma link C++ class Alieve::TRDHits+;
#pragma link C++ class Alieve::TRDHitsEditor+;
#pragma link C++ class Alieve::TRDClusters+;

// TOF
#pragma link C++ class Alieve::TOFDigitsInfo+;
#pragma link C++ class Alieve::TOFSector+;
#pragma link C++ class Alieve::TOFStrip+;

#pragma link C++ class Alieve::TOFDigitsInfoEditor+;
#pragma link C++ class Alieve::TOFSectorEditor+;
#pragma link C++ class Alieve::TOFStripEditor+;

// JetPlane
#pragma link C++ class Alieve::JetPlane+;
#pragma link C++ class Alieve::JetPlaneGL+;
#pragma link C++ class Alieve::JetPlaneEditor+;
#pragma link C++ class Alieve::JetPlaneEditor::StaticDataWindow+;
