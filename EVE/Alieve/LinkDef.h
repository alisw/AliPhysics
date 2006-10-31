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

#pragma link C++ class Alieve::VSDCreator+;

//================================
// detectors/
//================================

// ITS
#pragma link C++ class Alieve::ITSDigitsInfo+;
#pragma link C++ class Alieve::ITSModule+;

// MUON
#pragma link C++ class Alieve::MUONData+;
#pragma link C++ class Alieve::MUONChamber+;
#pragma link C++ class Alieve::MUONChamberData+;
#pragma link C++ class Alieve::MUONChamberEditor+;
#pragma link C++ class Alieve::MUONChamberGL+;

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
#pragma link C++ class Alieve::TRDLoader+;
#pragma link C++ class Alieve::TRDLoaderEditor+;
#pragma link C++ class Alieve::TRDModule+;
#pragma link C++ class Alieve::TRDChamber+;
#pragma link C++ class Alieve::TRDNode+;
#pragma link C++ class Alieve::TRDModuleEditor+;
#pragma link C++ class Alieve::TRDDigits+;
#pragma link C++ class Alieve::TRDHits+;
