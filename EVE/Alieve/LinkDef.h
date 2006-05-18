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

// ITS

#pragma link C++ class Alieve::ITSDigitsInfo+;
#pragma link C++ class Alieve::ITSModule+;

// TPC

#pragma link C++ class Alieve::TPCData+;

#pragma link C++ class Alieve::TPCSectorData+;
#pragma link C++ class Alieve::TPCSectorData::PadData;
#pragma link C++ class Alieve::TPCSectorData::PadIterator;
#pragma link C++ class Alieve::TPCSectorData::RowIterator;
#pragma link C++ class Alieve::TPCSectorData::SegmentInfo;

#pragma link C++ class Alieve::TPCSector2D+;
#pragma link C++ class Alieve::TPCSector2DEditor+;
#pragma link C++ class Alieve::TPCSector2DGL+;

// Almost ready
// #pragma link C++ class Alieve::TPCSector3D+;
// #pragma link C++ class Alieve::TPCSector3DEditor+;
