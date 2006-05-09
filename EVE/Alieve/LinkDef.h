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

#pragma link C++ class Alieve::ITSDigitsInfo+;
#pragma link C++ class Alieve::TPCDigitsInfo+;
#pragma link C++ class Alieve::TPCSeg+;

//================================
// g3d/
//================================

#pragma link C++ class Alieve::ITSModule+;
#pragma link C++ class Alieve::TPCSegment+;

//================================
// ged/
//================================

#pragma link C++ class Alieve::TPCSegmentEditor+;

//================================
// gl/
//================================

#pragma link C++ class Alieve::TPCSegmentGL+;
