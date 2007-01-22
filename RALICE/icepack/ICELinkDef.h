///////////////////////////////////////////////////////////////////////////
// All classes of ICE analysis software
// This class list is used to create the ICE dictionary via rootcint.
//
// Note : Headers have also to be entered into the list
//        contained in ICEHeaders.h
//
//--- NvE 23-jan-2003 Utrecht University
///////////////////////////////////////////////////////////////////////////
 
#ifdef __CINT__
 #pragma link off all globals;
 #pragma link off all classes;
 #pragma link off all functions;
 
 #pragma link C++ class IceEvent+;
 #pragma link C++ class IceGOM+;
 #pragma link C++ class IceAOM+;
 #pragma link C++ class IceDOM+;
 #pragma link C++ class IceIDOM+;
 #pragma link C++ class IceTDOM+;
 #pragma link C++ class IceXtalk+;
 #pragma link C++ class IceCalibrate+;
 #pragma link C++ class IceCleanHits+;
 #pragma link C++ class IceDwalk+;
 #pragma link C++ class IcePandel+;
 #pragma link C++ class IceLinefit+;
 #pragma link C++ class IceChi2+;
 #pragma link C++ class IceDwalkx+;
 #pragma link C++ class IceMakeHits+;
#endif
 
