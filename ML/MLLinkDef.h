#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliExternalBDT+;

/// classes working in  ROOT6 only
#ifdef __CLING__
#pragma link C++ class AliMLResponse+;
#pragma link C++ class AliMLModelHandler+;
#endif

#endif
