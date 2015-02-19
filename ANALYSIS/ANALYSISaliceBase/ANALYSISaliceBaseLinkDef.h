#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliMultiInputEventHandler+;

// Only if ROOT was configured with AliEn support
#ifdef WITHALIEN
#pragma link C++ class AliAnalysisAlien+;
#endif  // WITHALIEN

#endif  // __CINT__
