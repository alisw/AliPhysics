#if !defined(__CINT__) && !defined(__CLING__)
# error Not for compilation
#else 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliHLTEveBase+;
#pragma link C++ class AliHLTEveCalo+;
#pragma link C++ class AliHLTEvePhos+;
#pragma link C++ class AliHLTEveEmcal+;
#pragma link C++ class AliHLTEveTPC+;
#pragma link C++ class AliHLTEveHLT+;
#pragma link C++ class AliHLTEveITS+;
#pragma link C++ class AliHLTEveISSD+;
#pragma link C++ class AliHLTEveISPD+;
#pragma link C++ class AliHLTEveISDD+;
#pragma link C++ class AliHLTEveTRD+;
#pragma link C++ class AliHLTEveMuon+;
#pragma link C++ class AliHLTEveAny+;
#pragma link C++ class AliHLTEveMultCorr+;
#pragma link C++ class AliHLTEveHistoMerger+;

#endif // __CINT__
