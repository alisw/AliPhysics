#if !defined(__CINT__) && !defined(__CLING__)
# error Not for compilation
#else 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliHLTSampleComponent1+;
#pragma link C++ class AliHLTSampleComponent2+;
#pragma link C++ class AliHLTSampleCalibrationComponent+;
#pragma link C++ class AliHLTSampleESDAnalysisComponent+;
#pragma link C++ class AliHLTSampleRawAnalysisComponent+;
#pragma link C++ class AliHLTSampleMonitoringComponent+;
#pragma link C++ class AliHLTAgentSample+;
#pragma link C++ class AliHLTSamplePreprocessor+;
#pragma link C++ class AliHLTDummyComponent+;

#endif // __CINT__
