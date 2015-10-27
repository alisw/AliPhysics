#if !defined(__CINT__) && !defined(__CLING__)
# error Not for compilation
#else 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliHLTEMCALRawAnalyzerComponent+;
#pragma link C++ class AliHLTEMCALMapper+;
#pragma link C++ class AliHLTEMCALRawAnalyzerCrudeComponent+;
#pragma link C++ class AliHLTEMCALRawAnalyzerLMSComponent+;
#pragma link C++ class AliHLTEMCALRawAnalyzerPeakFinderComponent+;
#pragma link C++ class AliHLTEMCALRawAnalyzerFastFitComponent+;
#pragma link C++ class AliHLTEMCALRawAnalyzerNNComponent+;
#pragma link C++ class AliHLTEMCALConstants+;
#pragma link C++ class AliHLTEMCALDigitMakerComponent+;
#pragma link C++ class AliHLTEMCALClusterizerComponent+;
#pragma link C++ class AliHLTEMCALRecoParamHandler+;
#pragma link C++ class AliHLTEMCALClusterizerComponentNbyN+;
#pragma link C++ class AliHLTEMCALClusterMonitorComponent+;
#pragma link C++ class AliHLTEMCALClusterMonitor+;
#pragma link C++ class AliHLTEMCALAgent+;
#pragma link C++ class AliHLTEMCALRawHistoMaker+;
#pragma link C++ class AliHLTEMCALRawHistoMakerComponent+;

#endif // __CINT__
