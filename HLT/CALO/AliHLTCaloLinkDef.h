#if !defined(__CINT__) && !defined(__CLING__)
# error Not for compilation
#else 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliHLTCaloRawAnalyzerComponentv3+;
#pragma link C++ class AliHLTCaloUtilities+;
#pragma link C++ class AliHLTCaloMapper+;
#pragma link C++ class AliHLTCaloDefinitions+;
#pragma link C++ class AliHLTCaloConstants+;
#pragma link C++ class AliHLTCaloSanityInspector+;
#pragma link C++ class AliHLTCaloSharedMemoryInterfacev2+;
#pragma link C++ class AliHLTCaloFourier+;
#pragma link C++ class AliHLTCaloConstantsHandler+;
#pragma link C++ class AliHLTCaloClusterizer+;
#pragma link C++ class AliHLTCaloClusterizerNbyN+;
#pragma link C++ class AliHLTClusterFinder+;
#pragma link C++ class AliHLTCaloClusterizerComponent+;
#pragma link C++ class AliHLTCaloDigitMaker+;
#pragma link C++ class AliHLTCaloClusterAnalyser+;
#pragma link C++ class AliHLTCaloProcessor+;
#pragma link C++ class AliHLTCaloGeometry+;
#pragma link C++ class AliHLTCaloRecoParamHandler+;
#pragma link C++ class AliHLTCaloDigitPublisherComponent+;
#pragma link C++ class AliHLTCaloDigitHandler+;
#pragma link C++ class AliHLTEMCALDefinitions+;
#pragma link C++ class AliHLTEMCALDigitHandler+;
#pragma link C++ class AliHLTEMCALGeometry+;
#pragma link C++ class AliHLTPHOSDefinitions+;
#pragma link C++ class AliHLTPHOSDigitHandler+;
#pragma link C++ class AliHLTPHOSGeometry+;

#endif // __CINT__
