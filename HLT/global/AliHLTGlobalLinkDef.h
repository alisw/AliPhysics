#if !defined(__CINT__) && !defined(__CLING__)
# error Not for compilation
#else 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliHLTGlobalEsdConverterComponent+;
#pragma link C++ class AliHLTGlobalTrackMergerComponent+;
#pragma link C++ class AliHLTGlobalTrackMerger+;
#pragma link C++ class AliHLTGlobalAgent+;
#pragma link C++ class AliHLTGlobalPreprocessor+;
#pragma link C++ class AliHLTGlobalDCSPublisherComponent+;
#pragma link C++ class AliHLTGlobalVertexerComponent+;
#pragma link C++ class AliHLTGlobalOfflineVertexerComponent+;
#pragma link C++ class AliHLTGlobalTrackMatcher+;
#pragma link C++ class AliHLTGlobalTrackMatcherComponent+;
#pragma link C++ class AliHLTGlobalVertexerHistoComponent+;
#pragma link C++ class AliHLTGlobalHistoComponent+;
#pragma link C++ class AliHLTGlobalHistoCollector+;
#pragma link C++ class AliHLTVertexFinderBase+;
#pragma link C++ class AliHLTPrimaryVertexFinderComponent+;
#pragma link C++ class AliHLTV0FinderComponent+;
#pragma link C++ class AliHLTV0HistoComponent+;
#pragma link C++ class AliHLTCaloHistoComponent+;
#pragma link C++ class AliHLTCaloHistoProducer+;
#pragma link C++ class AliHLTCaloHistoInvMass+;
#pragma link C++ class AliHLTCaloHistoMatchedTracks+;
#pragma link C++ class AliHLTCaloHistoClusterEnergy+;
#pragma link C++ class AliHLTCaloHistoCellEnergy+;
#pragma link C++ class AliHLTMultiplicityCorrelations+;
#pragma link C++ class AliHLTMultiplicityCorrelationsComponent+;

#endif // __CINT__
