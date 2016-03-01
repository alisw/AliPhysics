#if !defined(__CINT__) && !defined(__CLING__)
# error Not for compilation
#else 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliHLTRootObjectMergerComponent+;
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
#pragma link C++ class AliHLTGlobalCompareFlatComponent+;
#pragma link C++ class AliHLTAsyncTestComponent+;
#pragma link C++ class AliHLTAsyncCalibrationComponent+;
#pragma link C++ class AliHLTZeroComponent+;
#pragma link C++ class AliHLTGlobalFlatEsdConverterComponent+;
#pragma link C++ class AliFlatESDFriendTrack+;
#pragma link C++ class AliFlatTPCseed+;
#pragma link C++ class AliFlatESDEvent+;
#pragma link C++ class AliHLTGlobalEsdToFlatConverterComponent+;
#pragma link C++ class AliFlatESDFriend+;
#pragma link C++ class AliFlatESDTrack+;
#pragma link C++ class AliHLTGlobalFlatEsdTestComponent+;
#pragma link C++ class AliHLTAnalysisManager+;
#pragma link C++ class AliHLTAnalysisManagerComponent+;
#pragma link C++ class AliHLTLumiRegComponent+;
#pragma link C++ class AliHLTGlobalPromptRecoQAComponent+;
#pragma link C++ class AliAnalysisTaskExampleV+;
#ifdef ZMQ
#pragma link C++ class AliOptionParser+;
#pragma link C++ class AliHLTZMQsink+;
#pragma link C++ class AliHLTZMQsource+;
#endif
#endif // __CINT__
