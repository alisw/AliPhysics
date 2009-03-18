#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class  AliAnalysisTaskSE+;
#pragma link C++ class  AliAnalysisTaskME+;
#pragma link C++ class  AliAnalysisTaskESDfilter+;
#pragma link C++ class  AliAnalysisTaskKineFilter+;
#pragma link C++ class  AliAnalysisTaskMCParticleFilter+;
#pragma link C++ class  AliAnalysisTaskTagCreator+;
#pragma link C++ class  AliKineTrackCuts+;
#pragma link C++ class  AliESDtrackCuts+;
#pragma link C++ class  AliESDv0Cuts+;

#pragma link C++ class  AliEventPoolOTF+;
#pragma link C++ class  AliEventPoolLoop+;

#pragma link C++ class AliEventPoolSparse+;

#pragma link C++ class AliMultiEventInputHandler+;

#ifdef WITHXML
#pragma link C++ class AliTagAnalysis+;
#pragma link C++ class AliXMLCollection+;
#pragma link C++ class AliAnalysisAlien+;
#endif

#endif
