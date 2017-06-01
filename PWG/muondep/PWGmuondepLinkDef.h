#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliAnalysisRunList+;

#pragma link C++ class AliAnalysisTaskESDMCLabelAddition+;
#pragma link C++ class AliAnalysisTaskPileup+;
#pragma link C++ class AliAnalysisTaskMuonCDBConnect+;
#pragma link C++ class AliAnalysisTaskMuonRefit+;
#pragma link C++ class AliAnalysisTaskMuonRefitVtx+;

#pragma link C++ class AliAnalysisMuMu+;
#pragma link C++ class AliAnalysisMuMuConfig+;
#pragma link C++ class AliAnalysisMuMuResult+;
#pragma link C++ class AliAnalysisMuMuJpsiResult+;
#pragma link C++ class AliAnalysisMuMuSpectraProcessor+;
#pragma link C++ class AliAnalysisMuMuSpectraProcessorPbPb+;
#pragma link C++ class AliAnalysisMuMuSpectraProcessorPbP+;
#pragma link C++ class AliAnalysisMuMuSpectraProcessorPP+;
#pragma link C++ class AliAnalysisMuMuFnorm+;
#pragma link C++ class AliAnalysisMuMuGraphUtil+;
#pragma link C++ class AliAnalysisMuMuSpectra+;
#pragma link C++ class AliAnalysisTriggerScalers+;
#pragma link C++ class AliAnalysisTriggerScalerItem+;

#pragma link C++ class AliMuonGridSubmitter+;
#pragma link C++ class AliMuonAccEffSubmitter+;
#pragma link C++ class AliMuonQAMergeSubmitter+;
#pragma link C++ class AliMuonOCDBSnapshotGenerator+;

#pragma link C++ class AliMuonCompactMapping+;
#pragma link C++ class AliMuonCompactEvent+;
#pragma link C++ class AliMuonCompactTrack+;
#pragma link C++ class AliMuonCompactCluster+;
#pragma link C++ class AliMuonCompactTreeMaker+;
#pragma link C++ class AliMuonCompactManuStatus+;
#pragma link C++ class AliMuonCompactQuickAccEff+;
#pragma link C++ class AliMuonCompactQuickAccEffChecker+;

#pragma link C++ function operator<<(std::ostream&, const AliMuonCompactMapping&);
#pragma link C++ function operator<<(std::ostream&, const AliMuonCompactTrack&);
#pragma link C++ function operator<<(std::ostream&, const AliMuonCompactEvent&);
#pragma link C++ function operator<<(std::ostream&, const AliMuonCompactCluster&);

#endif


