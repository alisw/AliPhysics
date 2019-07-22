#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliDxHFEParticleSelection+;
#pragma link C++ class AliDxHFEParticleSelectionD0+;
#pragma link C++ class AliDxHFEParticleSelectionEl+;
#pragma link C++ class AliDxHFEParticleSelectionMCD0+;
#pragma link C++ class AliDxHFEParticleSelectionMCEl+;
#pragma link C++ class AliDxHFEToolsMC+;
#pragma link C++ class AliDxHFECorrelation+;
#pragma link C++ class AliDxHFECorrelationMC+;
#pragma link C++ class AliAnalysisTaskDxHFEParticleSelection+;
#pragma link C++ class AliAnalysisTaskDxHFECorrelation+;
#pragma link C++ class AliHFAssociatedTrackCuts+;
#pragma link C++ class AliHFCorrelator+;
#pragma link C++ class AliHFOfflineCorrelator+;
#pragma link C++ class AliHFCorrelationBranchD+;
#pragma link C++ class AliHFCorrelationBranchTr+;
#pragma link C++ class AliReducedParticle+;
#pragma link C++ class AliD0hCutOptim+;
#pragma link C++ class AliDstarhCutOptim+;
#pragma link C++ class AliDPlushCutOptim+;
#pragma link C++ class AliAnalysisTaskDStarCorrelations+;
#pragma link C++ class AliAnalysisTaskSED0Correlations+;
#pragma link C++ class AliAnalysisTaskSED0CorrelationsVsMult+;
#pragma link C++ class AliAnalysisTaskSEDplusCorrelations+;
#pragma link C++ class AliAnalysisTaskSEmcCorr+;	
#pragma link C++ class AliAnalysisTaskSEHFCJqa+;
#pragma link C++ class AliHFDhadronCorrSystUnc+;	
#pragma link C++ class AliHFCorrelationFDsubtraction+;	
#pragma link C++ class AliHFDmesonCorrAverage+;
#pragma link C++ class AliHFCorrelationUtils+;
#pragma link C++ class AliAnalysisHFCorrOnFlySim+;
#pragma link C++ class AliHFCorrFitter+;
#pragma link C++ class AliHFCorrFitSystematics+;
#pragma link C++ class AliDhCorrelationExtraction+;

#if ROOT_VERSION_CODE > ROOT_VERSION(6,4,0)
#pragma link C++ namespace AliDHFeCorr+;
#pragma link C++ typedef AliDHFeCorr::AliDMeson+;
#pragma link C++ typedef AliDHFeCorr::AliElectron+;
#pragma link C++ typedef AliDHFeCorr::AliEventSelection+;
#pragma link C++ typedef AliDHFeCorr::AliElectronSelection+;
#pragma link C++ typedef AliDHFeCorr::AliPhotonSelection+;
#pragma link C++ typedef AliDHFeCorr::AliDMesonSelection+;
#pragma link C++ typedef AliDHFeCorr::AliEventQAHistograms+;
#pragma link C++ typedef AliDHFeCorr::AliElectronQAHistograms+;
#pragma link C++ typedef AliDHFeCorr::AliDMesonQAHistos+;
#pragma link C++ typedef AliDHFeCorr::AliConfigureElectronOpt+;
#pragma link C++ typedef AliDHFeCorr::AliConfigureDMesonOpt+;
#pragma link C++ class AliAnalysisTaskDHFeCorr+;
#endif

#endif
