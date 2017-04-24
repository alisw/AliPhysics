#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliAnalysisTaskSEDmesonsFilterCJ+;
#pragma link C++ class AliAnalysisTaskFlavourJetCorrelations+;

#pragma link C++ class AliPicoJet+;
#pragma link C++ class AliPicoJetHeader+;
#pragma link C++ class AliPicoHeaderCJ+;
#pragma link C++ class AliPicoV0Base+;
#pragma link C++ class AliPicoV0MC+;
#pragma link C++ class AliPicoV0RD+;
#pragma link C++ class AliAnalysisTaskEmcalJetV0CF+;
#pragma link C++ class AliAnalysisTaskEmcalJetV0Filter+;
#pragma link C++ class AliAnalysisTaskEmcalJetHF+;
#pragma link C++ class AliAnalysisTaskSEPicoV0Filter+;
#pragma link C++ class AliAnalysisTaskSEPicoV0Maker+;
#pragma link C++ class AliAnalysisTaskSEPicoV0MakerMC+;
#pragma link C++ class AliAnalysisTaskEmcalJetFlavourTagExample+;
#pragma link C++ class AliMCHFParticleSelector+;

#ifdef HAVE_FASTJET
#pragma link C++ class AliAnalysisTaskDmesonJets::AliJetInfo+;
#pragma link C++ class AliAnalysisTaskDmesonJets::AliDmesonJetInfo+;
#pragma link C++ class AliAnalysisTaskDmesonJets::AliDmesonInfoSummary+;
#pragma link C++ class AliAnalysisTaskDmesonJets::AliDmesonMCInfoSummary+;
#pragma link C++ class AliAnalysisTaskDmesonJets::AliD0InfoSummary+;
#pragma link C++ class AliAnalysisTaskDmesonJets::AliDStarInfoSummary+;
#pragma link C++ class AliAnalysisTaskDmesonJets::AliJetInfoSummary+;
#pragma link C++ class AliAnalysisTaskDmesonJets::AliJetInfoPbPbSummary+;
#pragma link C++ class AliAnalysisTaskDmesonJets::AnalysisEngine+;
#pragma link C++ class AliAnalysisTaskDmesonJets::AliHFJetDefinition+;
#pragma link C++ class std::pair<int,AliAnalysisTaskDmesonJets::AliDmesonJetInfo>+;
#pragma link C++ class std::pair<std::string,AliAnalysisTaskDmesonJets::AliJetInfo>+;
#pragma link C++ class AliAnalysisTaskDmesonJets+;
#pragma link C++ class AliAnalysisTaskDmesonJetsDetectorResponse::ResponseEngine+;
#pragma link C++ class AliAnalysisTaskDmesonJetsDetectorResponse::AliDmesonMatchInfoSummary+;
#pragma link C++ class AliAnalysisTaskDmesonJetsDetectorResponse::AliD0MatchInfoSummary+;
#pragma link C++ class AliAnalysisTaskDmesonJetsDetectorResponse::AliDStarMatchInfoSummary+;
#pragma link C++ class std::pair<AliAnalysisTaskDmesonJets::ECandidateType_t,AliAnalysisTaskDmesonJetsDetectorResponse::ResponseEngine>+;
#pragma link C++ class AliAnalysisTaskDmesonJetsDetectorResponse+;
#pragma link C++ class AliHFAODMCParticleContainer+;
#pragma link C++ class AliHFTrackContainer+;
#pragma link C++ class AliDJetTTreeReader+;
#pragma link C++ class AliDJetTHnReader+;
#endif

#endif
