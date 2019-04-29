#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliAnalysisTaskSEDmesonsFilterCJ+;
#pragma link C++ class AliAnalysisTaskFlavourJetCorrelations+;

#pragma link C++ class AliPicoJet+;
#pragma link C++ class AliPicoHeaderJet+;
#pragma link C++ class AliPicoHeaderV0+;
#pragma link C++ class AliPicoBase+;
#pragma link C++ class AliPicoV0+;
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
#pragma link C++ class AliAnalysisTaskDmesonJets::AliD0ExtendedInfoSummary+;
#pragma link C++ class AliAnalysisTaskDmesonJets::AliDStarInfoSummary+;
#pragma link C++ class AliAnalysisTaskDmesonJets::AliJetInfoSummary+;
#pragma link C++ class AliAnalysisTaskDmesonJets::AliJetInfoPbPbSummary+;
#pragma link C++ class AliAnalysisTaskDmesonJets::AliEventInfoSummary+;
#pragma link C++ class std::vector<AliAnalysisTaskDmesonJets::AliDmesonInfoSummary>+;
#pragma link C++ class std::vector<AliAnalysisTaskDmesonJets::AliDmesonMCInfoSummary>+;
#pragma link C++ class std::vector<AliAnalysisTaskDmesonJets::AliD0InfoSummary>+;
#pragma link C++ class std::vector<AliAnalysisTaskDmesonJets::AliD0ExtendedInfoSummary>+;
#pragma link C++ class std::vector<AliAnalysisTaskDmesonJets::AliDStarInfoSummary>+;
#pragma link C++ class std::vector<AliAnalysisTaskDmesonJets::AliJetInfoSummary>+;
#pragma link C++ class std::vector<AliAnalysisTaskDmesonJets::AliJetInfoPbPbSummary>+;
#pragma link C++ class AliAnalysisTaskDmesonJets::AnalysisEngine+;
#pragma link C++ class AliAnalysisTaskDmesonJets::AliHFJetDefinition+;
#pragma link C++ class std::pair<int,AliAnalysisTaskDmesonJets::AliDmesonJetInfo>+;
#pragma link C++ class std::pair<std::string,AliAnalysisTaskDmesonJets::AliJetInfo>+;
#pragma link C++ class AliAnalysisTaskDmesonJets+;
#pragma link C++ class AliAnalysisTaskHFJetIPQA+;

#pragma link C++ class AliAnalysisTaskDmesonJetsSub::AliJetInfo+;
#pragma link C++ class AliAnalysisTaskDmesonJetsSub::AliDmesonJetInfo+;
#pragma link C++ class AliAnalysisTaskDmesonJetsSub::AliDmesonInfoSummary+;
#pragma link C++ class AliAnalysisTaskDmesonJetsSub::AliDmesonMCInfoSummary+;
#pragma link C++ class AliAnalysisTaskDmesonJetsSub::AliD0InfoSummary+;
#pragma link C++ class AliAnalysisTaskDmesonJetsSub::AliD0ExtendedInfoSummary+;
#pragma link C++ class AliAnalysisTaskDmesonJetsSub::AliDStarInfoSummary+;
#pragma link C++ class AliAnalysisTaskDmesonJetsSub::AliJetInfoSummary+;
#pragma link C++ class AliAnalysisTaskDmesonJetsSub::AliJetInfoPbPbSummary+;
#pragma link C++ class AliAnalysisTaskDmesonJetsSub::AliEventInfoSummary+;
#pragma link C++ class std::vector<AliAnalysisTaskDmesonJetsSub::AliDmesonInfoSummary>+;
#pragma link C++ class std::vector<AliAnalysisTaskDmesonJetsSub::AliDmesonMCInfoSummary>+;
#pragma link C++ class std::vector<AliAnalysisTaskDmesonJetsSub::AliD0InfoSummary>+;
#pragma link C++ class std::vector<AliAnalysisTaskDmesonJetsSub::AliD0ExtendedInfoSummary>+;
#pragma link C++ class std::vector<AliAnalysisTaskDmesonJetsSub::AliDStarInfoSummary>+;
#pragma link C++ class std::vector<AliAnalysisTaskDmesonJetsSub::AliJetInfoSummary>+;
#pragma link C++ class std::vector<AliAnalysisTaskDmesonJetsSub::AliJetInfoPbPbSummary>+;
#pragma link C++ class AliAnalysisTaskDmesonJetsSub::AnalysisEngine+;
#pragma link C++ class AliAnalysisTaskDmesonJetsSub::AliHFJetDefinition+;
#pragma link C++ class std::pair<int,AliAnalysisTaskDmesonJetsSub::AliDmesonJetInfo>+;
#pragma link C++ class std::pair<std::string,AliAnalysisTaskDmesonJetsSub::AliJetInfo>+;
#pragma link C++ class AliAnalysisTaskDmesonJetsSub+;

#pragma link C++ class AliAnalysisTaskDmesonJetsDetectorResponse::ResponseEngine+;
#pragma link C++ class AliAnalysisTaskDmesonJetsDetectorResponse::AliDmesonMatchInfoSummary+;
#pragma link C++ class AliAnalysisTaskDmesonJetsDetectorResponse::AliD0MatchInfoSummary+;
#pragma link C++ class AliAnalysisTaskDmesonJetsDetectorResponse::AliDStarMatchInfoSummary+;
#pragma link C++ class std::pair<AliAnalysisTaskDmesonJets::ECandidateType_t,AliAnalysisTaskDmesonJetsDetectorResponse::ResponseEngine>+;
#pragma link C++ class std::pair<AliAnalysisTaskDmesonJetsSub::ECandidateType_t,AliAnalysisTaskDmesonJetsDetectorResponse::ResponseEngine>+;
#pragma link C++ class AliAnalysisTaskDmesonJetsDetectorResponse+;

#pragma link C++ class AliAnalysisTaskHFSubstructure+;


#pragma link C++ class AliHFAODMCParticleContainer+;
#pragma link C++ class AliHFTrackContainer+;
#pragma link C++ class AliDJetTTreeReader+;
#pragma link C++ class AliDJetTHnReader+;
#endif

#endif
