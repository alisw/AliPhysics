#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class std::vector < float> + ;
#pragma link C++ class std::pair < int, std::vector < float>> + ;
#pragma link C++ class std::pair < TString, int> + ;
#pragma link C++ class std::pair < TString, TString> + ;
#pragma link C++ class std::pair < float, float> + ;
#pragma link C++ unordered_map < int, std::pair < float, float>> + ;
#pragma link C++ enum ParticleType;
#pragma link C++ class std::vector < ParticleType> + ;
#pragma link C++ class MCParticleHists + ;
#pragma link C++ class std::pair < ParticleType, MCParticleHists> + ;
#pragma link C++ class std::map < ParticleType, MCParticleHists> + ;
#pragma link C++ class DataParticleHists + ;
#pragma link C++ class std::pair < ParticleType, DataParticleHists> + ;
#pragma link C++ class std::map < ParticleType, DataParticleHists> + ;

#pragma link C++ class AliAnalysisTaskCMEv2A + ;
#pragma link C++ class AliAnalysisTaskCMEPIDCVE + ;
#pragma link C++ class AliAnalysisTaskPIDCME + ;
#pragma link C++ class AliAnalysisTaskCVEPIDCMEDiff + ;
#pragma link C++ class AliAnalysisTaskCVEUtil + ;
#endif
