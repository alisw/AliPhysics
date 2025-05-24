#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class std::vector < float> + ;
#pragma link C++ class std::pair < int, std::vector < float>> + ;
#pragma link C++ class std::pair < TString, int> + ;
#pragma link C++ class std::pair < TString, TString> + ;

#pragma link C++ class AliAnalysisTaskCMEv2A + ;
#pragma link C++ class AliAnalysisTaskCMEPIDCVE + ;
#pragma link C++ class AliAnalysisTaskPIDCME + ;
#pragma link C++ class AliAnalysisTaskCVEPIDCMEDiff + ;
#pragma link C++ class AliAnalysisTaskCVEUtil + ;
#endif
