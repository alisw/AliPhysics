#ifdef __CINT__
 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// Base classes
#pragma link C++ class AliConversionPhotonBase++;
#pragma link C++ class AliAODConversionParticle+;
#pragma link C++ class AliAODConversionMother++;
#pragma link C++ class AliAODConversionPhoton++;
#pragma link C++ class AliKFConversionPhoton++;
#pragma link C++ class AliKFConversionMother++;
#pragma link C++ class AliConversionCuts++;
#pragma link C++ class AliConversionSelection++;
#pragma link C++ class AliV0ReaderV1+;
#pragma link C++ class AliConversionAODBGHandlerRP+;
#pragma link C++ class AliConversionTrackCuts++;

// User tasks
#pragma link C++ class AliAnalysisTaskGammaConvDalitz+;
#pragma link C++ class AliAnalysisTaskPi0v2+;

#pragma link C++ class AliAnaConvIsolation+;
#pragma link C++ class AliAnaConvCorrBase++;
#pragma link C++ class AliAnaConvCorrPion++;
#pragma link C++ class AliAnaConvCorrPhoton++;
#pragma link C++ class AliAnalysisTaskdPhi+;

// Old tasks
#pragma link C++ class AliAnalysisTaskGammaConversion+;
#pragma link C++ class AliAnalysisTaskCaloConv+;
#pragma link C++ class AliV0Reader+;
#pragma link C++ class AliGammaConversionHistograms+;
#pragma link C++ class AliGammaConversionBGHandler+;
#pragma link C++ class AliAnalysisTaskGCPartToPWG4Part+;



#endif
