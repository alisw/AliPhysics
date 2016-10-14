#ifdef __CINT__
 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliCFPI0+;
#pragma link C++ class AliFilteredTrack+;
#pragma link C++ class AliFilteredEvent+;
#pragma link C++ class AliFilteredEventInputHandler+;
#pragma link C++ class AliCorrelation3p_noQA+;
#pragma link C++ class AliCorrelation3p+;
#pragma link C++ class AliThreeParticleCorrelator<AliCorrelation3p_noQA>+;
#pragma link C++ class AliThreeParticleCorrelator<AliCorrelation3p>+;
#pragma link C++ class AliAnalysisTaskCorrelation3p+;
#pragma link C++ class AliAnalysisTaskCorrelation3p_lightefficiency+;
#pragma link C++ class AliAnalysisTaskBuildCorrTree+;
#endif
