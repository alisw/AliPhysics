#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliFemtoDreamEvent+;
#pragma link C++ class AliFemtoDreamEventHist+;
#pragma link C++ class AliFemtoDreamEventCuts+;
#pragma link C++ class AliFemtoDreamBasePart+;
#pragma link C++ class AliFemtoDreamTrack+;
#pragma link C++ class AliFemtoDreamTrackHist+;
#pragma link C++ class AliFemtoDreamTrackMCHist+;
#pragma link C++ class AliFemtoDreamTrackCuts+;
#pragma link C++ class AliFemtoDreamv0+;
#pragma link C++ class AliFemtoDreamv0Hist+;
#pragma link C++ class AliFemtoDreamv0MCHist+;
#pragma link C++ class AliFemtoDreamv0Cuts+;
#pragma link C++ class AliFemtoDreamControlSample+;
#pragma link C++ class AliFemtoDreamCascade+;
#pragma link C++ class AliFemtoDreamCascadeHist+;
#pragma link C++ class AliFemtoDreamCascadeCuts+;
#pragma link C++ class AliFemtoDreamPairCleanerHists+;
#pragma link C++ class AliFemtoDreamPairCleaner+;
#pragma link C++ class AliFemtoDreamCollConfig+;
#pragma link C++ class AliFemtoDreamCorrHists+;
#pragma link C++ class AliFemtoDreamPartContainer+;
#pragma link C++ class AliFemtoDreamZVtxMultContainer+;
#pragma link C++ class AliFemtoDreamPartCollection+;
#pragma link C++ class AliFemtoDreamAnalysis+;
#pragma link C++ class AliAnalysisTaskFemtoDream+;
#pragma link C++ class AliAnalysisTaskFemtoDreamDeuteron+;
#pragma link C++ class AliAnalysisTaskNanoXioton+;
#pragma link C++ class AliAnalysisTaskAODXioton+;
#pragma link C++ class AliAnalysisTaskNanoLoton+;
#pragma link C++ class AliAnalysisTaskAODLoton+;
#pragma link C++ class AliAnalysisTaskLeuteronAOD+;
#pragma link C++ class AliAnalysisTaskLeuteronNanoAOD+;
#pragma link C++ class AliAnalysisTaskNanoBBar+;
#pragma link C++ class AliAnalysisTaskFemtoDreamPhi+;
#pragma link C++ class AliAnalysisTaskDeuteronProtonEfficiency+;
#pragma link C++ class AliAnalysisTaskNanoAODFemtoDreamPhi+;
#pragma link C++ class AliAnalysisTaskNanoAODFemtoDreamLambdaPhi + ;
#pragma link C++ class AliAnalysisTaskFemtoDreamRho+;
#pragma link C++ class AliAnalysisTaskFemtoDreamSigPi+;
#pragma link C++ class AliAnalysisTaskNanoAODFemtoDreamSigPi+;
#pragma link C++ class AliAnalysisTaskFemtoDreamPion+;
#pragma link C++ class AliAnalysisTaskNanoAODFemtoDreamPion+;
#pragma link C++ class AliAnalysisTaskNanoAODSigma0Femto+;
#pragma link C++ class AliAnalysisTaskNanoSigmaPlus+;
#pragma link C++ class AliAnalysisTaskNanoLD+;
#pragma link C++ class AliAnalysisTaskLD+;
#pragma link C++ class AliAnalysisTaskGrandma+;
#pragma link C++ class AliAnalysisTaskOtonOmega+;
#pragma link C++ class AliOtonOmegaAnalysis+;
#pragma link C++ class AliOtonOmegaCascadeCuts+;
#pragma link C++ class AliOtonOmegaCascade+;
#pragma link C++ class AliFemtoDreamHigherPairMath+;
#pragma link C++ class AliAnalysisTaskNanoLX+;
#pragma link C++ class AliAnalysisTaskOtonOmegaNanoAOD+;
#pragma link C++ class AliAnalysisTaskOtonkd+;
#pragma link C++ class AliAnalysisTaskOtonkdAOD+;
#pragma link C++ class AliAnalysisTaskGeorgiosNTuple+;
#pragma link C++ class AliAnalysisTaskNanoPt+;
#pragma link C++ class AliAnalysisTaskPOmegaPenne+;
#pragma link C++ class AliFemtoDreamDump+;
#pragma link C++ class AliFemtoDreamEventDump+;
#pragma link C++ class AliFemtoDreamPairDump+;
#pragma link C++ class AliAnalysisTaskThreeBodyFemto+;
#pragma link C++ class AliAnalysisTaskNanoPPCoalescence+;
#pragma link C++ class AliAnalysisTaskFemtoSMI+;
#pragma link C++ class AliAnalysisTaskThreeBodyFemtoAOD+;
#pragma link C++ class AliAnalysisTaskThreeBodyFemtoAODPionProton+;
#pragma link C++ class AliAnalysisTaskThreeBodyFemtoPionProton+;
#pragma link C++ class AliAnalysisTaskThreeBodyProtonPrimary+;
#pragma link C++ class AliAnalysisTaskNanoTreeLPhi + ;
#pragma link C++ class AliAnalysisTaskNanoLambdaKaon + ;
#pragma link C++ class AliAnalysisTaskNanoLKr + ;
#pragma link C++ class AliAnalysisTaskThreeBodyFemtoMixedCharge+;
#pragma link C++ class AliAnalysisTaskNanoBenchmark + ;
#pragma link C++ class AliAnalysisTaskNanoFemtoProtonPion+;
#pragma link C++ class AliAnalysisTaskFemtoProtonPion+;

#pragma link C++ class AliSigma0AODPhotonMotherCuts+;
#pragma link C++ class AliSigma0PhotonCuts+;
#pragma link C++ class AliAnalysisTaskNanoXiPi+;
#pragma link C++ class AliAnalysisTaskNanoFemtoProtonKaonPlus+;
#pragma link C++ class AliAnalysisTaskFemtoProtonKaonPlus+;

#pragma link C++ class AliAnalysisTaskOtonXx+;


#endif

#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

//the classes below work only in ROOT6
#pragma link C++ class AliAnalysisTaskCharmingFemto+;
#pragma link C++ class AliAnalysisTaskCharmingFemto_og+;
#pragma link C++ class CustomQueue<std::vector<TLorentzVector>>+;
#pragma link C++ std::vector<CustomQueue<std::vector<TLorentzVector>>>+;
#pragma link C++ class CutContainer+;
#pragma link C++ class AliEasyFemto+;
#pragma link C++ class AliAnalysisTaskPionDeuteronMC+;
#pragma link C++ class AliAnalysisTaskPionDeuteron+;

#endif
