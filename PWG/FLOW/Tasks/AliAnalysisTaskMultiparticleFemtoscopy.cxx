/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  * 
**************************************************************************/

 /******************************** 
 * femtoscopy with multiparticle *
 *           technology          * 
 *                               * 
 * author: Ante Bilandzic        * 
 *        (abilandzic@gmail.com) *
 ********************************/ 
  
#include "Riostream.h"
#include "TDatabasePDG.h"
#include "AliAnalysisTaskMultiparticleFemtoscopy.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliStack.h"
#include "TFile.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskMultiparticleFemtoscopy)

//================================================================================================================

AliAnalysisTaskMultiparticleFemtoscopy::AliAnalysisTaskMultiparticleFemtoscopy(const char *name, Bool_t useParticleWeights): 
 AliAnalysisTaskSE(name), 
 fHistList(NULL),
 fAnalysisType(NULL),
 fPIDResponse(NULL),
 fMaxNoGlobalTracksAOD(5), // TBI this is landmine
 fProcessBothKineAndReco(kFALSE),
 fProcessOnlyKine(kFALSE),
 fProcessOnlyReco(kFALSE),
 fRejectFakeTracks(kTRUE),
 fMC(NULL),
 // 1.) Control histograms:
 fControlHistogramsList(NULL),
 fControlHistogramsFlagsPro(NULL),
 fFillControlHistograms(kFALSE),
 fControlHistogramsEventList(NULL),
 fControlHistogramsEventFlagsPro(NULL),
 fFillControlHistogramsEvent(kFALSE),
 fGetNumberOfTracksHist(NULL),
 fGetNumberOfGlobalTracksHist(NULL),
 fGetNumberOfV0sHist(NULL),
 fGetNumberOfCascadesHist(NULL),
 fGetMagneticFieldHist(NULL),
 fGetEventTypeHist(NULL),
 fGetCentralityHist(NULL),
 fGetNContributorsHist(NULL),
 fGetChi2perNDFHist(NULL),
 fGetNDaughtersHist(NULL),
 fControlHistogramsNonIdentifiedParticlesList(NULL),
 fControlHistogramsNonIdentifiedParticlesFlagsPro(NULL),
 fFillControlHistogramsNonIdentifiedParticles(kFALSE),
 fChargeHist(NULL),
 fGetTPCNclsHist(NULL),
 fGetTPCsignalNHist(NULL),
 fGetITSNclsHist(NULL),
 fdEdxVsPtHist(NULL),
 fPtHist(NULL),
 fEtaHist(NULL),
 fPhiHist(NULL),
 fMassHist(NULL),
 fGetFilterMap(NULL),
 fGetPdgCode(NULL),
 fControlHistogramsNonIdentifiedParticlesFTSFList(NULL),
 fControlHistogramsNonIdentifiedParticlesFTSFFlagsPro(NULL),
 fFillControlHistogramsNonIdentifiedParticlesFTSF(kFALSE),
 fFilterBitFTSF(128),
 fChargeFTSFHist(NULL),
 fGetTPCNclsFTSFHist(NULL),
 fGetTPCsignalNFTSFHist(NULL),
 fGetITSNclsFTSFHist(NULL),
 fdEdxVsPtFTSFHist(NULL),
 fPtFTSFHist(NULL),
 fEtaFTSFHist(NULL),
 fPhiFTSFHist(NULL),
 fMassFTSFHist(NULL),
 fGetFilterMapFTSF(NULL),
 fGetPdgCodeFTSF(NULL),
 fControlHistogramsIdentifiedParticlesList(NULL),
 fControlHistogramsIdentifiedParticlesFlagsPro(NULL),
 fFillControlHistogramsIdentifiedParticles(kFALSE),
 fFillControlHistogramsWithGlobalTrackInfo(kFALSE),
 fInclusiveSigmaCutsPro(NULL),
 fExclusiveSigmaCutsPro(NULL),
 fUseDefaultInclusiveSigmaCuts(kTRUE),
 fUseDefaultExclusiveSigmaCuts(kTRUE),
 fControlHistogramsV0sList(NULL),
 fControlHistogramsV0sFlagsPro(NULL),
 fFillControlHistogramsV0s(kFALSE),
 fGetNProngsHist(NULL),
 fMassK0ShortHist(NULL), 
 fMassLambdaHist(NULL),
 fMassAntiLambdaHist(NULL),
 fOpenAngleV0Hist(NULL),
 fRadiusV0Hist(NULL),
 fDcaV0ToPrimVertexHist(NULL),
 fMomV0XHist(NULL),
 fMomV0YHist(NULL),
 fMomV0ZHist(NULL),
 fPtV0Hist(NULL),
 fPseudoRapV0Hist(NULL),
 fPAHist(NULL),
 // 2.) Event-by-event histograms:
 fEBEHistogramsList(NULL),
 fEBEObjectsFlagsPro(NULL),
 //fFillEBEHistograms(kTRUE),
 fUniqueIDHistEBE(NULL), 
 // 3.) Correlation functions:
 fCorrelationFunctionsList(NULL),
 fCorrelationFunctionsFlagsPro(NULL),
 f2pCorrelationFunctionsFlagsPro(NULL),
 f3pCorrelationFunctionsFlagsPro(NULL),
 f4pCorrelationFunctionsFlagsPro(NULL),
 fFillCorrelationFunctions(kFALSE),
 fNormalizeCorrelationFunctions(kFALSE),
 fCorrelationFunctionsIndices(NULL),
 fFill3pCorrelationFunctions(kFALSE),
 fFill4pCorrelationFunctions(kFALSE),
 fNormalizationOption(0),
 fnMergedBins(-44),
 // 4.) Background:
 fBackgroundList(NULL),
 fBackgroundFlagsPro(NULL),
 f2pBackgroundFlagsPro(NULL),
 f3pBackgroundFlagsPro(NULL),
 f4pBackgroundFlagsPro(NULL),
 fBackgroundOption(0),
 fEstimate2pBackground(kFALSE),
 fEstimate3pBackground(kFALSE),
 fEstimate4pBackground(kFALSE),
 fMaxBufferSize1(10),
 // 5.) Buffers:
 fBuffersList(NULL),
 fBuffersFlagsPro(NULL),
 fFillBuffers(kFALSE),
 fMaxBuffer(-44),
 // 6.) QA:
 fQAList(NULL),
 fQAFlagsPro(NULL),
 fFillQA(kFALSE),
 fBailOutAfterQA(kFALSE),
 fFillQAEvents(kFALSE),
 fFillQAParticles(kFALSE),
 fQAEventsList(NULL),
 fQAParticlesList(NULL),
 fQAFilterBitScan(NULL),
 fQAIDvsFilterBit(NULL),
 // 7.) Common event cuts:
 fRejectEventsWithoutPrimaryVertex(kTRUE),
 fMinMagneticField(0.001),
 fCutOnNumberOfTracks(kFALSE),
 fMinNumberOfTracks(-44),
 fMaxNumberOfTracks(-44),
 fCutOnNumberOfGlobalTracks(kFALSE),
 fMinNumberOfGlobalTracks(-44),
 fMaxNumberOfGlobalTracks(-44),
 fCutOnNumberOfV0s(kFALSE),
 fMinNumberOfV0s(-44),
 fMaxNumberOfV0s(-44),
 fCutOnNumberOfCascades(kFALSE),
 fMinNumberOfCascades(-44),
 fMaxNumberOfCascades(-44),
 fCutOnVertexX(kFALSE),
 fMinVertexX(-44.),
 fMaxVertexX(-44.),
 fCutOnVertexY(kFALSE),
 fMinVertexY(-44.),
 fMaxVertexY(-44.),
 fCutOnVertexZ(kFALSE),
 fMinVertexZ(-44.),
 fMaxVertexZ(-44.),
 fCutOnNContributors(kFALSE),
 fMinNContributors(-44),
 fMaxNContributors(-44),
 // 8.) Global track cuts:
 fGlobalTrackCutsList(NULL),
 fGlobalTrackCutsFlagsPro(NULL),
 fApplyGlobalTrackCuts(kFALSE),
 // *.) Testing new ways to calculate correlation functions:
 fCorrelationFunctionsTESTList(NULL),
 fCorrelationFunctionsTESTFlagsPro(NULL),
 fBoost(kFALSE),
 fBoostVelocity(0.,0.,0.),
 fnQ2bins(100),
 fnQ2min(0.),
 fnQ2max(10.),
 fnQ3bins(100),
 fnQ3min(0.),
 fnQ3max(10.),
 // *.) Testing new ways to calculate background:
 fBackgroundTESTList(NULL),
 fBackgroundTESTFlagsPro(NULL),
 // *.) 'hybrid approach':
 fHybridApproachList(NULL),
 fHybridApproachFlagsPro(NULL),
 fDoHybridApproach(kFALSE),
 // *.) Online monitoring:
 fOnlineMonitoring(kFALSE),
 fUpdateOutputFile(kFALSE),
 fUpdateFrequency(-44),
 fUpdateWhichOutputFile(NULL),
 fMaxNumberOfEvents(-44),
 // *.) Debugging:
 fDoSomeDebugging(kFALSE),
 fWaitForSpecifiedEvent(kFALSE),
 fRun(0),
 fBunchCross(0),
 fOrbit(0),
 fPeriod(0)
 {
  // Constructor.
 
  AliDebug(2,"AliAnalysisTaskMultiparticleFemtoscopy::AliAnalysisTaskMultiparticleFemtoscopy(const char *name, Bool_t useParticleWeights)");

  // Base list:
  fHistList = new TList();
  fHistList->SetName("MPF"); // MultiParticle Femtoscopy
  fHistList->SetOwner(kTRUE);

  // Initialize all arrays:
  this->InitializeArrays();
  this->InitializeArraysForControlHistograms();
  this->InitializeArraysForEBEObjects();
  this->InitializeArraysForCorrelationFunctions();
  this->InitializeArraysForBackground();
  this->InitializeArraysForBuffers();
  this->InitializeArraysForQA();
  this->InitializeArraysForGlobalTrackCuts();
  this->InitializeArraysForCorrelationFunctionsTEST();
  this->InitializeArraysForBackgroundTEST();

  // Define input and output slots here
  // Input slot #0 works with an AliFlowEventSimple
  //DefineInput(0, AliFlowEventSimple::Class());  
  // Input slot #1 is needed for the weights input file:
  //if(useParticleWeights)
  //{
  // DefineInput(1, TList::Class());   
  //}  
  // Output slot #0 is reserved              
  // Output slot #1 writes into a TList container
  
  DefineOutput(1, TList::Class());  

  if(useParticleWeights)
  {
   // TBI add support eventually for particle weights
  }

} // AliAnalysisTaskMultiparticleFemtoscopy::AliAnalysisTaskMultiparticleFemtoscopy(const char *name, Bool_t useParticleWeights): 

//================================================================================================================

AliAnalysisTaskMultiparticleFemtoscopy::AliAnalysisTaskMultiparticleFemtoscopy(): 
 AliAnalysisTaskSE(),
 fHistList(NULL),
 fAnalysisType(NULL),
 fPIDResponse(NULL),
 fMaxNoGlobalTracksAOD(5), // TBI this is landmine
 fProcessBothKineAndReco(kFALSE),
 fProcessOnlyKine(kFALSE),
 fProcessOnlyReco(kFALSE),
 fRejectFakeTracks(kTRUE),
 fMC(NULL),
 // 1.) Control histograms:
 fControlHistogramsList(NULL),
 fControlHistogramsFlagsPro(NULL),
 fFillControlHistograms(kFALSE),
 fControlHistogramsEventList(NULL),
 fControlHistogramsEventFlagsPro(NULL),
 fFillControlHistogramsEvent(kFALSE),
 fGetNumberOfTracksHist(NULL),
 fGetNumberOfGlobalTracksHist(NULL),
 fGetNumberOfV0sHist(NULL),
 fGetNumberOfCascadesHist(NULL),
 fGetMagneticFieldHist(NULL),
 fGetEventTypeHist(NULL),
 fGetCentralityHist(NULL),
 fGetNContributorsHist(NULL),
 fGetChi2perNDFHist(NULL),
 fGetNDaughtersHist(NULL),
 fControlHistogramsNonIdentifiedParticlesList(NULL),
 fControlHistogramsNonIdentifiedParticlesFlagsPro(NULL),
 fFillControlHistogramsNonIdentifiedParticles(kFALSE),
 fChargeHist(NULL),
 fGetTPCNclsHist(NULL),
 fGetTPCsignalNHist(NULL),
 fGetITSNclsHist(NULL),
 fdEdxVsPtHist(NULL),
 fPtHist(NULL),
 fEtaHist(NULL),
 fPhiHist(NULL),
 fMassHist(NULL),
 fGetFilterMap(NULL),
 fGetPdgCode(NULL),
 fControlHistogramsNonIdentifiedParticlesFTSFList(NULL),
 fControlHistogramsNonIdentifiedParticlesFTSFFlagsPro(NULL),
 fFillControlHistogramsNonIdentifiedParticlesFTSF(kFALSE),
 fFilterBitFTSF(128),
 fChargeFTSFHist(NULL),
 fGetTPCNclsFTSFHist(NULL),
 fGetTPCsignalNFTSFHist(NULL),
 fGetITSNclsFTSFHist(NULL),
 fdEdxVsPtFTSFHist(NULL),
 fPtFTSFHist(NULL),
 fEtaFTSFHist(NULL),
 fPhiFTSFHist(NULL),
 fMassFTSFHist(NULL),
 fGetFilterMapFTSF(NULL),
 fGetPdgCodeFTSF(NULL),
 fControlHistogramsIdentifiedParticlesList(NULL),
 fControlHistogramsIdentifiedParticlesFlagsPro(NULL),
 fFillControlHistogramsIdentifiedParticles(kFALSE),
 fFillControlHistogramsWithGlobalTrackInfo(kFALSE),
 fInclusiveSigmaCutsPro(NULL),
 fExclusiveSigmaCutsPro(NULL),
 fUseDefaultInclusiveSigmaCuts(kFALSE),
 fUseDefaultExclusiveSigmaCuts(kFALSE),
 fControlHistogramsV0sList(NULL),
 fControlHistogramsV0sFlagsPro(NULL),
 fFillControlHistogramsV0s(kFALSE),
 fGetNProngsHist(NULL),
 fMassK0ShortHist(NULL), 
 fMassLambdaHist(NULL),
 fMassAntiLambdaHist(NULL),
 fOpenAngleV0Hist(NULL),
 fRadiusV0Hist(NULL),
 fDcaV0ToPrimVertexHist(NULL), 
 fMomV0XHist(NULL),
 fMomV0YHist(NULL),
 fMomV0ZHist(NULL),
 fPtV0Hist(NULL),
 fPseudoRapV0Hist(NULL),
 fPAHist(NULL),
 // 2.) Event-by-event histograms:
 fEBEHistogramsList(NULL),
 fEBEObjectsFlagsPro(NULL),
 //fFillEBEHistograms(kFALSE),
 fUniqueIDHistEBE(NULL),
 // 3.) Correlation functions:
 fCorrelationFunctionsList(NULL),
 fCorrelationFunctionsFlagsPro(NULL),
 f2pCorrelationFunctionsFlagsPro(NULL),
 f3pCorrelationFunctionsFlagsPro(NULL),
 f4pCorrelationFunctionsFlagsPro(NULL),
 fFillCorrelationFunctions(kFALSE),
 fNormalizeCorrelationFunctions(kFALSE),
 fCorrelationFunctionsIndices(NULL),
 fFill3pCorrelationFunctions(kFALSE),
 fFill4pCorrelationFunctions(kFALSE),
 fNormalizationOption(0),
 fnMergedBins(-44),
 // 4.) Background:
 fBackgroundList(NULL),
 fBackgroundFlagsPro(NULL),
 f2pBackgroundFlagsPro(NULL),
 f3pBackgroundFlagsPro(NULL),
 f4pBackgroundFlagsPro(NULL),
 fBackgroundOption(0),
 fEstimate2pBackground(kFALSE),
 fEstimate3pBackground(kFALSE),
 fEstimate4pBackground(kFALSE),
 fMaxBufferSize1(0),
 // 5.) Buffers:
 fBuffersList(NULL),
 fBuffersFlagsPro(NULL),
 fFillBuffers(kFALSE),
 fMaxBuffer(-44),
 // 6.) QA:
 fQAList(NULL),
 fQAFlagsPro(NULL),
 fFillQA(kFALSE),
 fBailOutAfterQA(kFALSE),
 fFillQAEvents(kFALSE),
 fFillQAParticles(kFALSE),
 fQAEventsList(NULL),
 fQAParticlesList(NULL),
 fQAFilterBitScan(NULL),
 fQAIDvsFilterBit(NULL),
 // 7.) Common event cuts:
 fRejectEventsWithoutPrimaryVertex(kFALSE),
 fMinMagneticField(0.001),
 fCutOnNumberOfTracks(kFALSE),
 fMinNumberOfTracks(-44),
 fMaxNumberOfTracks(-44),
 fCutOnNumberOfGlobalTracks(kFALSE),
 fMinNumberOfGlobalTracks(-44),
 fMaxNumberOfGlobalTracks(-44),
 fCutOnNumberOfV0s(kFALSE),
 fMinNumberOfV0s(-44),
 fMaxNumberOfV0s(-44),
 fCutOnNumberOfCascades(kFALSE),
 fMinNumberOfCascades(-44),
 fMaxNumberOfCascades(-44),
 fCutOnVertexX(kFALSE),
 fMinVertexX(-44.),
 fMaxVertexX(-44.),
 fCutOnVertexY(kFALSE),
 fMinVertexY(-44.),
 fMaxVertexY(-44.),
 fCutOnVertexZ(kFALSE),
 fMinVertexZ(-44.),
 fMaxVertexZ(-44.),
 fCutOnNContributors(kFALSE),
 fMinNContributors(-44),
 fMaxNContributors(-44),
 // 8.) Global track cuts:
 fGlobalTrackCutsList(NULL),
 fGlobalTrackCutsFlagsPro(NULL),
 fApplyGlobalTrackCuts(kFALSE),
 // *.) Testing new ways to calculate correlation functions:
 fCorrelationFunctionsTESTList(NULL),
 fCorrelationFunctionsTESTFlagsPro(NULL),
 fBoost(kFALSE),
 fBoostVelocity(0.,0.,0.),
 fnQ2bins(-44),
 fnQ2min(-44.),
 fnQ2max(-44.),
 fnQ3bins(-44),
 fnQ3min(-44.),
 fnQ3max(-44.),
 // *.) Testing new ways to calculate background:
 fBackgroundTESTList(NULL),
 fBackgroundTESTFlagsPro(NULL),
 // *.) 'hybrid approach':
 fHybridApproachList(NULL),
 fHybridApproachFlagsPro(NULL),
 fDoHybridApproach(kFALSE),
 // *.) Online monitoring:
 fOnlineMonitoring(kFALSE),
 fUpdateOutputFile(kFALSE),
 fUpdateFrequency(-44),
 fUpdateWhichOutputFile(NULL),
 fMaxNumberOfEvents(-44),
 // *.) Debugging:
 fDoSomeDebugging(kFALSE),
 fWaitForSpecifiedEvent(kFALSE),
 fRun(0),
 fBunchCross(0),
 fOrbit(0),
 fPeriod(0)
{
  // Dummy constructor.
 
  AliDebug(2,"AliAnalysisTaskMultiparticleFemtoscopy::AliAnalysisTaskMultiparticleFemtoscopy()");

} // AliAnalysisTaskMultiparticleFemtoscopy::AliAnalysisTaskMultiparticleFemtoscopy():

//================================================================================================================

AliAnalysisTaskMultiparticleFemtoscopy::~AliAnalysisTaskMultiparticleFemtoscopy()
{
 // Destructor.

 if(fHistList) delete fHistList;
 if(fPIDResponse) delete fPIDResponse;

 for(Int_t index=0;index<fMaxNoGlobalTracksAOD;index++)
 {
  if(fGlobalTracksAOD[index]) delete fGlobalTracksAOD[index];
 }

} // AliAnalysisTaskMultiparticleFemtoscopy::~AliAnalysisTaskMultiparticleFemtoscopy()

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::UserCreateOutputObjects() 
{
 // Called at every worker node to initialize.

 // a) Insanity checks;
 // b) Trick to avoid name clashes, part 1;
 // c) Set all flags which are not directly set via setters;
 // d) Book and nest all lists;
 // e) Book all objects;
 // f) Trick to avoid name clashes, part 2.

 // a) Insanity checks:
 this->InsanityChecksUserCreateOutputObjects();

 // b) Trick to avoid name clashes, part 1:
 Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
 TH1::AddDirectory(kFALSE);

 // c) Set all flags which are not directly set via setters:
 fFillQA = fFillQAEvents || fFillQAParticles;
 fFillControlHistograms = fFillControlHistogramsEvent || fFillControlHistogramsNonIdentifiedParticles || fFillControlHistogramsNonIdentifiedParticlesFTSF || fFillControlHistogramsIdentifiedParticles || fFillControlHistogramsV0s;

 // d) Book and nest all lists:
 this->BookAndNestAllLists();

 // e) Book all objects:
 this->BookEverything();
 this->BookEverythingForControlHistograms();
 this->BookEverythingForEBEObjects();
 this->BookEverythingForCorrelationFunctions();
 this->BookEverythingForBackground();
 this->BookEverythingForBuffers();
 this->BookEverythingForQA();
 this->BookEverythingForGlobalTrackCuts();
 this->BookEverythingForCorrelationFunctionsTEST();
 this->BookEverythingForBackgroundTEST();
 this->BookEverythingForHybridApproach();

 // f) Trick to avoid name clashes, part 2:
 TH1::AddDirectory(oldHistAddStatus);

 PostData(1,fHistList);

} // void AliAnalysisTaskMultiparticleFemtoscopy::UserCreateOutputObjects() 

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::UserExec(Option_t *) 
{
 // Main loop (called for each event).

 // a) Determine Ali{MC,ESD,AOD}Event and the analysis type;
 // b) Insanity checks; // TBI disabled at the moment (see implementation)
 // c) QA;
 // d) Start analysis over MC, AOD, ESD, MC_AOD or MC_ESD;
 // e) Reset event-by-event objects;
 // f) PostData;
 // g) Online monitoring.

 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::UserExec(Option_t *)";

 // a) Determine Ali{MC,ESD,AOD}Event and the analysis type:
 AliMCEvent *aMC = MCEvent();                                  // from TaskSE
 AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(InputEvent()); // from TaskSE
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(InputEvent()); // from TaskSE
 if(1 != (Int_t)fProcessBothKineAndReco+(Int_t)fProcessOnlyKine+(Int_t)fProcessOnlyReco)
 {
  Fatal(sMethodName.Data(),"One (and only one!) of fProcessBothKineAndReco, fProcessOnlyKine and fProcessOnlyReco must be kTRUE in AddTask* macro!!!!");
 }
 else if(fProcessBothKineAndReco)
 {
  if(aMC && aAOD){*fAnalysisType="MC_AOD";}
  else if(aMC && aESD){*fAnalysisType="MC_ESD";}
  else{Fatal(sMethodName.Data(),"fProcessBothKineAndReco is kTRUE, but TBI ...");}
 } // if(fProcessBothKineAndReco)
 else if(fProcessOnlyKine)
 {
  if(aMC){*fAnalysisType="MC";}
  else{Fatal(sMethodName.Data(),"fProcessOnlyKine is kTRUE, but TBI ...");}
 }
 else if(fProcessOnlyReco)
 {
  if(aAOD){*fAnalysisType="AOD";}
  else if(aESD){*fAnalysisType="ESD";}
  else{Fatal(sMethodName.Data(),"fProcessOnlyReco is kTRUE, but TBI ...");}
 }
 //cout<<Form("\n=> AATMF::UserExec(): fAnalysisType = \"%s\"\n",fAnalysisType->Data())<<endl;

 // b) Insanity checks:
 InsanityChecksUserExec();

 // c) QA:
 if(fFillQA)
 {
  if(fAnalysisType->EqualTo("MC")&&aMC){QA(aMC);}
  else if(fAnalysisType->EqualTo("AOD")&&aAOD){QA(aAOD);}
  else if(fAnalysisType->EqualTo("ESD")&&aESD){QA(aESD);}
  // TBI Do I have to support also cases "MC_AOD" and"MC_ESD"?
 } // if(fFillQA)
 if(fBailOutAfterQA){return;} // by default, this flag is set to kFALSE

 // d) Start analysis over MC, AOD, ESD, MC_AOD or MC_ESD:
 if(fAnalysisType->Contains("MC"))
 {
  fMC = aMC;
 }
 if(fAnalysisType->EqualTo("MC"))
 {
  if(aMC){MC(aMC);}
 }
 else if(fAnalysisType->EqualTo("AOD"))
 {
  if(aAOD){AOD(aAOD);}
 }
 else if(fAnalysisType->EqualTo("ESD"))
 {
  if(aESD){ESD(aESD);}
 }
 else if(fAnalysisType->EqualTo("MC_AOD"))
 {
  if(aMC&&aAOD){AOD(aAOD);}
 }
 else if(fAnalysisType->EqualTo("MC_ESD"))
 {
  if(aMC&&aESD){ESD(aESD);}
 }

 // e) Reset event-by-event objects:
 this->ResetEBEObjects();

 // f) PostData:
 PostData(1,fHistList);

 // g) Online monitoring:
 if(fOnlineMonitoring){OnlineMonitoring();}

} // void AliAnalysisTaskMultiparticleFemtoscopy::UserExec(Option_t *)

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::Terminate(Option_t *)
{
 // At the end of the day...

 // a) Get pointer to the grandmother of all lists "fHistList";
 // b) Get pointers to all objects in "fHistList";
 // c) Normalize correlation functions;
 // d) Dump the results.

 // a) Get pointer to the grandmother of all lists:
 fHistList = (TList*)GetOutputData(1); // TBI doesn't work locally :'(
 if(!fHistList){exit(1);}

 // b) Get pointers to all objects in "fHistList":
 this->GetOutputHistograms(fHistList);

 // c) Normalize correlation functions:
 if(fNormalizeCorrelationFunctions){this->NormalizeCorrelationFunctions();}

 // d) Dump the results:
 //TDirectoryFile *df = new TDirectoryFile("outputMPFanalysis","");
 //df->Add(fHistList);
 TFile *f = new TFile("AnalysisResults.root","RECREATE");
 fHistList->Write(fHistList->GetName(),TObject::kSingleKey);

 delete f;

} // end of void AliAnalysisTaskMultiparticleFemtoscopy::Terminate(Option_t *)

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::QA(AliVEvent *ave)
{
 // Local Quality Assurance. TBI

 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::QA(AliVEvent *ave)";

 // a) Determine Ali{MC,ESD,AOD}Event;

 // a) Determine Ali{MC,ESD,AOD}Event:
 AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 if(aMC)
 {
  // TBI
 }
 else if(aESD)
 {
  // TBI
 }
 else if(aAOD)
 {
  // a) Loop over 'aTracks'
  // b) Cut directly on aTracks:
  //    1.) e.g. is it TPC-only => via filter bits
  //    2.) kinematics => implement a function which from aTrack access its kinematics, and it cuts on it directly
  //    3.) PID => implement a function which from aTrack access its gTracks, and from it the PID, and then it cuts on it directly
  //    4.) global track parameters => implement a function which from aTrack access its gTracks, and then it cuts on its parameters directly
  // c)


  // a) Loop over 'aTracks':
  Int_t nTracks = aAOD->GetNumberOfTracks();
  for(Int_t iTrack=0;iTrack<nTracks;iTrack++)
  {
   AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack));
   Int_t nFilterBits = 14;
   for(Int_t fb=0;fb<nFilterBits;fb++) // 'fb' is a 'left shifting opetator', i.e. 1<<fb = filterbit
   {
    if(atrack->TestFilterBit(1<<fb))
    {
     fQAFilterBitScan->Fill(fb);
     fQAIDvsFilterBit->Fill(fb,atrack->GetID());
    } // if(atrack->TestFilterBit(1<<fb))
   } // for(Int_t fb=1;fb<nFilterBits;fb++)



   // t.b.c.



  } // for(Int_t iTrack=0;iTrack<nTracks;iTrack++)







/*
A) PARTICLE CUTS:

- Assuming two sets of histograms for track cuts:
    - fQAParticleHist[0][distribution_index][cut_index] => "before rain"
    - fQAParticleHist[1][distribution_index][cut_index] => "after rain"
- Used in conjunction with four set of cuts:
    - fParticleCuts[0] => direct cuts on aTrack
    - fParticleCuts[1] => given aTrack, these are the cuts on corresponding gTrack
    - fParticleCutsQA[0] => QA direct cuts on aTrack
    - fParticleCutsQA[1] => given aTrack, these are the QA cuts on corresponding gTrack
- Example 1:
    - Given aTrack
        - fParticleCuts[0] && fParticleCuts[1] => Fill QA histos "before"
        - fParticleCutsQA[0] && fParticleCutsQA[1] => Fill QA histos "after"


B) EVENT CUTS:

- Assuming two sets of histograms for event cuts:
    - fQAEventHist[0][distribution_index][cut_number] => "before rain"
    - fQAEventHist[1][distribution_index][cut_number] => "after rain"
- Used in conjunction with two set of cuts:
    - fEventCuts[0] => "default" cuts on event variables
    - fEventCuts[1] => given the "default" cuts on event variables, cut event further

- Example 1:
    - Given aAOD
        - if fEventCuts[0] => Fill QA histos "before"
        - if fEventCuts[1] => Fill QA histos "after"

*/


 }

 return;

} // void AliAnalysisTaskMultiparticleFemtoscopy::QA(AliVEvent *ave)

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::NormalizeCorrelationFunctions()
{
 // Normalize correlation functions with the background.

 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::NormalizeCorrelationFunctions()";

 // Say something nice to the audience...
 cout<<endl;
 cout<<"=> Normalization option:"<<endl;
 switch(fNormalizationOption)
 {
  case 0:
   cout<<"\"just scale\""<<endl;
  break;

  case 1:
   cout<<Form("   \"use concrete interval\": min = %f, max = %f",fNormalizationInterval[0],fNormalizationInterval[1])<<endl;
   cout<<"TBI: not implemented yet"<<endl;
   exit(0);
  break;

  default:
   cout<<Form("And the fatal 'fNormalizationOption' value is... %d. Congratulations!!",fNormalizationOption)<<endl;
   Fatal(sMethodName.Data(),"switch(fNormalizationOption)");
 } // switch(fNormalizationOption)
 if(fnMergedBins!=-44)
 {
  cout<<Form("=> Histograms will be rebinned: fnMergedBins = %d",fnMergedBins)<<endl;
 }
 cout<<endl;

 cout<<"TBI: not finalized yet"<<endl;

 for(Int_t pid1=0;pid1<10;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pid2=0;pid2<10;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   if(fCorrelationFunctions[pid1][pid2] && fCorrelationFunctions[pid1][pid2]->GetEntries() > 0 &&
      f2pBackground[pid1][pid2] && f2pBackground[pid1][pid2]->GetEntries() > 0)
   {
    fCorrelationFunctions[pid1][pid2]->Divide(f2pBackground[pid1][pid2]);
   } // if(..)
  } // for(Int_t pid2=0;pid2<10;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 } // for(Int_t pid1=0;pid1<10;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]


 for(Int_t pid1=0;pid1<10;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pid2=0;pid2<10;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   for(Int_t pid3=0;pid3<10;pid3++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
   {

    if(f3pCorrelationFunctions[pid1][pid2][pid3] && f3pCorrelationFunctions[pid1][pid2][pid3]->GetEntries() > 0 &&
      f3pBackground[pid1][pid2][pid3] && f3pBackground[pid1][pid2][pid3]->GetEntries() > 0)
    {
     f3pCorrelationFunctions[pid1][pid2][pid3]->Divide(f3pBackground[pid1][pid2][pid3]);
    } // if(..)

   } // for(Int_t pid3=0;pid3<10;pid3++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  } // for(Int_t pid2=0;pid2<10;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 } // for(Int_t pid1=0;pid1<10;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]


 for(Int_t pid1=0;pid1<10;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pid2=0;pid2<10;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   for(Int_t pid3=0;pid3<10;pid3++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
   {
    for(Int_t pid4=0;pid4<10;pid4++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
    {

     if(f4pCorrelationFunctions[pid1][pid2][pid3][pid4] && f4pCorrelationFunctions[pid1][pid2][pid3][pid4]->GetEntries() > 0 &&
       f4pBackground[pid1][pid2][pid3][pid4] && f4pBackground[pid1][pid2][pid3][pid4]->GetEntries() > 0)
     {
      f4pCorrelationFunctions[pid1][pid2][pid3][pid4]->Divide(f4pBackground[pid1][pid2][pid3][pid4]);
     } // if(..)
    } // for(Int_t pid4=0;pid4<10;pid4++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
   } // for(Int_t pid3=0;pid3<10;pid3++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  } // for(Int_t pid2=0;pid2<10;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 } // for(Int_t pid1=0;pid1<10;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]


} // void AliAnalysisTaskMultiparticleFemtoscopy::NormalizeCorrelationFunctions()

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsEvent(AliVEvent *ave)
{
 // Fill control histograms for global event observables.

 // TBI: not really validated, synchronize all if statements eventually with the analogous ones in UserExec(...)

 // To do:
 // 1) Add support for MC and ESD.

 // a) Determine Ali{MC,ESD,AOD}Event;
 // b) ...

 // a) Determine Ali{MC,ESD,AOD}Event:
 AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 if(aMC && fProcessOnlyKine) // TBI not validated
 {
  // MC event:
  fGetNumberOfTracksHist->Fill(aMC->GetNumberOfTracks());
  // ... TBI fill the rest
 } // if(aMC)
 else if(aESD)
 {
  // TBI
 } // else if(aESD)
 else if(aAOD) // TBI not validated, nevertheless works at the moment only for fProcessBothKineAndReco = kTRUE and fProcessOnlyReco = kTRUE
 {
  // AOD event:
  fGetNumberOfTracksHist->Fill(aAOD->GetNumberOfTracks()); // not all tracks are unique, see comments in function GlobalTracksAOD(AliAODEvent *aAOD, Int_t index)
  fGetNumberOfGlobalTracksHist->Fill(fGlobalTracksAOD[0]->GetSize()); // this is then my multiplicity. It shall be OK, since in fGlobalTracksAOD[0] I have filled only tracks with ID >= 0
  fGetNumberOfV0sHist->Fill(aAOD->GetNumberOfV0s()); // some V0s share the daughter, in this histogram they are not separated
  fGetNumberOfCascadesHist->Fill(aAOD->GetNumberOfCascades()); // TBI not validated
  fGetMagneticFieldHist->Fill(aAOD->GetMagneticField()); // TBI not validated
  fGetEventTypeHist->Fill(aAOD->GetEventType()); // TBI not validated
  if(aAOD->GetCentrality())
  {
   fGetCentralityHist->Fill(aAOD->GetCentrality()->GetCentralityPercentile("V0M")); // TBI not validated
  }
  // AOD primary vertex:
  AliAODVertex *avtx = (AliAODVertex*)aAOD->GetPrimaryVertex();
  fVertexXYZ[0]->Fill(avtx->GetX());
  fVertexXYZ[1]->Fill(avtx->GetY());
  fVertexXYZ[2]->Fill(avtx->GetZ());
  fGetNContributorsHist->Fill(avtx->GetNContributors()); // TBI not validated
  fGetChi2perNDFHist->Fill(avtx->GetChi2perNDF()); // TBI not validated
  fGetNDaughtersHist->Fill(avtx->GetNDaughters()); // TBI not validated
 } // else if(aAOD)
 
} // void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsEvent(AliVEvent *ave)

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsParticle(AliVEvent *ave)
{
 // Fill control histograms for particles.
 // Remark: The idea is to have one loop over the particles and to fill everything which needs to be filled.

 // To do:
 // a) Not really validated for fProcessOnlyKine = kTRUE
 // b) ESD not supported yet

 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsParticle(AliVEvent *ave)";

 // a) Determine Ali{MC,ESD,AOD}Event:
 AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 if(aMC && fProcessOnlyKine) // TBI not validated
 {
  // MC tracks:
  for(Int_t iTrack=0;iTrack<aMC->GetNumberOfTracks();iTrack++)
  {
   // a) Determine "amctrack";
   // b) Insanity checks for "amctrack";
   // c) Fill the control histograms for non-identified particles;
   // d) Fill the control histograms for identified particles.

   // a) Determine "amctrack":
   AliAODMCParticle *amcparticle = (AliAODMCParticle*)aMC->GetTrack(iTrack);

   // b) Insanity checks for "amctrack":
   if(!amcparticle){Fatal(sMethodName.Data(),"!amctrack");}

   // c) Fill the control histograms for non-identified particles:
   if(fFillControlHistogramsNonIdentifiedParticles){this->FillControlHistogramsNonIdentifiedParticles(amcparticle);}

   // d) Fill the control histograms for identified particles:
   if(fFillControlHistogramsIdentifiedParticles){this->FillControlHistogramsIdentifiedParticles(amcparticle);}

  } // for(Int_t iTrack=0;iTrack<aMC->GetNumberOfTracks();iTrack++)
 } // if(aMC)

 //=================================================================================================================

 else if(aESD)
 {
  // TBI
 } // else if(aESD)

 //=================================================================================================================

 else if(aAOD) // TBI: Works for fProcessBothKineAndReco = kTRUE and fProcessOnlyReco = kTRUE
 {
  // AOD tracks:
  for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++)
  {
   // a) Determine "atrack" (a.k.a. "any track in AOD");
   // b) Insanity checks for "atrack";
   // c) Determine the corresponding "gtrack" (a.k.a. "normal global" track);
   // d) Insanity checks for "gtrack";
   // e) Fill the control histograms for non-identified particles; TBI: not validated
   // f) Fill the control histograms for non-identified particles (f.t.s.f); // for the specific filterbit TBI: not validated
   // g) Fill the control histograms for identified particles (validated, works both for fProcessBothKineAndReco = kTRUE and fProcessOnlyReco = kTRUE, shall be used for purities)

   // a) Determine "atrack" (a.k.a. "any track in AOD"):
   AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack));

   // b) Insanity checks for "atrack":
   if(!atrack){Fatal(sMethodName.Data(),"!atrack");}
   Int_t n = InsanityChecksForTracks(atrack); // return values are documented in this function
   if(0!=n){Fatal(sMethodName.Data(),"InsanityChecksForTracks(atrack), n = %d",n);}

   // c) Determine the corresponding "gtrack" (a.k.a. "normal global" track):
   if(0 == fGlobalTracksAOD[0]->GetSize()){GlobalTracksAOD(aAOD,0);} // TBI most likely, this is an overkill
   if(0 == fGlobalTracksAOD[0]->GetSize()){return;}
   Int_t id = atrack->GetID();
   AliAODTrack *gtrack = dynamic_cast<AliAODTrack*>(id>=0 ? aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(id)) : aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(-(id+1))));

   // d) Insanity checks for "gtrack":
   if(!gtrack){Fatal(sMethodName.Data(),"!gtrack");}
   Int_t ng = InsanityChecksForGlobalTracks(gtrack); // return values are documented in this function
   if(0!=ng){Fatal(sMethodName.Data(),"InsanityChecksForGlobalTracks(gtrack), ng = %d",ng);}

   // e) Fill the control histograms for non-identified particles: TBI is there a more elegant solution to avoid double counting, besides checking for ID?
   if(fFillControlHistogramsNonIdentifiedParticles && atrack->GetID()>=0){this->FillControlHistogramsNonIdentifiedParticles(gtrack);} // for AOD, these are 'normal global' tracks

   // f) Fill the control histograms for non-identified particles (f.t.s.f): // for the specific filterbit TBI validated at the moment only for TPC-only, due to h.w. atrack->GetID()<0
   if(fFillControlHistogramsNonIdentifiedParticlesFTSF && atrack->GetID()<0 && atrack->TestFilterBit(fFilterBitFTSF)){this->FillControlHistogramsNonIdentifiedParticlesFTSF(atrack);}

   // g) Fill the control histograms for identified particles:
   if(fFillControlHistogramsIdentifiedParticles){this->FillControlHistogramsIdentifiedParticles(atrack,gtrack);}

  } // for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++)

 } // else if(aAOD)

} // void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsParticle(AliVEvent *ave)

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsNonIdentifiedParticles(AliAODTrack *gtrack)
{
 // Fill control histograms for non-identified particles.

 // a) Insanity checks;
 // b) Check out common selection criteria for 'normal global' tracks;
 // c) Fill control histograms.

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsNonIdentifiedParticles(AliAODTrack *gtrack)";
 if(!gtrack){Fatal(sMethodName.Data(),"!gtrack");} // TBI keep this for some time, eventually just continue

 // b) Check out common selection criteria for 'normal global' tracks:
 if(!PassesGlobalTrackCuts(gtrack)){return;} // TBI in the method PassesGlobalTrackCuts track is hardwired to AliAODTrack. Try to generalize

 // c) Fill control histograms:
 fChargeHist->Fill(gtrack->Charge()+0.5); // see how this histogram was booked for convention used
 fGetTPCNclsHist->Fill(gtrack->GetTPCNcls());
 fGetTPCsignalNHist->Fill(gtrack->GetTPCsignalN());
 fGetITSNclsHist->Fill(gtrack->GetITSNcls());
 fdEdxVsPtHist->Fill(gtrack->GetTPCmomentum(),gtrack->GetTPCsignal());
 fPtHist->Fill(gtrack->Pt());
 fEtaHist->Fill(gtrack->Eta());
 fPhiHist->Fill(gtrack->Phi());
 fMassHist->Fill(gtrack->M());

} // void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsNonIdentifiedParticles(AliAODTrack *gtrack)

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsNonIdentifiedParticles(AliAODMCParticle *amcparticle)
{
 // Fill control histograms for all charged MC particles.

 // TBI: re-validate

 // a) Insanity checks;
 // b) Check cut selection criteria;
 // c) Fill control histograms.

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsNonIdentifiedParticles(AliAODMCParticle *amcparticle)";
 if(!amcparticle){Fatal(sMethodName.Data(),"!amcparticle");}

 // b) Check cut selection criteria:
 //if(!PassesGlobalTrackCuts(gtrack)){return;} // TBI in the method PassesGlobalTrackCuts track is hardwired to AliAODTrack. Try to generalize

 // c) Fill control histograms:
 fChargeHist->Fill(amcparticle->Charge()+0.5); // see how this histogram was booked for convention used
 fPtHist->Fill(amcparticle->Pt());
 fEtaHist->Fill(amcparticle->Eta());
 fPhiHist->Fill(amcparticle->Phi());
 fMassHist->Fill(amcparticle->M());
 fGetPdgCode->Fill(amcparticle->GetPdgCode());

} // void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsNonIdentifiedParticles(AliAODMCParticle *amcparticle)

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsNonIdentifiedParticlesFTSF(AliAODTrack *atrack)
{
 // Fill control histograms for non identified particles, for the specified filterbit.

 // a) Insanity checks;
 // b) Check out the selection criteria for atracks;
 // c) Fill control histograms.

 // a) Insanity checks:
 TString sMethodName = "AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsNonIdentifiedParticlesFTSF(AliAODTrack *atrack)";
 if(!atrack){Fatal(sMethodName.Data(),"!atrack");} // TBI keep this for some time, eventually just continue

 // b) Check out the selection criteria for atracks:
 if(!PassesCommonTrackCuts(atrack)){return;}

 // c) Fill control histograms.
 fChargeFTSFHist->Fill(atrack->Charge()+0.5); // see how this histogram was booked for convention used
 fGetTPCNclsFTSFHist->Fill(atrack->GetTPCNcls());
 fGetTPCsignalNFTSFHist->Fill(atrack->GetTPCsignalN());
 fGetITSNclsFTSFHist->Fill(atrack->GetITSNcls());
 fdEdxVsPtFTSFHist->Fill(atrack->GetTPCmomentum(),atrack->GetTPCsignal());
 fPtFTSFHist->Fill(atrack->Pt());
 fEtaFTSFHist->Fill(atrack->Eta());
 fPhiFTSFHist->Fill(atrack->Phi());
 fMassFTSFHist->Fill(atrack->M());

} // void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsNonIdentifiedParticlesFTSF(AliAODTrack *atrack)

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsIdentifiedParticles(AliAODTrack *atrack, AliAODTrack *gtrack)
{
 // Fill control histograms for identified particles.

 // Remark 0: if fFillControlHistogramsWithGlobalTrackInfo = kTRUE track parameters are taken from 'gtrack', otherwise (and by default) from 'atrack'
 // Remark 1: if fFillControlHistogramsWithGlobalTrackInfo = kTRUE, only PassesGlobalTrackCuts(gtrack) is checked
 //           if fFillControlHistogramsWithGlobalTrackInfo = kFALSE, only PassesCommonTrackCuts(atrack) is checked.
 //           In both cases, PID functions Pion(...), etc. are called for 'gtrack', as only it has a PID info

 // a) Insanity checks;
 // b) Check cut selection criteria;
 // c) Fill control histograms.

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsIdentifiedParticles(AliAODTrack *atrack, AliAODTrack *gtrack)";
 if(!atrack){Fatal(sMethodName.Data(),"!atrack");} // TBI keep this for some time, eventually just continue
 if(!gtrack){Fatal(sMethodName.Data(),"!gtrack");} // TBI keep this for some time, eventually just continue

 // b) Check cut selection criteria:
 if(fFillControlHistogramsWithGlobalTrackInfo)
 {
  if(atrack->GetID()<0){return;} // this is needed to avoid double counting, I think. TBI well, re-think...
  if(!PassesGlobalTrackCuts(gtrack)){return;}
 }
 else
 {
  if(!PassesCommonTrackCuts(atrack)){return;}
 } // if(fFillControlHistogramsWithGlobalTrackInfo)

 // c) Fill control histograms:
 AliAODTrack *agtrack = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
 if(fFillControlHistogramsWithGlobalTrackInfo){agtrack = gtrack;}
 else {agtrack = atrack;}
 if(!agtrack){Fatal(sMethodName.Data(),"!agtrack");}
 Int_t nCharge = 1;
 Bool_t bPrimary = kTRUE;
 for(Int_t charge=0;charge<2;charge++) // 0 = +q, 1 = -q
 {
  if(1==charge){nCharge = -1;} // TBI this is just ugly
  for(Int_t ps=0;ps<2;ps++) // 0 = kPrimary, 1 = kFromDecayVtx
  {
   if(1==ps){bPrimary = kFALSE;}
   // c2) Pions:
   if(Pion(gtrack,nCharge,bPrimary))
   {
    fPtPIDHist[2][charge][ps]->Fill(agtrack->Pt());
    fPPIDHist[2][charge][ps][0]->Fill(agtrack->Px());
    fPPIDHist[2][charge][ps][1]->Fill(agtrack->Py());
    fPPIDHist[2][charge][ps][2]->Fill(agtrack->Pz());
    fMassPIDHist[2][charge][ps]->Fill(agtrack->M());
    fEtaPIDHist[2][charge][ps]->Fill(agtrack->Eta());
    fPhiPIDHist[2][charge][ps]->Fill(agtrack->Phi());
   }
   // c3) Kaons:
   if(Kaon(gtrack,nCharge,bPrimary))
   {
    fPtPIDHist[3][charge][ps]->Fill(agtrack->Pt());
    fPPIDHist[3][charge][ps][0]->Fill(agtrack->Px());
    fPPIDHist[3][charge][ps][1]->Fill(agtrack->Py());
    fPPIDHist[3][charge][ps][2]->Fill(agtrack->Pz());
    fMassPIDHist[3][charge][ps]->Fill(agtrack->M());
    fEtaPIDHist[3][charge][ps]->Fill(agtrack->Eta());
    fPhiPIDHist[3][charge][ps]->Fill(agtrack->Phi());
   }
   // c4) Protons:
   if(Proton(gtrack,nCharge,bPrimary))
   {
    fPtPIDHist[4][charge][ps]->Fill(agtrack->Pt());
    fPPIDHist[4][charge][ps][0]->Fill(agtrack->Px());
    fPPIDHist[4][charge][ps][1]->Fill(agtrack->Py());
    fPPIDHist[4][charge][ps][2]->Fill(agtrack->Pz());
    fMassPIDHist[4][charge][ps]->Fill(agtrack->M());
    fEtaPIDHist[4][charge][ps]->Fill(agtrack->Eta());
    fPhiPIDHist[4][charge][ps]->Fill(agtrack->Phi());
   }
  } // for(Int_t ps=0;ps<2;ps++) // 0 = kPrimary, 1 = kFromDecayVtx
 } // for(Int_t charge=0;charge<2;charge++) // 0 = +q, 1 = -q

} // void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsIdentifiedParticles(AliAODTrack *atrack, AliAODTrack *gtrack)

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsIdentifiedParticles(AliAODMCParticle *amcparticle)
{
 // Fill control histograms for identified Monte Carlo particles.

 // a) Insanity checks;
 // b) Check cut selection criteria;
 // c) Fill control histograms.

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsIdentifiedParticles(AliAODMCParticle *amcparticle)";
 if(!amcparticle){Fatal(sMethodName.Data(),"!amcparticle");}

 // b) Check cut selection criteria:
 if(!PassesCommonTrackCuts(amcparticle)){return;} // TBI have a look at implementation

 // c) Fill control histograms:
 Int_t index = -44;
 //  c2) Pions:
 if(211==TMath::Abs(amcparticle->GetPdgCode()))
 {
  index = 2;
 } // if(211==TMath::Abs(amcparticle->GetPdgCode()))
 //  c3) Kaons:
 else if(321==TMath::Abs(amcparticle->GetPdgCode()))
 {
  index = 3;
 } // if(321==TMath::Abs(amcparticle->GetPdgCode()))
 //  c4) Protons:
 else if(2212==TMath::Abs(amcparticle->GetPdgCode()))
 {
  index = 4;
 } // if(2212==TMath::Abs(amcparticle->GetPdgCode()))
 if(-44 != index)
 {
  Int_t charge = ((amcparticle->GetPdgCode()>0. ? 0 : 1)); // Okay... TBI
  Int_t isPhysicalPrimary = ((amcparticle->IsPhysicalPrimary() ? 0 : 1)); // Okay... TBI
  fPtPIDHist[index][charge][isPhysicalPrimary]->Fill(amcparticle->Pt());
  fPPIDHist[index][charge][isPhysicalPrimary][0]->Fill(amcparticle->Px());
  fPPIDHist[index][charge][isPhysicalPrimary][1]->Fill(amcparticle->Py());
  fPPIDHist[index][charge][isPhysicalPrimary][2]->Fill(amcparticle->Pz());
  fMassPIDHist[index][charge][isPhysicalPrimary]->Fill(amcparticle->M()); // in this context, clearly this is just a control histogram
  fEtaPIDHist[index][charge][isPhysicalPrimary]->Fill(amcparticle->Eta());
  fPhiPIDHist[index][charge][isPhysicalPrimary]->Fill(amcparticle->Phi());
 } // if(-44 != index)

} // void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsIdentifiedParticles(AliAODMCParticle *amcparticle)

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::EstimateBackground(AliVEvent *ave)
{
 // Estimate background.

 // To do:
 // 1) Add support for MC and ESD.

 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::EstimateBackground(AliVEvent *ave)";

 // a) Determine Ali{MC,ESD,AOD}Event;
 // b) ...

 // a) Determine Ali{MC,ESD,AOD}Event:
 AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 if(aMC)
 {
  TString pattern = ".11.13.211.321.2212."; // TBI this can be done programatically TBI only for these particles I need to estimate background

  return; // TBI re-validate the lines below, within if statement, by following the analogy with AOD case

  if(!PassesMixedEventCuts(aMC)){return;} // TBI this is empty at the moment for MC
  Int_t nTracks = aMC->GetNumberOfTracks();
  if(!fMixedEvents0[0])
  {
   TClonesArray ca0("AliAODMCParticle"); // temporary holder, its entries will be copied to fMixedEvents0[0]
   Int_t ca0Counter = 0;
   for(Int_t iTrack=0;iTrack<nTracks;iTrack++)
   {
    AliAODMCParticle *amcparticle = dynamic_cast<AliAODMCParticle*>(aMC->GetTrack(iTrack));
    if(!amcparticle){Fatal(sMethodName.Data(),"!amcparticle");}
    if(!PassesCommonTrackCuts(amcparticle)){continue;} // TBI re-think, see implemntation of this method
    if(!(pattern.Contains(Form(".%d.",TMath::Abs(amcparticle->GetPdgCode()))))){continue;}
    // Fill fMixedEvents0[0] with tracks which survived all checks:
    ca0[ca0Counter++] = amcparticle;
   } // for(Int_t iTrack=0;iTrack<nTracks;iTrack++)
   fMixedEvents0[0] = (TClonesArray*)ca0.Clone();
   if(!fMixedEvents0[0]){Fatal(sMethodName.Data(),"!fMixedEvents0[0]");}
  }
  else if(!fMixedEvents0[1])
  {
   TClonesArray ca1("AliAODMCParticle"); // temporary holder, its entries will be copied to fMixedEvents0[1]
   Int_t ca1Counter = 0;
   for(Int_t iTrack=0;iTrack<nTracks;iTrack++)
   {
    AliAODMCParticle *amcparticle = dynamic_cast<AliAODMCParticle*>(aMC->GetTrack(iTrack));
    if(!amcparticle){Fatal(sMethodName.Data(),"!amcparticle");}
    if(!PassesCommonTrackCuts(amcparticle)){continue;} // TBI re-think, see implemntation of this method
    if(!(pattern.Contains(Form(".%d.",TMath::Abs(amcparticle->GetPdgCode()))))){continue;}
    // Fill fMixedEvents0[0] with tracks which survived all checks:
    ca1[ca1Counter++] = amcparticle;
   } // for(Int_t iTrack=0;iTrack<nTracks;iTrack++)
   fMixedEvents0[1] = (TClonesArray*)ca1.Clone();
   if(!fMixedEvents0[1]){Fatal(sMethodName.Data(),"!fMixedEvents0[1]");}
  }
  // Shall we do something?
  if(fMixedEvents0[0] && fMixedEvents0[1])
  {
   Calculate2pBackground(fMixedEvents0[0],fMixedEvents0[1],kTRUE);
  }
 } // if(aMC)
 else if(aESD)
 {
  // TBI
 } // else if(aESD)
 else if(aAOD)
 {
  if(!PassesMixedEventCuts(aAOD)){return;} // TBI returns always true at the moment...

  switch(fBackgroundOption)
  {
   // 0 : shifting
   // 1 : permutations + binning in vertex z-position
   //     Key objects are fMixedEvents1[10][50]. First index indicates the vertex-z range, when the all entries in the second index are filled, then:
   //     a) Make all possible permutations of 2-events, 3-events, etc.
   //     b) For each permutation call Calculate2pBackground(...), Calculate3pBackground(...), etc.

   case 0: // shifting
   if(0==fMixedEvents0[0]->GetEntries())
   {
    if(fGlobalTracksAOD[1]) fGlobalTracksAOD[1]->Delete();
    GlobalTracksAOD(aAOD,1);
    if(0 == fGlobalTracksAOD[1]->GetSize()){return;}
    *fMixedEvents0[0] = *(aAOD->GetTracks());
    if(0==fMixedEvents0[0]->GetEntries()){return;}
   }
   else if(0==fMixedEvents0[1]->GetEntries())
   {
    if(fGlobalTracksAOD[2]) fGlobalTracksAOD[2]->Delete();
    GlobalTracksAOD(aAOD,2);
    if(0 == fGlobalTracksAOD[2]->GetSize()){return;}
    *fMixedEvents0[1] = *(aAOD->GetTracks());
    if(0==fMixedEvents0[1]->GetEntries()){return;}
   }
   else if(0==fMixedEvents0[2]->GetEntries())
   {
    if(fGlobalTracksAOD[3]) fGlobalTracksAOD[3]->Delete();
    GlobalTracksAOD(aAOD,3);
    if(0 == fGlobalTracksAOD[3]->GetSize()){return;}
    *fMixedEvents0[2] = *(aAOD->GetTracks());
    if(0==fMixedEvents0[2]->GetEntries()){return;}
   }
   // Shall we do something?
   if(0!=fMixedEvents0[0]->GetEntries() && 0!=fMixedEvents0[1]->GetEntries() && 0!=fMixedEvents0[2]->GetEntries())
   {
    if(fEstimate2pBackground) Calculate2pBackground(fMixedEvents0[0],fMixedEvents0[1]); // TBI rename
    if(fEstimate3pBackground) Calculate3pBackground(fMixedEvents0[0],fMixedEvents0[1],fMixedEvents0[2]);
    //if(fEstimate4pBackground) Calculate4pBackground(fMixedEvents0[0],fMixedEvents0[1],fMixedEvents0[2],fMixedEvents0[3]);
    // Shift mixed events:
    // [1] -> [0]
    fMixedEvents0[0]->Delete();
    *fMixedEvents0[0] = *fMixedEvents0[1];
    *fGlobalTracksAOD[1] = *fGlobalTracksAOD[2];
    // [2] -> [1]
    fMixedEvents0[1]->Delete();
    *fMixedEvents0[1] = *fMixedEvents0[2];
    *fGlobalTracksAOD[2] = *fGlobalTracksAOD[3];
    // Clean [2]:
    fMixedEvents0[2]->Delete();
    fGlobalTracksAOD[3]->Delete();
   } // if(0!=fMixedEvents0[0]->GetEntries() && 0!=fMixedEvents0[1]->GetEntries() && 0!=fMixedEvents0[2]->GetEntries())
   break; // case 0: // shifting

   //=============================================================

   case 1: // permutations + binning in vertex z-position
   AliAODVertex *avtx = (AliAODVertex*)aAOD->GetPrimaryVertex();
   if(!avtx){return;}
   Int_t nBinsVertexZrange = sizeof(fMixedEvents1)/sizeof(fMixedEvents1[0]); // extracts numbers of rows in TClonesArray *fMixedEvents1[10][5];
   if(nBinsVertexZrange<=0){Fatal(sMethodName.Data(),"nBinsVertexZrange<=0");}
   Int_t indexX = -44; // this is eventually the first index in TClonesArray *fMixedEvents1[10][5];
   Float_t flBinWidthVertexZrange = (fMaxVertexZ-fMinVertexZ)/(1.*nBinsVertexZrange);
   for(Int_t b=0;b<nBinsVertexZrange;b++)
   {
    //cout<<Form("%f - %f : %f",fMinVertexZ+1.*b*flBinWidthVertexZrange,fMinVertexZ+1.*(b+1)*flBinWidthVertexZrange,avtx->GetZ())<<endl;
    if(avtx->GetZ()>=fMinVertexZ+1.*b*flBinWidthVertexZrange && avtx->GetZ()<fMinVertexZ+1.*(b+1)*flBinWidthVertexZrange)
    {
     indexX = b;
     break;
    }
   } // for(Int_t b=0;b<nBinsVertexZrange;b++)
   if(-44==indexX){Fatal(sMethodName.Data(),"-44==indexX");}
   Int_t indexY = -44; // this is eventually the second index in TClonesArray *fMixedEvents1[10][5];
   for(Int_t bufferNo=0;bufferNo<fMaxBufferSize1;bufferNo++) // determining the 1st empty buffer
   {
    if(0==fMixedEvents1[indexX][bufferNo]->GetEntries()) // empty buffer
    {
     indexY = bufferNo;
     break;
    }
   } // for(Int_t bufferNo=0;bufferNo<5;bufferNo++)
   // Save in the buffer this event:
   *fMixedEvents1[indexX][indexY] = *(aAOD->GetTracks());
   GlobalTracksAOD(aAOD,indexX,indexY);
   // If the buffer is full, calculate background:
   if(fMaxBufferSize1-1==indexY)
   {
    cout<<"Flushing the buffer for background..."<<endl;
    cout<<Form("fMaxBufferSize1 = %d",fMaxBufferSize1)<<endl;
    cout<<Form("fMixedEvents1[indexX][indexY] = fMixedEvents1[%d][%d]",indexX,indexY)<<endl;
    // Shall we do something for 2p background?
    if(fEstimate2pBackground) // this buffer is full TBI hardwired 4. Making all possible combinations of 2 events
    {
     for(Int_t bufferNo1=0;bufferNo1<fMaxBufferSize1;bufferNo1++) // 1st event in the mixed events
     {
      for(Int_t bufferNo2=bufferNo1+1;bufferNo2<fMaxBufferSize1;bufferNo2++) // 2nd event in the mixed events
      {
       Calculate2pBackground(fMixedEvents1[indexX][bufferNo1],fMixedEvents1[indexX][bufferNo2],fGlobalTracksAOD1[indexX][bufferNo1],fGlobalTracksAOD1[indexX][bufferNo2]);
      } // for(Int_t bufferNo2=bufferNo1+1;bufferNo2<5;bufferNo2++) // 2nd event in the mixed events
     } // for(Int_t bufferNo1=0;bufferNo1<5;bufferNo1++) // first event in the mixed events
    } // if(fEstimate2pBackground)

    // Shall we do something for 3p background?
    if(fEstimate3pBackground) // this buffer is full TBI hardwired 4. Making all possible combinations of 3 events
    {
     for(Int_t bufferNo1=0;bufferNo1<fMaxBufferSize1;bufferNo1++) // 1st event in the mixed events
     {
      for(Int_t bufferNo2=0;bufferNo2<fMaxBufferSize1;bufferNo2++) // 2nd event in the mixed events
      {
       if(bufferNo2<=bufferNo1){continue;}
       for(Int_t bufferNo3=0;bufferNo3<fMaxBufferSize1;bufferNo3++) // 3rd event in the mixed events
       {
        if(bufferNo2<=bufferNo1 || bufferNo3<=bufferNo1 || bufferNo3<=bufferNo2){continue;}
        //cout<<Form("%d %d %d",bufferNo1,bufferNo2,bufferNo3)<<endl;
        Calculate3pBackground(fMixedEvents1[indexX][bufferNo1],fMixedEvents1[indexX][bufferNo2],fMixedEvents1[indexX][bufferNo3],fGlobalTracksAOD1[indexX][bufferNo1],fGlobalTracksAOD1[indexX][bufferNo2],fGlobalTracksAOD1[indexX][bufferNo3]);
       } // for(Int_t bufferNo3=0;bufferNo3<fMaxBufferSize1;bufferNo3++) // 3rd event in the mixed events
      } // for(Int_t bufferNo2=bufferNo1+1;bufferNo2<5;bufferNo2++) // 2nd event in the mixed events
     } // for(Int_t bufferNo1=0;bufferNo1<5;bufferNo1++) // first event in the mixed events
    } // if(fEstimate3pBackground)

    // Shall we do something for 4p background?
    if(fEstimate4pBackground) // this buffer is full TBI hardwired 4. Making all possible combinations of 4 events
    {
     // ...
    } // if(fEstimate4pBackground)

    // Clean the buffer:
    for(Int_t bufferNo=0;bufferNo<fMaxBufferSize1;bufferNo++)
    {
     fMixedEvents1[indexX][bufferNo]->Delete();
     fGlobalTracksAOD1[indexX][bufferNo]->Delete();
    } // for(Int_t bufferNo=0;bufferNo<5;bufferNo++)

   } // if(fMaxBufferSize1-1==indexY)

   break; // permutations + binning in vertex z-position

   //=============================================================

/*
   default:
    cout<<Form("And the fatal 'fBackgroundOption' value is... %d. Congratulations!!",fBackgroundOption)<<endl;
    Fatal(sMethodName.Data(),"switch(fBackgroundOption)");
*/

  } // switch(fBackgroundOption)

 } // else if(aAOD)

} // void AliAnalysisTaskMultiparticleFemtoscopy::EstimateBackground(AliVEvent *ave)

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::EstimateBackgroundTEST(AliVEvent *ave)
{
 // Estimate background for TEST, using only 'shifting' at the moment.

 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::EstimateBackgroundTEST(AliVEvent *ave)";

 // a) Determine Ali{MC,ESD,AOD}Event;
 // b) ...

 // a) Determine Ali{MC,ESD,AOD}Event:
 AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 if(aMC)
 {
  // TBI
 } // if(aMC)
 else if(aESD)
 {
  // TBI
 } // else if(aESD)
 else if(aAOD)
 {
  if(!PassesMixedEventCuts(aAOD)){return;} // TBI returns always true at the moment...

   if(0==fMixedEventsTEST[0]->GetEntries())
   {
    if(fGlobalTracksAODTEST[0]) fGlobalTracksAODTEST[0]->Delete();
    GlobalTracksAODTEST(aAOD,0);
    if(0 == fGlobalTracksAODTEST[0]->GetSize()){return;}
    *fMixedEventsTEST[0] = *(aAOD->GetTracks());
    if(0==fMixedEventsTEST[0]->GetEntries()){return;}
   }
   else if(0==fMixedEventsTEST[1]->GetEntries())
   {
    if(fGlobalTracksAODTEST[1]) fGlobalTracksAODTEST[1]->Delete();
    GlobalTracksAODTEST(aAOD,1);
    if(0 == fGlobalTracksAODTEST[1]->GetSize()){return;}
    *fMixedEventsTEST[1] = *(aAOD->GetTracks());
    if(0==fMixedEventsTEST[1]->GetEntries()){return;}
   }
   else if(0==fMixedEventsTEST[2]->GetEntries())
   {
    if(fGlobalTracksAODTEST[2]) fGlobalTracksAODTEST[2]->Delete();
    GlobalTracksAODTEST(aAOD,2);
    if(0 == fGlobalTracksAODTEST[2]->GetSize()){return;}
    *fMixedEventsTEST[2] = *(aAOD->GetTracks());
    if(0==fMixedEventsTEST[2]->GetEntries()){return;}
   }
   // Shall we do something?
   if(0!=fMixedEventsTEST[0]->GetEntries() && 0!=fMixedEventsTEST[1]->GetEntries() && 0!=fMixedEventsTEST[2]->GetEntries())
   {

    if(fEstimate2pBackground) Calculate2pBackgroundTEST(fMixedEventsTEST[0],fMixedEventsTEST[1],fGlobalTracksAODTEST[0],fGlobalTracksAODTEST[1]); // TBI rename
    if(fEstimate3pBackground) Calculate3pBackgroundTEST(fMixedEventsTEST[0],fMixedEventsTEST[1],fMixedEventsTEST[2],fGlobalTracksAODTEST[0],fGlobalTracksAODTEST[1],fGlobalTracksAODTEST[2]);

    //if(fEstimate4pBackground) Calculate4pBackground(fMixedEventsTEST[0],fMixedEventsTEST[1],fMixedEventsTEST[2],fMixedEventsTEST[3]);
    // Shift mixed events:
    // [1] -> [0]
    fMixedEventsTEST[0]->Delete();
    *fMixedEventsTEST[0] = *fMixedEventsTEST[1];
    *fGlobalTracksAODTEST[0] = *fGlobalTracksAODTEST[1];
    // [2] -> [1]
    fMixedEventsTEST[1]->Delete();
    *fMixedEventsTEST[1] = *fMixedEventsTEST[2];
    *fGlobalTracksAODTEST[1] = *fGlobalTracksAODTEST[2];
    // Clean [2]:
    fMixedEventsTEST[2]->Delete();
    fGlobalTracksAODTEST[2]->Delete();
   } // if(0!=fMixedEventsTEST[0]->GetEntries() && 0!=fMixedEventsTEST[1]->GetEntries() && 0!=fMixedEventsTEST[2]->GetEntries())

 } // else if(aAOD)

} // void AliAnalysisTaskMultiparticleFemtoscopy::EstimateBackgroundTEST(AliVEvent *ave)

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::DoHybridApproach(AliVEvent *ave)
{
 // Estimate background for TEST, using only 'shifting' at the moment.

 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::DoHybridApproach(AliVEvent *ave)";

 // a) Determine Ali{MC,ESD,AOD}Event;
 // b) ...

 // a) Determine Ali{MC,ESD,AOD}Event:
 AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 if(aMC)
 {
  // TBI
 } // if(aMC)
 else if(aESD)
 {
  // TBI
 } // else if(aESD)
 else if(aAOD)
 {
  if(!PassesMixedEventCuts(aAOD)){return;} // TBI returns always true at the moment...

   // a) Calculate 1st term:
   HybridApproach1stTerm(aAOD); // calculate N3(p1,p2,p3)

   // b) Mixed events buffers:
   if(0==fMixedEventsHA[0]->GetEntries())
   {
    if(fGlobalTracksAODHA[0]) fGlobalTracksAODHA[0]->Delete();
    GlobalTracksAODHA(aAOD,0);
    if(0 == fGlobalTracksAODHA[0]->GetSize()){return;}
    *fMixedEventsHA[0] = *(aAOD->GetTracks());
    if(0==fMixedEventsHA[0]->GetEntries()){return;}
   }
   else if(0==fMixedEventsHA[1]->GetEntries())
   {
    if(fGlobalTracksAODHA[1]) fGlobalTracksAODHA[1]->Delete();
    GlobalTracksAODHA(aAOD,1);
    if(0 == fGlobalTracksAODHA[1]->GetSize()){return;}
    *fMixedEventsHA[1] = *(aAOD->GetTracks());
    if(0==fMixedEventsHA[1]->GetEntries()){return;}
   }
   else if(0==fMixedEventsHA[2]->GetEntries())
   {
    if(fGlobalTracksAODHA[2]) fGlobalTracksAODHA[2]->Delete();
    GlobalTracksAODHA(aAOD,2);
    if(0 == fGlobalTracksAODHA[2]->GetSize()){return;}
    *fMixedEventsHA[2] = *(aAOD->GetTracks());
    if(0==fMixedEventsHA[2]->GetEntries()){return;}
   }

   // c) If all three buffers for mixed events are full, do something for the remaining terms:
   if(0!=fMixedEventsHA[0]->GetEntries() && 0!=fMixedEventsHA[1]->GetEntries() && 0!=fMixedEventsHA[2]->GetEntries())
   {
    // d) Cross-terms, e.g. N_{2}^A(p1,p2)N^B_{1}(p3). TBI Since I am focused in this exercise only on 3 identical pions, because of symmetry I need to calculate only one cross-term
    //    Remark: By design, the current event can go only to the last buffer ([2]), so it's safe to mix it both with [0] and [1], to gain statistics
    HybridApproach2ndTerm(aAOD,fMixedEventsHA[0],fGlobalTracksAODHA[0]); // calculate N_{2}^A(p1,p2)N^B_{1}(p3), where B is [0]
    HybridApproach2ndTerm(aAOD,fMixedEventsHA[1],fGlobalTracksAODHA[1]); // calculate N_{2}^A(p1,p2)N^B_{1}(p3), where B is [1]

    // e) Genuine mixed-event terms:
    HybridApproach5thTerm(fMixedEventsHA[0],fMixedEventsHA[1],fMixedEventsHA[2],fGlobalTracksAODHA[0],fGlobalTracksAODHA[1],fGlobalTracksAODHA[2]); // calculate N^A_{1}(p1)N^B_{1}(p2)N^C_{1}(p3)

    // Shift mixed events:
    // [1] -> [0]
    fMixedEventsHA[0]->Delete();
    *fMixedEventsHA[0] = *fMixedEventsHA[1];
    *fGlobalTracksAODHA[0] = *fGlobalTracksAODHA[1];
    // [2] -> [1]
    fMixedEventsHA[1]->Delete();
    *fMixedEventsHA[1] = *fMixedEventsHA[2];
    *fGlobalTracksAODHA[1] = *fGlobalTracksAODHA[2];
    // Clean [2]:
    fMixedEventsHA[2]->Delete();
    fGlobalTracksAODHA[2]->Delete();

   } // if(0!=fMixedEventsHA[0]->GetEntries() && 0!=fMixedEventsHA[1]->GetEntries() && 0!=fMixedEventsHA[2]->GetEntries())

 } // else if(aAOD)

} // void AliAnalysisTaskMultiparticleFemtoscopy::DoHybridApproach(AliVEvent *ave)

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::DoSomeDebugging(AliVEvent *ave)
{
 // Do all debugging in this function.

 // TBI: add support for fProcessBothKineAndReco

 // a) Determine Ali{MC,ESD,AOD}Event;
 // b) Wait for specified event.

 // a) Determine Ali{MC,ESD,AOD}Event:
 AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 // b) Wait for specified event:
 if(fWaitForSpecifiedEvent)
 {
  if(aMC)
  {
   // TBI
  } // if(aMC)
  else if(aESD)
  {
   // TBI
  } // else if(aESD)
  else if(aAOD)
  {
   cout<<Form("aAOD->GetRunNumber() = %d",aAOD->GetRunNumber())<<endl;
   cout<<Form("aAOD->GetBunchCrossNumber() = %d",aAOD->GetBunchCrossNumber())<<endl;
   cout<<Form("aAOD->GetOrbitNumber() = %d",aAOD->GetOrbitNumber())<<endl;
   cout<<Form("aAOD->GetPeriodNumber() = %d",aAOD->GetPeriodNumber())<<endl;
   if(!SpecifiedEvent(aAOD->GetRunNumber(),aAOD->GetBunchCrossNumber(),aAOD->GetOrbitNumber(),aAOD->GetPeriodNumber())){return;}
  } // else if(aAOD)
 } // if(fWaitForSpecifiedEvent)

} // void AliAnalysisTaskMultiparticleFemtoscopy::DoSomeDebugging(AliVEvent *ave)

//================================================================================================================

Int_t AliAnalysisTaskMultiparticleFemtoscopy::CurrentEventNumber()
{
 // Determine the current event number.

 Int_t currentEventNumber = -44;

 if(fAnalysisType->EqualTo("MC"))
 {
  currentEventNumber = (Int_t)fGetNumberOfTracksHist->GetEntries(); // TBI there is an issue clearly if fFillControlHistogramsEvent is not enabled. Yes, this is really shaky...
 }

 else if(fAnalysisType->EqualTo("ESD"))
 {
  // TBI
 }

 else if(fAnalysisType->Contains("AOD")) // TBI: Okay...?
 {
  currentEventNumber = (Int_t)fGetNumberOfTracksHist->GetEntries(); // TBI there is an issue clearly if fFillControlHistogramsEvent is not enabled. Yes, this is really shaky...
 }

 return currentEventNumber;

} // Int_t AliAnalysisTaskMultiparticleFemtoscopy::CurrentEventNumber()

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::OnlineMonitoring()
{
 // Per request, do some online monitoring.

 // a) Update regularly the output file;
 // b) Bail out after specified number of events.

 // a) Update regularly the output file:
 if(fUpdateOutputFile)
 {
  Int_t currentEventNo = this->CurrentEventNumber(); // TBI not supported yet for MC and ESD
  if(0 == currentEventNo % fUpdateFrequency)
  {
   cout<<Form("nEvts: %d",currentEventNo)<<endl;
   TFile *f = new TFile(fUpdateWhichOutputFile->Data(),"recreate");
   fHistList->Write(fHistList->GetName(),TObject::kSingleKey);
   f->Close();
  }
 } // if(fUpdateOutputFile)

 // b) Bail out after specified number of events:
 if(fMaxNumberOfEvents > 0)
 {
  Int_t currentEventNo = this->CurrentEventNumber(); // TBI not supported yet for MC and ESD
  if(fMaxNumberOfEvents == currentEventNo)
  {
   cout<<Form("\nPer request, bailing out after %d events in the file %s .\n",fMaxNumberOfEvents,fUpdateWhichOutputFile->Data())<<endl;
   TFile *f = new TFile(fUpdateWhichOutputFile->Data(),"recreate");
   fHistList->Write(fHistList->GetName(),TObject::kSingleKey);
   f->Close();
   exit(0);
  }
 } // if(fMaxNumberOfEvents > 0)

} // void AliAnalysisTaskMultiparticleFemtoscopy::OnlineMonitoring()

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::InsanityChecksUserCreateOutputObjects()
{
 // Insanity checks for UserCreateOutputObjects() method.

 // a) Insanity checks for global track cuts;
 // b) Disabled temporarily 4-p stuff. TBI a) fMixedEvents0[3] needs to be upgraded to fMixedEvents0[4], and b) EstimateBackground(AliVEvent *ave) needs to be re-implemented
 // c) Insanity checks for background.

 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::InsanityChecksUserCreateOutputObjects()";

 // a) Insanity checks for global track cuts:
 Int_t returnValueICFGTC = this->InsanityChecksForGlobalTrackCuts();
 if(0!=returnValueICFGTC)
 {
  cout<<Form("\n\nSomething is fundamentally wrong with global track cuts!!!!")<<endl;
  cout<<Form("InsanityChecksForGlobalTrackCuts() returns %d, check its implementation for an explanation of return values.\n\n",returnValueICFGTC)<<endl;
  Fatal(sMethodName.Data(),"if(0!=returnValueICFGTC)");
 } // if(0!=returnValueICFGTC)

 // b) Disabled temporarily 4-p stuff. TBI
 if(fFill4pCorrelationFunctions || fEstimate4pBackground)
 {
  cout<<Form("\n\nCalculation of 4-p stuff is not validated yet!!!! \n\n")<<endl; // TBI
  Fatal(sMethodName.Data(),"fFill4pCorrelationFunctions || fEstimate4pBackground");
 }

 // c) Insanity checks for background:
 if(!(fBackgroundOption == 0 || fBackgroundOption == 1))
 {
  cout<<Form("\n\nThis option for calculating background is not supported yet!!!! \n\n")<<endl; // TBI
  Fatal(sMethodName.Data(),"!(fBackgroundOption == 0 || fBackgroundOption == 1)");
 } // if(fBackgroundOption !=0 || fBackgroundOption != 1)
 if(fMaxBufferSize1>=50)
 {
  cout<<Form("\n\nfMaxBufferSize1 can be at maximum 50 !!!! \n\n")<<endl; // TBI
  Fatal(sMethodName.Data(),"fMaxBufferSize1>=50");
 }

} // void AliAnalysisTaskMultiparticleFemtoscopy::InsanityChecksUserCreateOutputObjects()

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::InsanityChecksUserExec()
{
 // Insanity checks for UserExec() method.

 // a) Insanity checks specific for all analyses types;
 // b) Insanity checks specific only for MC analyses;
 // c) Insanity checks specific only for ESD analyses;
 // d) Insanity checks specific only for AOD analyses.

 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::InsanityChecksUserExec()";

 // a) Insanity checks specific for all analyses types: // TBI is this really specific for all analyses types
 if(fBailOutAfterQA && !fFillQA){Fatal(sMethodName.Data(),"fBailOutAfterQA && !fFillQA");}
 if(fOnlineMonitoring && !fFillControlHistogramsEvent){Fatal(sMethodName.Data(),"fOnlineMonitoring && !fFillControlHistogramsEvent.\n\nAt the moment, fOnlineMonitoring can be used only if fFillControlHistogramsEvent is enabled.");}
 if(fFillControlHistogramsNonIdentifiedParticlesFTSF && !(fFilterBitFTSF==128)){Fatal(sMethodName.Data(),"FillControlHistogramsNonIdentifiedParticlesFTSF && !(fFilterBitFTSF==128)\n\nAt the moment, fFillControlHistogramsNonIdentifiedParticlesFTSF is validated only for TPC-only tracks.");}



 return; // TBI re-think and re-validate the rest



 // a2) TBI:
 for(Int_t pid=0;pid<5;pid++) // [0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pa=0;pa<2;pa++) // particle(+q)/antiparticle(-q)
  {
   for(Int_t ps=0;ps<2;ps++) // kPrimary/kFromDecayVtx
   {
    if(!fPIDCA[pid][pa][ps]) {Fatal(sMethodName.Data(),"!fPIDCA[pid][pa][ps]");}
    if(0 != fPIDCA[pid][pa][ps]->GetEntriesFast()){Fatal(sMethodName.Data(),"0 != fPIDCA[pid][pa][ps]->GetEntriesFast()"); }
   }
  }
 }
 for(Int_t pid=0;pid<1;pid++) // [0=Lambda,1=...]   // TBI is this really specific for all analyses types
 {
  if(!fPIDV0sCA[pid]){Fatal(sMethodName.Data(),"!fPIDV0sCA[pid]");}
  if(0 != fPIDV0sCA[pid]->GetEntriesFast()){Fatal(sMethodName.Data(),"0 != fPIDV0sCA[pid]->GetEntriesFast()");}
 }

 // b) Insanity checks specific only for MC analyses:
 if(fAnalysisType->EqualTo("MC"))
 {
  // TBI
 }

 // c) Insanity checks specific only for ESD analyses:
 if(fAnalysisType->EqualTo("ESD")) // TBI 'else if' ?
 {
  // TBI
 }

 // d) Insanity checks specific only for AOD analyses:
 if(fAnalysisType->Contains("AOD")) // TBI 'else if' ?
 {
  if(!fGlobalTracksAOD[0]){Fatal(sMethodName.Data(),"!fGlobalTracksAOD[0]");}
  if(0 != fGlobalTracksAOD[0]->GetSize()){Fatal(sMethodName.Data(),"0 != fGlobalTracksAOD[0]->GetSize()");}
  //if(fProcessBothKineAndReco && fEstimateBackground){Fatal(sMethodName.Data(),"fProcessBothKineAndReco && fEstimateBackground");} TBI re-think and re-enable
  // Remark: Note that the similar check is not needed for instance for fGlobalTracksAOD[1] and fGlobalTracksAOD[2],
  // which are used to estimate background, since some events can be stored in a buffer, until a suitable co-event is found,
  // to calculate background.
 } // if(fAnalysisType->Contains("AOD")

} // void AliAnalysisTaskMultiparticleFemtoscopy::InsanityChecksUserExec()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::MC(AliMCEvent *aMC)
{
 // Monte Carlo analysis is performed in this method.

 // a) Debugging;
 // b) Common event selection criteria;
 // c) Fill control histogram for global event observables;
 // d) Fill control histogram for particles;
 // e) Calculate correlation functions;
 // f) Calculate background;
 // g) V0s.

 // a) Debugging:
 if(fDoSomeDebugging){this->DoSomeDebugging(aMC);}

 // b) Common event selection criteria:
 if(!this->PassesCommonEventCuts(aMC)){return;}

 // c) Fill control histogram for global event observables:
 if(fFillControlHistogramsEvent){this->FillControlHistogramsEvent(aMC);}

 // d) Fill control histograms for particles:
 if(fFillControlHistogramsNonIdentifiedParticles || fFillControlHistogramsIdentifiedParticles){this->FillControlHistogramsParticle(aMC);}

 // e) Calculate correlation functions:
 if(fFillCorrelationFunctions){this->CalculateCorrelationFunctions(aMC);}

 // f) Calculate background:
 if(fEstimate2pBackground || fEstimate3pBackground || fEstimate4pBackground){this->EstimateBackground(aMC);}

 return;

 // g) V0s:
 // TBI

} // void AliAnalysisTaskMultiparticleFemtoscopy::MC(AliMCEvent *aMC)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::ESD(AliESDEvent *aESD)
{
 // ESD analysis is performed in this method.

 // ...
 if(aESD)
 {
  aESD = NULL; // TBI eliminating warnings temporarily
 }

} // void AliAnalysisTaskMultiparticleFemtoscopy::ESD(AliESDEvent *aESD)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::AOD(AliAODEvent *aAOD)
{
 // AOD analysis is performed in this method.

 // a) Debugging;
 // b) Common event selection criteria;
 // c) Filter out "normal global" tracks for default analysis and cut on their number;
 // d) Fill control histogram for global event observables;
 // e) Fill control histogram for particles;
 // f) Calculate correlation functions;
 // g) Calculate background;
 // h) Calculate test correlation functions;
 // i) Calculate test background;
 // j) 'Hybrid appriach';
 // *) V0s.

 // a) Debugging:
 if(fDoSomeDebugging){this->DoSomeDebugging(aAOD);}

 // b) Common event selection criteria:
 if(!this->PassesCommonEventCuts(aAOD)){return;}

 // c) Filter out "normal global" tracks for default analysis and cut on their number:
 // Remark 1: 'TPC-only' tracks and 'global constrained to vertex' come with negative ID, and are therefore not stored in fGlobalTracksAOD[*]
 // Remark 2: Default analysis means that it is not used for mixed events
 this->GlobalTracksAOD(aAOD,0); // [0] stands for default analysis
 if(0 == fGlobalTracksAOD[0]->GetSize()) return; // yes, go to next event
 if(fCutOnNumberOfGlobalTracks)
 {
  if(fGlobalTracksAOD[0]->GetSize() < fMinNumberOfGlobalTracks) return;
  if(fGlobalTracksAOD[0]->GetSize() > fMaxNumberOfGlobalTracks) return;
 }

 // d) Fill control histogram for global event observables:
 if(fFillControlHistogramsEvent){this->FillControlHistogramsEvent(aAOD);}

 // e) Fill control histograms for particles:
 if(fFillControlHistogramsNonIdentifiedParticles || fFillControlHistogramsNonIdentifiedParticlesFTSF || fFillControlHistogramsIdentifiedParticles){this->FillControlHistogramsParticle(aAOD);}

 // f) Calculate correlation functions:
 if(fFillCorrelationFunctions){this->CalculateCorrelationFunctions(aAOD);}

 // g) Calculate background:
 if(fEstimate2pBackground || fEstimate3pBackground || fEstimate4pBackground){this->EstimateBackground(aAOD);}

 // h) Calculate test correlation functions:
 Bool_t bFillCorrelationFunctionsTEST = kFALSE;
 for(Int_t t=0;t<10;t++)
 {
  bFillCorrelationFunctionsTEST = bFillCorrelationFunctionsTEST || fFillCorrelationFunctionsTEST[t];
 }
 if(bFillCorrelationFunctionsTEST){this->CalculateCorrelationFunctionsTEST(aAOD);}

 // i) Calculate test background:
 Bool_t bFillBackgroundTEST = kFALSE;
 for(Int_t t=0;t<10;t++)
 {
  bFillBackgroundTEST = bFillBackgroundTEST || fFillBackgroundTEST[t];
 }
 if(bFillBackgroundTEST){this->EstimateBackgroundTEST(aAOD);}

 // j) 'Hybrid appriach':
 if(fDoHybridApproach){this->DoHybridApproach(aAOD);}

 return;


 // *) V0s:
 V0s(aAOD); // TBI implement flag to enable/disable this call

} // void AliAnalysisTaskMultiparticleFemtoscopy::AOD(AliAODEvent *aAOD)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::V0s(AliVEvent *ave)
{
 // Analysis with V0s.

 // a) Determine Ali{MC,ESD,AOD}Event;

 // a) Determine Ali{MC,ESD,AOD}Event:
 AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 if(aMC)
 {
  // TBI
 }
 else if(aESD)
 {
  // TBI
 }
 else if(aAOD)
 {
  // AOD V0s:
  TClonesArray *caV0s = aAOD->GetV0s(); 
  Int_t index = 0;
  AliAODv0 *aAODv0 = NULL;
  Int_t nProngs = -44;
  while(caV0s->At(index))
  {
   aAODv0 = (AliAODv0*) caV0s->At(index++);
   if(!aAODv0) break;


   //AliAODv0 *temp = (AliAODv0*)fPIDV0sCA[0]->ConstructedAt(index-1);
   //temp = (AliAODv0*)aAODv0->Clone();
   
   //cout<<Form("indddex = %d, fPIDV0sCA[0]->GetEntries() = %d",index,fPIDV0sCA[0]->GetEntries())<<endl;

   // TBI...
   continue;

   nProngs = aAODv0->GetNProngs(); 
   fGetNProngsHist->Fill(nProngs);
   fMassK0ShortHist->Fill(aAODv0->MassK0Short());
   fMassLambdaHist->Fill(aAODv0->MassLambda());
   fMassAntiLambdaHist->Fill(aAODv0->MassAntiLambda());
   fOpenAngleV0Hist->Fill(aAODv0->OpenAngleV0());
   fRadiusV0Hist->Fill(aAODv0->RadiusV0());
   fDcaV0ToPrimVertexHist->Fill(aAODv0->DcaV0ToPrimVertex());
   fMomV0XHist->Fill(aAODv0->MomV0X());
   fMomV0YHist->Fill(aAODv0->MomV0Y());
   fMomV0ZHist->Fill(aAODv0->MomV0Z());
   fPtV0Hist->Fill(pow(aAODv0->Pt2V0(),0.5));
   fPseudoRapV0Hist->Fill(aAODv0->PseudoRapV0());
   fPAHist->Fill(aAODv0->Alpha(),aAODv0->PtArmV0());
   // Check sharing:

   fUniqueIDHistEBE->Fill(aAODv0->GetPosID());
   fUniqueIDHistEBE->Fill(aAODv0->GetNegID());

   /*
   //cout<<Form("V0PosID: %d , V0NegID: %d",aAODv0->GetPosID(),aAODv0->GetNegID())<<endl;
   Int_t trackPos = aAODv0->GetPosID()>=0 ? fGlobalTracksAOD[0]->GetValue(aAODv0->GetPosID()) : fGlobalTracksAOD[0]->GetValue(-(aAODv0->GetPosID()+1));
   Int_t trackNeg = aAODv0->GetNegID()>=0 ? fGlobalTracksAOD[0]->GetValue(aAODv0->GetNegID()) : fGlobalTracksAOD[0]->GetValue(-(aAODv0->GetNegID()+1));
   //cout<<Form("global : %d, global : %d", trackPos, trackNeg)<<endl;
   */

   if(-1 != fUniqueIDHistEBE->FindFirstBinAbove(1.44,1)) // TBI
   {
    cout<<Form("fUniqueIDHistEBE->FindFirstBinAbove(1.44) %d:",(Int_t)fUniqueIDHistEBE->FindFirstBinAbove(1.44,1))<<endl; // 
   }
  } // while(caV0s->At(index))


  /*
  Int_t index2 = 0;
  while(caV0s->At(index2))
  {
   cout<<((AliAODv0*) caV0s->At(index2))->Pt()<<endl;
   cout<<((AliAODv0*) fPIDV0sCA[0]->At(index2))->Pt()<<endl;
   index2++;
   cout<<endl;
  }
  */

   //   AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,globaltrack);
   //   AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,globaltrack);

 }

} // void AliAnalysisTaskMultiparticleFemtoscopy::V0s(AliVEvent *ave)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::ResetEBEObjects()
{
 // Reset all event-by-event objects.

 // a) Reset event-by-event objects specific for all analyses types;
 // b) Reset event-by-event objects specific only for MC analyses;
 // c) Reset event-by-event objects specific only for ESD analyses;
 // d) Reset event-by-event objects specific only for AOD analyses.

 // a) Reset event-by-event objects specific for all analyses types:
 if(fUniqueIDHistEBE) fUniqueIDHistEBE->Reset();
 // TBI add comment
 for(Int_t pid=0;pid<5;pid++) // [0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pa=0;pa<2;pa++) // particle(+q)/antiparticle(-q)
  {
   for(Int_t ps=0;ps<2;ps++) // kPrimary/kFromDecayVtx
   {
    if(fPIDCA[pid][pa][ps]) fPIDCA[pid][pa][ps]->Delete();
   }
  }
 }
 // TBI add comment
 for(Int_t pid=0;pid<1;pid++) // [0=Lambda,1=...] 
 {
  if(fPIDV0sCA[pid]) fPIDV0sCA[pid]->Delete();
 }

 // b) Reset event-by-event objects specific only for MC analyses:
 if(fAnalysisType->EqualTo("MC"))
 {
  // TBI
 }

 // c) Reset event-by-event objects specific only for ESD analyses:
 if(fAnalysisType->EqualTo("ESD")) // TBI 'else if' ?
 {
  // TBI
 }

 // d) Reset event-by-event objects specific only for AOD analyses:
 if(fAnalysisType->Contains("AOD")) // TBI 'else if' ?
 {
  if(fGlobalTracksAOD[0]) fGlobalTracksAOD[0]->Delete();
  // Remark: Note that fGlobalTracksAOD[1] and fGlobalTracksAOD[2] used for event mixing do not have
  // to be reset e-b-e. Where then I reset them? TBI
 }

} // void AliAnalysisTaskMultiparticleFemtoscopy::ResetEBEObjects()

//=======================================================================================================================

Bool_t AliAnalysisTaskMultiparticleFemtoscopy::Pion(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)
{
 // Is this primary/secondary (anti-)pion?

 // a) Insanity checks;
 // b) Trivial checks; // TBI
 // c) Track quality cuts;
 // d) PID;
 // e) PID MC (if available).

 // a) Insanity checks
 TString sMethodName = "Bool_t AliAnalysisTaskMultiparticleFemtoscopy::Pion(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)";
 if(!(1 == charge || -1 == charge)){Fatal(sMethodName.Data(),"!(1 == charge || -1 == charge)");}
 if(!gtrack){Fatal(sMethodName.Data(),"!gtrack");}

 // b) Trivial checks:
 if(charge != (Int_t)gtrack->Charge()) return kFALSE;
 if(bPrimary && gtrack->GetType() != AliAODTrack::kPrimary)
 {
  return kFALSE;
 }
 else if(!bPrimary && gtrack->GetType() != AliAODTrack::kFromDecayVtx)
 {
  return kFALSE;
 }

 // c) Track quality cuts:
/*
 if(bPrimary)
 {
  // TBI  
 }
*/

 // d) PID:
 // For pT < 0.75 use only TPC
 if(gtrack->GetTPCmomentum() < 0.75) // TBI hardwired 0.75
 {
  AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,gtrack);
  if(!statusTPC) return kFALSE;
  // sigmas:
  Double_t dSigmaTPCPion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(gtrack,(AliPID::kPion)));
  Double_t dSigmaTPCKaon = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(gtrack,(AliPID::kKaon)));
  Double_t dSigmaTPCProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(gtrack,(AliPID::kProton)));
  // inclusive cuts and exclusive cuts:
  if(!(dSigmaTPCPion < fInclusiveSigmaCuts[2] && dSigmaTPCProton > fExclusiveSigmaCuts[2][4] && dSigmaTPCKaon > fExclusiveSigmaCuts[2][3])) return kFALSE;
 } // if(gtrack->GetTPCmomentum() <= 0.75) // TBI hardwired 0.75
 else if(gtrack->GetTPCmomentum() >= 0.75 )
 {
  // TBI use combined TPC and TOf
  return kFALSE; // TBI remove eventually
 }

 // e) PID MC (if available):
 if(fProcessBothKineAndReco)
 {
  if(!fMC){Fatal(sMethodName.Data(),"!fMC");}
  AliAODMCParticle *mcParticle = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(TMath::Abs(gtrack->GetLabel())));
  if(!mcParticle){cout<<"WARNING: mcParticle is NULL in Pion(...)! gtrack = "<<gtrack<<endl; return kFALSE; } // TBI well, it remains undetermined, investigate further why in some very rare cases I cannot get this pointer. if I get too many warnings, that purity estimation will be affected
  if(charge < 0 && mcParticle->GetPdgCode() == -211) return kTRUE;
  else if(charge > 0 && mcParticle->GetPdgCode() == 211) return kTRUE;
  else return kFALSE;
 } // if(fProcessBothKineAndReco)

 return kTRUE;

} // Bool_t AliAnalysisTaskMultiparticleFemtoscopy::Pion(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)


//=======================================================================================================================

Bool_t AliAnalysisTaskMultiparticleFemtoscopy::Kaon(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)
{
 // Is this primary/secondary (anti-)kaon?

 // a) Insanity checks;
 // b) Trivial checks;
 // c) Track quality cuts;
 // d) PID;
 // e) PID MC (if available).

 // a) Insanity checks:
 TString sMethodName = "Bool_t AliAnalysisTaskMultiparticleFemtoscopy::Kaon(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)";
 if(!(1 == charge || -1 == charge)){Fatal(sMethodName.Data(),"!(1 == charge || -1 == charge)");}
 if(!gtrack){Fatal(sMethodName.Data(),"!gtrack");}

 // b) Trivial checks:
 if(charge != (Int_t)gtrack->Charge()) return kFALSE;
 if(bPrimary && gtrack->GetType() != AliAODTrack::kPrimary)
 {
  return kFALSE;
 }
 else if(!bPrimary && gtrack->GetType() != AliAODTrack::kFromDecayVtx)
 {
  return kFALSE;
 }

 // c) Track quality cuts:
 /*
 if(bPrimary)
 {
  // TBI  
 }
 */

 // d) PID:
 // For pT < 0.75 use only TPC
 if(gtrack->GetTPCmomentum() < 0.75) // TBI hardwired 0.75
 {
  AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,gtrack);
  if(!statusTPC) return kFALSE;
  // sigmas:
  Double_t dSigmaTPCKaon = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(gtrack,(AliPID::kKaon)));
  Double_t dSigmaTPCProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(gtrack,(AliPID::kProton)));
  Double_t dSigmaTPCPion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(gtrack,(AliPID::kPion)));
  // inclusive and exclusive cuts:
  if(!(dSigmaTPCKaon < fInclusiveSigmaCuts[3] && dSigmaTPCProton > fExclusiveSigmaCuts[3][4] && dSigmaTPCPion > fExclusiveSigmaCuts[3][2])) return kFALSE;
 } // if(gtrack->GetTPCmomentum() <= 0.75) // TBI hardwired 0.75
 else if(gtrack->GetTPCmomentum() >= 0.75 )
 {
  // TBI use combined TPC and TOf
  return kFALSE; // TBI remove eventually
 }

 // e) PID MC (if available):
 if(fProcessBothKineAndReco)
 {
  if(!fMC){Fatal(sMethodName.Data(),"!fMC");}
  AliAODMCParticle *mcParticle = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(TMath::Abs(gtrack->GetLabel())));
  if(!mcParticle){cout<<"WARNING: mcParticle is NULL in Kaon(...)! gtrack = "<<gtrack<<endl; return kFALSE; } // TBI well, it remains undetermined, investigate further why in some very rare cases I cannot get this pointer. if I get too many warnings, that purity estimation will be affected
  if(charge < 0 && mcParticle->GetPdgCode() == -321) return kTRUE;
  else if(charge > 0 && mcParticle->GetPdgCode() == 321) return kTRUE;
  else return kFALSE;
 } // if(fProcessBothKineAndReco)

 return kTRUE;

} // Bool_t AliAnalysisTaskMultiparticleFemtoscopy::Kaon(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)

//=======================================================================================================================

Bool_t AliAnalysisTaskMultiparticleFemtoscopy::Proton(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)
{
 // Is this primary/secondary (anti-)proton?

 // a) Insanity checks;
 // b) Trivial checks;
 // c) Track quality cuts;
 // d) PID;
 // e) PID MC (if available).

 // a) Insanity checks:
 TString sMethodName = "Bool_t AliAnalysisTaskMultiparticleFemtoscopy::Proton(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)";
 if(!(1 == charge || -1 == charge)){Fatal(sMethodName.Data(),"!(1 == charge || -1 == charge)");}
 if(!gtrack){Fatal(sMethodName.Data(),"!gtrack");}

 // b) Trivial checks:
 if(charge != (Int_t)gtrack->Charge()) return kFALSE;
 if(bPrimary && gtrack->GetType() != AliAODTrack::kPrimary)
 {
  return kFALSE;
 }
 else if(!bPrimary && gtrack->GetType() != AliAODTrack::kFromDecayVtx)
 {
  return kFALSE;
 }
 
 // c) Track quality cuts:
/*
 if(bPrimary)
 {
  // TBI  
 }
*/

 // d) PID:
 // For pT < 0.75 use only TPC
 if(gtrack->GetTPCmomentum() < 0.75) // TBI hardwired 0.75
 {
  AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,gtrack);
  if(!statusTPC) return kFALSE;
  // sigmas:
  Double_t dSigmaTPCProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(gtrack,(AliPID::kProton)));
  Double_t dSigmaTPCPion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(gtrack,(AliPID::kPion)));
  Double_t dSigmaTPCKaon = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(gtrack,(AliPID::kKaon)));
  if(!(dSigmaTPCProton < fInclusiveSigmaCuts[4] && dSigmaTPCPion > fExclusiveSigmaCuts[4][2] && dSigmaTPCKaon > fExclusiveSigmaCuts[4][3])) return kFALSE;
 } // if(gtrack->GetTPCmomentum() <= 0.75) // TBI hardwired 0.75
 else if(gtrack->GetTPCmomentum() >= 0.75 )
 {
  // TBI use combined TPC and TOf
  return kFALSE; // TBI remove eventually
 }

 // e) PID MC (if available):
 if(fProcessBothKineAndReco)
 {
  if(!fMC){Fatal(sMethodName.Data(),"!fMC");}
  AliAODMCParticle *mcParticle = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(TMath::Abs(gtrack->GetLabel())));
  if(!mcParticle){cout<<"WARNING: mcParticle is NULL in Proton(...)! gtrack = "<<gtrack<<endl; return kFALSE; } // TBI well, it remains undetermined, investigate further why in some very rare cases I cannot get this pointer. if I get too many warnings, that purity estimation will be affected
  if(charge < 0 && mcParticle->GetPdgCode() == -2212) return kTRUE;
  else if(charge > 0 && mcParticle->GetPdgCode() == 2212) return kTRUE;
  else return kFALSE;
 } // if(fProcessBothKineAndReco)

 return kTRUE;

} // Bool_t AliAnalysisTaskMultiparticleFemtoscopy::Proton(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArrays()
{
 // Initialize arrays for all objects not classified yet.

 for(Int_t index=0;index<10;index++) // [0] is used in the default analysis, [1] and [2] for event mixing, etc.
 {
  fGlobalTracksAOD[index] = NULL;
 }

} // void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArrays()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForControlHistograms()
{
 // Initialize all arrays for control histograms.

 for(Int_t xyz=0;xyz<3;xyz++)
 {
  fVertexXYZ[xyz] = NULL;               
 }

 for(Int_t pid=0;pid<5;pid++) // [0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pa=0;pa<2;pa++) // particle(+q)/antiparticle(-q)
  {
   for(Int_t ps=0;ps<2;ps++) // kPrimary/kFromDecayVtx
   {
    fMassPIDHist[pid][pa][ps] = NULL;
    fPtPIDHist[pid][pa][ps] = NULL;
    for(Int_t xyz=0;xyz<3;xyz++) // kPrimary/kFromDecayVtx
    {
     fPPIDHist[pid][pa][ps][xyz] = NULL;
    }
    fEtaPIDHist[pid][pa][ps] = NULL;
    fPhiPIDHist[pid][pa][ps] = NULL;
   }
  }
 } // for(Int_t pid=0;pid<4;pid++) // [0=e,1=mu,2=pi,3=K,4=p]

 // Inclusive sigma cuts:
 for(Int_t pf=0;pf<5;pf++) // PID function [0=Electron(...),1=Muon(...),2=Pion(...),3=Kaon(...),4=Proton(...)]
 {
  fInclusiveSigmaCuts[pf] = 0.;
 }

 // Exclusive sigma cuts:
 for(Int_t pf=0;pf<5;pf++) // PID function [0=Electron(...),1=Muon(...),2=Pion(...),3=Kaon(...),4=Proton(...)]
 {
  for(Int_t pid=0;pid<5;pid++) // [0=e,1=mu,2=pi,3=K,4=p]
  {
   fExclusiveSigmaCuts[pf][pid] = 0.;
  }
 }

 // Normalization:
 fNormalizationInterval[0] = -0.44;
 fNormalizationInterval[1] = -0.44;

} // void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForControlHistograms()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForEBEObjects()
{
 // Initialize all arrays for e-b-e objects.

 for(Int_t pid=0;pid<5;pid++) // [0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pa=0;pa<2;pa++) // particle(+q)/antiparticle(-q)
  {
   for(Int_t ps=0;ps<2;ps++) // kPrimary/kFromDecayVtx
   {
    fPIDCA[pid][pa][ps] = NULL;
   }
  }
 } // for(Int_t pid=0;pid<4;pid++) // [0=e,1=mu,2=pi,3=K,4=p]


 for(Int_t pid=0;pid<1;pid++)
 {
  fPIDV0sCA[pid] = NULL;               
 }

} // void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForEBEObjects()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForCorrelationFunctions()
{
 // Initialize all arrays for correlation functions.

 // a) Initialize sublists;
 // b) Initialize 2p correlation functions;
 // c) Initialize 3p correlation functions;
 // d) Initialize 4p correlation functions.

 // a) Initialize sublists:
 for(Int_t cfs=0;cfs<3;cfs++)
 {
  fCorrelationFunctionsSublist[cfs] = NULL;
 }

 // b) Initialize 2p correlation functions:
 for(Int_t pid1=0;pid1<10;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pid2=0;pid2<10;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   fCorrelationFunctions[pid1][pid2] = NULL;
  }
 }

 // c) Initialize 3p correlation functions:
 for(Int_t pid1=0;pid1<10;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pid2=0;pid2<10;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   for(Int_t pid3=0;pid3<10;pid3++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
   {
    f3pCorrelationFunctions[pid1][pid2][pid3] = NULL;
   }
  }
 }

 // c) Initialize 4p correlation functions:
 for(Int_t pid1=0;pid1<10;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pid2=0;pid2<10;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   for(Int_t pid3=0;pid3<10;pid3++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
   {
    for(Int_t pid4=0;pid4<10;pid4++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
    {
     f4pCorrelationFunctions[pid1][pid2][pid3][pid4] = NULL;
    }
   }
  }
 }

} // void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForCorrelationFunctions()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForBackground()
{
 // Initialize all arrays for background.

 for(Int_t bs=0;bs<3;bs++)
 {
  fBackgroundSublist[bs] = NULL;        // lists to hold all background correlations, for 2p [0], 3p [1], 4p [2], etc., separately
 }

 for(Int_t pid1=0;pid1<10;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pid2=0;pid2<10;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   f2pBackground[pid1][pid2] = NULL;
  }
 }

 for(Int_t pid1=0;pid1<10;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pid2=0;pid2<10;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   for(Int_t pid3=0;pid3<10;pid3++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
   {
    f3pBackground[pid1][pid2][pid3] = NULL;
   }
  }
 }

 for(Int_t me=0;me<3;me++) // [0] is buffer for 1st event; [1] for 2nd, etc.
 {
  fMixedEvents0[me] = NULL;
 }

 for(Int_t vzr=0;vzr<10;vzr++) // vertez-z range
 {
  for(Int_t me=0;me<fMaxBufferSize1;me++) // mixed events
  {
   fMixedEvents1[vzr][me] = NULL;
   fGlobalTracksAOD1[vzr][me] = NULL;
  }
 }

} // void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForBackground()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForBuffers()
{
 // Initialize all arrays for buffers.

 for(Int_t am=0;am<2;am++) // analysis method [0=AOD||ESD,1=MC]
 {
  for(Int_t e=0;e<10;e++) // number of events to buffer, max = apparently 10
  {
   for(Int_t p=0;p<10000;p++) // charged particles, max = apparently 10000
   {
    fChargedParticlesCA[am][e][p] = NULL;
    // Example: fChargedParticlesCA[0][4][6] is 6th charged particles in 5th buffered event, in AOD||ESD analysis
   } // for(Int_t p=0;p<10000;p++) // charged particles, max = apparently 10000
  } // for(Int_t e=0;e<10;e++) // number of events to buffer, max = apparently 10
 } // for(Int_t am=0;am<2;am++) // analysis method [0=AOD||ESD,1=MC]

 for(Int_t e=0;e<10;e++) // number of events to buffer, max = apparently 10
 {
  fChargedParticlesEM[e] = NULL;
 } // for(Int_t e=0;e<10;e++) // number of events to buffer, max = apparently 10

} // void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForBuffers()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForQA()
{
 // Initialize all arrays for QA.

 for(Int_t ba=0;ba<2;ba++) // 0="before rain",1="after rain"
 {
  for(Int_t di=0;di<10;di++) // number of distinct distributions, maximum is 10 (at the moment, as it seems...)
  {
   for(Int_t ci=0;ci<10;ci++) // number of distinct track cuts, maximum is 10 (at the moment, as it seems...)
   {
    fQAParticleHist[ba][di][ci] = NULL;
    // Example: fQAParticleHist[0][2][4] is distribution of 3rd observable BEFORE the 5th cut was applied
    //          fQAParticleHist[1][2][4] is distribution of 3rd observable AFTER the 5th cut was applied
   } // for(Int_t ci=0;ci<10;ci++) // number of distinct track cuts, maximum is 10 (at the moment)
  } // for(Int_t di=0;di<10;di++) // number of distinct distributions, maximum is 10 (at the moment)
 } // for(Int_t ba=0;ba<2;ba++) // 0="before rain",1="after rain"

} // void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForQA()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForGlobalTrackCuts()
{
 // Initialize all arrays for common global tracks cuts. In essence, these are the default cuts for "normal" global tracks in AOD.

 // The default values, can be overwritten with the setters SetPtRange(...), etc.
 fPtRange[0] = 0.2; // pt min
 fPtRange[1] = 5.0; // pt max
 fEtaRange[0] = -0.8; // eta min
 fEtaRange[1] = 0.8; // eta max
 fPhiRange[0] = 0.; // phi min
 fPhiRange[1] = TMath::TwoPi(); // phi max

} // void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForGlobalTrackCuts()

//=======================================================================================================================

Int_t AliAnalysisTaskMultiparticleFemtoscopy::InsanityChecksForGlobalTrackCuts()
{
 // Isanity checks for global track cuts.

 // Return values:
 // 1) pt range;
 // 2) eta range;
 // 3) phi range;

 if(fPtRange[0]<0.||fPtRange[1]<0.){cout<<Form("fPtRange[0] = %f, fPtRange[1] = %f",fPtRange[0],fPtRange[1])<<endl; return 1;} // pt is positive
 if(fPtRange[0]>fPtRange[1]){cout<<Form("fPtRange[0] = %f, fPtRange[1] = %f",fPtRange[0],fPtRange[1])<<endl; return 1;} // ptmin > ptmax

 if(fEtaRange[0]>fEtaRange[1]){cout<<Form("fEtaRange[0] = %f, fEtaRange[1] = %f",fEtaRange[0],fEtaRange[1])<<endl; return 2;} // etamin > etamax

 if(fPhiRange[0]<0.||fPhiRange[1]<0.){cout<<Form("fPhiRange[0] = %f, fPhiRange[1] = %f",fPhiRange[0],fPhiRange[1])<<endl; return 3;} // phi is positive
 //if(fPhiRange[0]>TMath::TwoPi()||fPhiRange[1]>TMath::TwoPi()){return 3;} // phi shall be smaller than 2\pi TBI
 if(fPhiRange[0]>fPhiRange[1]){cout<<Form("fPhiRange[0] = %f, fPhiRange[1] = %f",fPhiRange[0],fPhiRange[1])<<endl; return 3;} // phimin > phimax

 return 0; // Life is like a wheel...

} // Int_t AliAnalysisTaskMultiparticleFemtoscopy::InsanityChecksForGlobalTrackCuts()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForCorrelationFunctionsTEST()
{
 // Initialize all arrays for test correlation functions.

 const Int_t nTestsMax = 10; // see what is hardwired in .h file

 for(Int_t t=0;t<nTestsMax;t++)
 {
  fCorrelationFunctionsTESTSublist[t] = NULL;
  fFillCorrelationFunctionsTEST[t] = kFALSE;
 }

 for(Int_t t=0;t<nTestsMax;t++)
 {
  for(Int_t q=0;q<2;q++) // Q2 and Q3
  {
   for(Int_t tt=0;tt<7;tt++) // 7 distinct terms in the expression of 3p cumulant. TBI for 2p, this is an overkill, but for 4p this will need to be enlarged
   {
    for(Int_t d=0;d<10;d++) // differential options, see what is hardwired in .h file
    {
     fCorrelationFunctionsTEST[t][q][tt][d] = NULL;
    }
   } // for(Int_t t=0;t<nTests;t++)
  } // for(Int_t q=0;q<2;q++) // Q2 and Q3
 } // for(Int_t t=0;t<nTests;t++)

 for(Int_t t=0;t<nTestsMax;t++)
 {
  for(Int_t q=0;q<2;q++) // Q2 and Q3
  {
   for(Int_t tt=0;tt<4;tt++) // 4 distinct cumulants, three 2p and one 3p. TBI for 2p, this is an overkill, but for 4p this will need to be enlarged
   {
    for(Int_t d=0;d<10;d++) // differential options, see what is hardwired in .h file
    {
     fSignalCumulantsTEST[t][q][tt][d] = NULL;
    }
   } // for(Int_t t=0;t<nTests;t++)
  } // for(Int_t q=0;q<2;q++) // Q2 and Q3
 } // for(Int_t t=0;t<nTests;t++)

 for(Int_t t=0;t<2;t++)
 {
  fSignalYieldTEST[t] = NULL;
  fEab_TEST6[t] = NULL;
 }

} // void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForCorrelationFunctionsTEST()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForBackgroundTEST()
{
 // Initialize all arrays for test background.

 const Int_t nTestsMax = 10; // see what is hardwired in .h file

 for(Int_t t=0;t<nTestsMax;t++)
 {
  fBackgroundTESTSublist[t] = NULL;
  fFillBackgroundTEST[t] = kFALSE;
 }

 for(Int_t t=0;t<nTestsMax;t++)
 {
  for(Int_t q=0;q<2;q++) // Q2 and Q3
  {
   for(Int_t tt=0;tt<7;tt++) // 7 distinct terms in the expression of 3p cumulant. TBI for 2p, this is an overkill, but for 4p this will need to be enlarged
   {
    for(Int_t d=0;d<10;d++) // differential options, see what is hardwired in .h file
    {
     fBackgroundTEST[t][q][tt][d] = NULL;
    }
   } // for(Int_t tt=0;tt<7;tt++)
  } // for(Int_t q=0;q<2;q++) // Q2 and Q3
 } // for(Int_t t=0;t<nTestsMax;t++)

 for(Int_t me=0;me<3;me++) // [0] is buffer for 1st event; [1] for 2nd, etc.
 {
  fMixedEventsTEST[me] = NULL;
  fGlobalTracksAODTEST[me] = NULL;
 }

 for(Int_t t=0;t<2;t++)
 {
  fBackgroundYieldTEST[t] = NULL;
 }

} // void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForBackgroundTEST()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForHybridApproach()
{
 // Initialize all arrays for 'hybrid approach'.

 for(Int_t t=0;t<5;t++) // five distinct terms in the numerator of Eq. (18)
 {
  fDistributionsVsQ3[t] = NULL;
 }

 for(Int_t me=0;me<3;me++) // [0] is buffer for 1st event; [1] for 2nd, etc.
 {
  fMixedEventsHA[me] = NULL;
  fGlobalTracksAODHA[me] = NULL;
 }

} // void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForHybridApproach()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForCorrelationFunctionsTEST()
{
 // Book everything for test correlation functions.

 // a) Book the profile holding all the flags for test correlation functions;
 // b) Book correlation functions for all tests;
 // c) Store temporarely also the yield.

 // a) Book the profile holding all the flags for test correlation functions:
 const Int_t nTests = 10;
 fCorrelationFunctionsTESTFlagsPro = new TProfile("fCorrelationFunctionsTESTFlagsPro","Flags and settings for test correlation functions",nTests,0,nTests);
 fCorrelationFunctionsTESTFlagsPro->SetTickLength(-0.01,"Y");
 fCorrelationFunctionsTESTFlagsPro->SetMarkerStyle(25);
 fCorrelationFunctionsTESTFlagsPro->SetLabelSize(0.04);
 fCorrelationFunctionsTESTFlagsPro->SetLabelOffset(0.02,"Y");
 fCorrelationFunctionsTESTFlagsPro->SetStats(kFALSE);
 fCorrelationFunctionsTESTFlagsPro->SetFillColor(kGray);
 fCorrelationFunctionsTESTFlagsPro->SetLineColor(kBlack);
 for(Int_t t=0;t<nTests;t++)
 {
  fCorrelationFunctionsTESTFlagsPro->GetXaxis()->SetBinLabel(t+1,Form("fFillCorrelationFunctionsTEST[%d]",t)); fCorrelationFunctionsTESTFlagsPro->Fill(t+0.5,fFillCorrelationFunctionsTEST[t]);
 }
 fCorrelationFunctionsTESTList->Add(fCorrelationFunctionsTESTFlagsPro);

 // b) Book correlation functions for all tests:
 const Int_t nCumulants = 4; // TBI this applies only to vs. Q3
 const Int_t n3pCumulantTerms = 7;
 const Int_t nQs = 2; // Q2 and Q3 at the moment, if you change it here, change it also in exprs. like for(Int_t q=0;q<2;q++)
 const Int_t n2pCumulantTerms = 3;
 TString s2pCumulantTerms[n2pCumulantTerms] = {"#LTX_{1}#GT","#LTX_{2}#GT","#LTX_{1}X_{2}#GT"};
 TString s3pCumulantTerms[n3pCumulantTerms] = {"#LTX_{1}#GT","#LTX_{2}#GT","#LTX_{3}#GT","#LTX_{1}X_{2}#GT","#LTX_{1}X_{3}#GT","#LTX_{2}X_{3}#GT","#LTX_{1}X_{2}X_{3}#GT"};
 TString sCumulants[nCumulants] = {"#LTX_{1}X_{2}#GT_{c}","#LTX_{1}X_{3}#GT_{c}","#LTX_{2}X_{3}#GT_{c}","#LTX_{1}X_{2}X_{3}#GT_{c}"}; // TBI this applies only to vs. Q3
 TString sXYZ[3] = {"x","y","z"};
 TString sQs[nQs] = {"Q_{2}","Q_{3}"};
 for(Int_t t=0;t<nTests;t++)
 {
  if(fFillCorrelationFunctionsTEST[t])
  {

   // Correlations vs. Q2:
   for(Int_t ct=0;ct<n2pCumulantTerms;ct++)
   {
    for(Int_t xyz=0;xyz<3;xyz++)
    {
     fCorrelationFunctionsTEST[t][0][ct][xyz] = new TProfile(Form("fCorrelationFunctionsTEST[%d][0][%d][%d]",t,ct,xyz),Form("TEST: %d",t),fnQ2bins,fnQ2min,fnQ2max);
     fCorrelationFunctionsTEST[t][0][ct][xyz]->SetStats(kFALSE);
     fCorrelationFunctionsTEST[t][0][ct][xyz]->SetMarkerStyle(kFullSquare);
     fCorrelationFunctionsTEST[t][0][ct][xyz]->SetMarkerColor(kBlue);
     fCorrelationFunctionsTEST[t][0][ct][xyz]->SetLineColor(kBlue);
     fCorrelationFunctionsTEST[t][0][ct][xyz]->GetXaxis()->SetTitle(sQs[0].Data());
     fCorrelationFunctionsTEST[t][0][ct][xyz]->GetYaxis()->SetTitle(Form("%s_{%s}",s2pCumulantTerms[ct].Data(),sXYZ[xyz].Data()));
     fCorrelationFunctionsTESTSublist[t]->Add(fCorrelationFunctionsTEST[t][0][ct][xyz]);
    } // for(Int_t xyz=0;xyz<3;xyz++)
   } // for(Int_t ct=0;ct<n2pCumulantTerms;ct++)

   // Cumulants vs. Q2:
   for(Int_t xyz=0;xyz<3;xyz++)
   {
    fSignalCumulantsTEST[t][0][0][xyz] = new TProfile(Form("fSignalCumulantsTEST[%d][0][%d][%d]",t,0,xyz),Form("TEST: %d",t),fnQ2bins,fnQ2min,fnQ2max);
    fSignalCumulantsTEST[t][0][0][xyz]->SetStats(kFALSE);
    fSignalCumulantsTEST[t][0][0][xyz]->SetMarkerStyle(kFullSquare);
    fSignalCumulantsTEST[t][0][0][xyz]->SetMarkerColor(kGreen+2);
    fSignalCumulantsTEST[t][0][0][xyz]->SetLineColor(kGreen+2);
    fSignalCumulantsTEST[t][0][0][xyz]->GetXaxis()->SetTitle(sQs[0].Data());
    fSignalCumulantsTEST[t][0][0][xyz]->GetYaxis()->SetTitle(Form("#LTX_{1}X_{2}#GT_{c,%s}",sXYZ[xyz].Data()));
    fCorrelationFunctionsTESTSublist[t]->Add(fSignalCumulantsTEST[t][0][0][xyz]);
   } // for(Int_t xyz=0;xyz<3;xyz++)

   // Correlations vs. Q3: TBI this clearly can be unified with the above.
   for(Int_t ct=0;ct<n3pCumulantTerms;ct++)
   {
    for(Int_t xyz=0;xyz<3;xyz++)
    {
     fCorrelationFunctionsTEST[t][1][ct][xyz] = new TProfile(Form("fCorrelationFunctionsTEST[%d][1][%d][%d]",t,ct,xyz),Form("TEST: %d",t),fnQ3bins,fnQ3min,fnQ3max);
     fCorrelationFunctionsTEST[t][1][ct][xyz]->SetStats(kFALSE);
     fCorrelationFunctionsTEST[t][1][ct][xyz]->SetMarkerStyle(kFullSquare);
     fCorrelationFunctionsTEST[t][1][ct][xyz]->SetMarkerColor(kBlue);
     fCorrelationFunctionsTEST[t][1][ct][xyz]->SetLineColor(kBlue);
     fCorrelationFunctionsTEST[t][1][ct][xyz]->GetXaxis()->SetTitle(sQs[1].Data());
     fCorrelationFunctionsTEST[t][1][ct][xyz]->GetYaxis()->SetTitle(Form("%s_{%s}",s3pCumulantTerms[ct].Data(),sXYZ[xyz].Data()));
     fCorrelationFunctionsTESTSublist[t]->Add(fCorrelationFunctionsTEST[t][1][ct][xyz]);
    } // for(Int_t xyz=0;xyz<3;xyz++)
   } // for(Int_t ct=0;ct<n3pCumulantTerms;ct++)

   // Cumulants vs. Q3:
   for(Int_t ct=0;ct<nCumulants;ct++) // TBI this applies only to vs. Q3
   {
    for(Int_t xyz=0;xyz<3;xyz++)
    {
     fSignalCumulantsTEST[t][1][ct][xyz] = new TProfile(Form("fSignalCumulantsTEST[%d][1][%d][%d]",t,ct,xyz),Form("TEST: %d",t),fnQ3bins,fnQ3min,fnQ3max);
     fSignalCumulantsTEST[t][1][ct][xyz]->SetStats(kFALSE);
     fSignalCumulantsTEST[t][1][ct][xyz]->SetMarkerStyle(kFullSquare);
     fSignalCumulantsTEST[t][1][ct][xyz]->SetMarkerColor(kGreen+2);
     fSignalCumulantsTEST[t][1][ct][xyz]->SetLineColor(kGreen+2);
     fSignalCumulantsTEST[t][1][ct][xyz]->GetXaxis()->SetTitle(sQs[1].Data());
     fSignalCumulantsTEST[t][1][ct][xyz]->GetYaxis()->SetTitle(Form("%s_{%s}",sCumulants[ct].Data(),sXYZ[xyz].Data()));
     fCorrelationFunctionsTESTSublist[t]->Add(fSignalCumulantsTEST[t][1][ct][xyz]);
    } // for(Int_t xyz=0;xyz<3;xyz++)
   } // for(Int_t ct=0;ct<nCumulants;ct++)

  } // if(fFillCorrelationFunctionsTEST[t])

 } // for(Int_t t=0;t<nTests;t++)

 // c) Store temporarely also the yield:
 // vs. Q2:
 fSignalYieldTEST[0] = new TH1F("fSignalYieldTEST[0]","dN(X_{1}X_{2})/Q_{2} for signal (TEST_0)",fnQ2bins,fnQ2min,fnQ2max);
 fSignalYieldTEST[0]->SetStats(kFALSE);
 fSignalYieldTEST[0]->SetLineColor(kBlue);
 fSignalYieldTEST[0]->SetFillColor(kBlue-10);
 fSignalYieldTEST[0]->GetXaxis()->SetTitle("Q_{2}");
 fSignalYieldTEST[0]->GetYaxis()->SetTitle("counts");
 fCorrelationFunctionsTESTList->Add(fSignalYieldTEST[0]);
 // vs. Q3:
 fSignalYieldTEST[1] = new TH1F("fSignalYieldTEST[1]","dN(X_{1}X_{2}X_{3})/Q_{3} for signal (TEST_1)",fnQ3bins,fnQ3min,fnQ3max);
 fSignalYieldTEST[1]->SetStats(kFALSE);
 fSignalYieldTEST[1]->SetLineColor(kBlue);
 fSignalYieldTEST[1]->SetFillColor(kBlue-10);
 fSignalYieldTEST[1]->GetXaxis()->SetTitle("Q_{3}");
 fSignalYieldTEST[1]->GetYaxis()->SetTitle("counts");
 fCorrelationFunctionsTESTList->Add(fSignalYieldTEST[1]);

 // TBI temporary here for "Test 6" (checking E_a^b as L.I. observable):
 // vs. Q2:
 fEab_TEST6[0] = new TH1F("fEab_TEST6[0]","",10*fnQ2bins,fnQ2min,fnQ2max);
 fEab_TEST6[0]->SetStats(kFALSE);
 fEab_TEST6[0]->SetLineColor(kBlue);
 fEab_TEST6[0]->SetFillColor(kBlue-10);
 fEab_TEST6[0]->GetXaxis()->SetTitle("Q_{2}");
 fEab_TEST6[0]->GetYaxis()->SetTitle("dN/dE_{a}^{b}");
 fCorrelationFunctionsTESTList->Add(fEab_TEST6[0]);

 fEab_TEST6[1] = new TH1F("fEab_TEST6[1]","",10*fnQ2bins,fnQ2min,fnQ2max);
 fEab_TEST6[1]->SetStats(kFALSE);
 fEab_TEST6[1]->SetLineColor(kRed);
 fEab_TEST6[1]->SetFillColor(kRed-10);
 fEab_TEST6[1]->GetXaxis()->SetTitle("Q_{2}");
 fEab_TEST6[1]->GetYaxis()->SetTitle("dN/dE_{a}^{b}");
 fCorrelationFunctionsTESTList->Add(fEab_TEST6[1]);

} // void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForCorrelationFunctionsTEST()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForBackgroundTEST()
{
 // Book everything for test background.

 // a) Book the profile holding all the flags for test background;
 // b) Book background objects for all tests;
 // c) Book fMixedEventsTEST[3] and fGlobalTracksAODTEST[3];
 // d) Store temporary the yield as well.

 // a) Book the profile holding all the flags for test background:
 const Int_t nTests = 10;
 fBackgroundTESTFlagsPro = new TProfile("fBackgroundTESTFlagsPro","Flags and settings for test background",nTests,0,nTests);
 fBackgroundTESTFlagsPro->SetTickLength(-0.01,"Y");
 fBackgroundTESTFlagsPro->SetMarkerStyle(25);
 fBackgroundTESTFlagsPro->SetLabelSize(0.04);
 fBackgroundTESTFlagsPro->SetLabelOffset(0.02,"Y");
 fBackgroundTESTFlagsPro->SetStats(kFALSE);
 fBackgroundTESTFlagsPro->SetFillColor(kGray);
 fBackgroundTESTFlagsPro->SetLineColor(kBlack);
 for(Int_t t=0;t<nTests;t++)
 {
  fBackgroundTESTFlagsPro->GetXaxis()->SetBinLabel(t+1,Form("fFillBackgroundTEST[%d]",t)); fBackgroundTESTFlagsPro->Fill(t+0.5,fFillBackgroundTEST[t]);
 }
 fBackgroundTESTList->Add(fBackgroundTESTFlagsPro);

 // b) Book background objects for all tests:
 const Int_t nCumulants = 4; // TBI this applies only to vs. Q3
 const Int_t n3pCumulantTerms = 7;
 const Int_t nQs = 2; // Q2 and Q3 at the moment, if you change it here, change it also in exprs. like for(Int_t q=0;q<2;q++)
 const Int_t n2pCumulantTerms = 3;
 TString s2pCumulantTerms[n2pCumulantTerms] = {"#LTX_{1}#GT","#LTX_{2}#GT","#LTX_{1}X_{2}#GT"};
 TString s3pCumulantTerms[n3pCumulantTerms] = {"#LTX_{1}#GT","#LTX_{2}#GT","#LTX_{3}#GT","#LTX_{1}X_{2}#GT","#LTX_{1}X_{3}#GT","#LTX_{2}X_{3}#GT","#LTX_{1}X_{2}X_{3}#GT"};
 TString sCumulants[nCumulants] = {"#LTX_{1}X_{2}#GT_{c}","#LTX_{1}X_{3}#GT_{c}","#LTX_{2}X_{3}#GT_{c}","#LTX_{1}X_{2}X_{3}#GT_{c}"}; // TBI this applies only to vs. Q3
 TString sXYZ[3] = {"x","y","z"};
 TString sQs[nQs] = {"Q_{2}","Q_{3}"};
 for(Int_t t=0;t<nTests;t++)
 {
  if(fFillBackgroundTEST[t])
  {

   // Background correlations vs. Q2:
   for(Int_t ct=0;ct<n2pCumulantTerms;ct++)
   {
    for(Int_t xyz=0;xyz<3;xyz++)
    {
     fBackgroundTEST[t][0][ct][xyz] = new TProfile(Form("fBackgroundTEST[%d][0][%d][%d]",t,ct,xyz),Form("TEST: %d",t),fnQ2bins,fnQ2min,fnQ2max);
     fBackgroundTEST[t][0][ct][xyz]->SetStats(kFALSE);
     fBackgroundTEST[t][0][ct][xyz]->SetMarkerStyle(kFullSquare);
     fBackgroundTEST[t][0][ct][xyz]->SetMarkerColor(kBlue);
     fBackgroundTEST[t][0][ct][xyz]->SetLineColor(kBlue);
     fBackgroundTEST[t][0][ct][xyz]->GetXaxis()->SetTitle(sQs[0].Data());
     fBackgroundTEST[t][0][ct][xyz]->GetYaxis()->SetTitle(Form("%s_{%s}",s2pCumulantTerms[ct].Data(),sXYZ[xyz].Data()));
     fBackgroundTESTSublist[t]->Add(fBackgroundTEST[t][0][ct][xyz]);
    } // for(Int_t xyz=0;xyz<3;xyz++)
   } // for(Int_t ct=0;ct<n2pCumulantTerms;ct++)

   // Background cumulants vs. Q2:
   for(Int_t xyz=0;xyz<3;xyz++)
   {
    fBackgroundCumulantsTEST[t][0][0][xyz] = new TProfile(Form("fBackgroundCumulantsTEST[%d][0][%d][%d]",t,0,xyz),Form("TEST: %d",t),fnQ2bins,fnQ2min,fnQ2max);
    fBackgroundCumulantsTEST[t][0][0][xyz]->SetStats(kFALSE);
    fBackgroundCumulantsTEST[t][0][0][xyz]->SetMarkerStyle(kFullSquare);
    fBackgroundCumulantsTEST[t][0][0][xyz]->SetMarkerColor(kGreen+2);
    fBackgroundCumulantsTEST[t][0][0][xyz]->SetLineColor(kGreen+2);
    fBackgroundCumulantsTEST[t][0][0][xyz]->GetXaxis()->SetTitle(sQs[0].Data());
    fBackgroundCumulantsTEST[t][0][0][xyz]->GetYaxis()->SetTitle(Form("#LTX_{1}X_{2}#GT_{c,%s}",sXYZ[xyz].Data()));
    fBackgroundTESTSublist[t]->Add(fBackgroundCumulantsTEST[t][0][0][xyz]);
   } // for(Int_t xyz=0;xyz<3;xyz++)

   // Background correlations vs. Q3: TBI this clearly can be unified with the above.
   for(Int_t ct=0;ct<n3pCumulantTerms;ct++)
   {
    for(Int_t xyz=0;xyz<3;xyz++)
    {
     fBackgroundTEST[t][1][ct][xyz] = new TProfile(Form("fBackgroundTEST[%d][1][%d][%d]",t,ct,xyz),Form("TEST: %d",t),fnQ3bins,fnQ3min,fnQ3max);
     fBackgroundTEST[t][1][ct][xyz]->SetStats(kFALSE);
     fBackgroundTEST[t][1][ct][xyz]->SetMarkerStyle(kFullSquare);
     fBackgroundTEST[t][1][ct][xyz]->SetMarkerColor(kBlue);
     fBackgroundTEST[t][1][ct][xyz]->SetLineColor(kBlue);
     fBackgroundTEST[t][1][ct][xyz]->GetXaxis()->SetTitle(sQs[1].Data());
     fBackgroundTEST[t][1][ct][xyz]->GetYaxis()->SetTitle(Form("%s_{%s}",s3pCumulantTerms[ct].Data(),sXYZ[xyz].Data()));
     fBackgroundTESTSublist[t]->Add(fBackgroundTEST[t][1][ct][xyz]);
    } // for(Int_t xyz=0;xyz<3;xyz++)
   } // for(Int_t ct=0;ct<n3pCumulantTerms;ct++)

   // Background cumulants vs. Q3:
   for(Int_t ct=0;ct<nCumulants;ct++) // TBI this applies only to vs. Q3
   {
    for(Int_t xyz=0;xyz<3;xyz++)
    {
     fBackgroundCumulantsTEST[t][1][ct][xyz] = new TProfile(Form("fBackgroundCumulantsTEST[%d][1][%d][%d]",t,ct,xyz),Form("TEST: %d",t),fnQ3bins,fnQ3min,fnQ3max);
     fBackgroundCumulantsTEST[t][1][ct][xyz]->SetStats(kFALSE);
     fBackgroundCumulantsTEST[t][1][ct][xyz]->SetMarkerStyle(kFullSquare);
     fBackgroundCumulantsTEST[t][1][ct][xyz]->SetMarkerColor(kGreen+2);
     fBackgroundCumulantsTEST[t][1][ct][xyz]->SetLineColor(kGreen+2);
     fBackgroundCumulantsTEST[t][1][ct][xyz]->GetXaxis()->SetTitle(sQs[1].Data());
     fBackgroundCumulantsTEST[t][1][ct][xyz]->GetYaxis()->SetTitle(Form("%s_{%s}",sCumulants[ct].Data(),sXYZ[xyz].Data()));
     fBackgroundTESTSublist[t]->Add(fBackgroundCumulantsTEST[t][1][ct][xyz]);
    } // for(Int_t xyz=0;xyz<3;xyz++)
   } // for(Int_t ct=0;ct<nCumulants;ct++)
  } // if(fFillCorrelationFunctionsTEST[t])

 } // for(Int_t t=0;t<nTests;t++)


 // c) Book fMixedEventsTEST[3] and fGlobalTracksAODTEST[3]:
 for(Int_t me=0;me<3;me++) // [0] is buffer for 1st event; [1] for 2nd, etc.
 {
  fMixedEventsTEST[me] = new TClonesArray("AliAODTrack",10000);
  fGlobalTracksAODTEST[me] = new TExMap();
 }

 // d) Store temporary the yield as well:
 // vs. Q2:
 fBackgroundYieldTEST[0] = new TH1F("fBackgroundYieldTEST[0]","dN(X_{1}X_{2})/Q_{2} for mixed-events (TEST_0)",fnQ2bins,fnQ2min,fnQ2max);
 fBackgroundYieldTEST[0]->SetStats(kFALSE);
 fBackgroundYieldTEST[0]->SetLineColor(kRed);
 fBackgroundYieldTEST[0]->SetFillColor(kRed-10);
 fBackgroundYieldTEST[0]->GetXaxis()->SetTitle("Q_{2}");
 fBackgroundYieldTEST[0]->GetYaxis()->SetTitle("counts");
 fBackgroundTESTList->Add(fBackgroundYieldTEST[0]);
 // vs. Q3:
 fBackgroundYieldTEST[1] = new TH1F("fBackgroundYieldTEST[1]","dN(X_{1}X_{2}X_{3})/Q_{3} for mixed-events (TEST_1)",fnQ3bins,fnQ3min,fnQ3max);
 fBackgroundYieldTEST[1]->SetStats(kFALSE);
 fBackgroundYieldTEST[1]->SetLineColor(kRed);
 fBackgroundYieldTEST[1]->SetFillColor(kRed-10);
 fBackgroundYieldTEST[1]->GetXaxis()->SetTitle("Q_{3}");
 fBackgroundYieldTEST[1]->GetYaxis()->SetTitle("counts");
 fBackgroundTESTList->Add(fBackgroundYieldTEST[1]);

} // void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForBackgroundTEST()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForHybridApproach()
{
 // Book everything for 'hybrid approach'.

 // a) Book the profile holding all the flags;
 // b) Book fMixedEventsHA[3] and fGlobalTracksAODHA[3];
 // c) ...

 // a) Book the profile holding all the flags:
 fHybridApproachFlagsPro = new TProfile("fHybridApproachFlagsPro","Flags and settings for 'hybrid approach'",1,0,1);
 fHybridApproachFlagsPro->SetTickLength(-0.01,"Y");
 fHybridApproachFlagsPro->SetMarkerStyle(25);
 fHybridApproachFlagsPro->SetLabelSize(0.04);
 fHybridApproachFlagsPro->SetLabelOffset(0.02,"Y");
 fHybridApproachFlagsPro->SetStats(kFALSE);
 fHybridApproachFlagsPro->SetFillColor(kGray);
 fHybridApproachFlagsPro->SetLineColor(kBlack);
 fHybridApproachFlagsPro->GetXaxis()->SetBinLabel(1,"fDoHybridApproach"); fHybridApproachFlagsPro->Fill(0.5,fDoHybridApproach);
 fHybridApproachList->Add(fHybridApproachFlagsPro);

 if(!fDoHybridApproach){return;}

 TString sTerm[5] = {"N^{A}_{3}(p_{1},p_{2},p_{3})","N^{A}_{2}(p_{1},p_{2})N^{B}_{1}(p_{3})","N^{A}_{2}(p_{2},p_{3})N^{B}_{1}(p_{1})","N^{A}_{2}(p_{3},p_{1})N^{B}_{1}(p_{2})","N^{A}_{1}(p_{1})N^{B}_{1}(p_{2})N^{C}_{1}(p_{3})"};
 for(Int_t t=0;t<5;t++) // five distinct terms in the numerator of Eq. (18)
 {
  fDistributionsVsQ3[t] = new TH1F(Form("fDistributionsVsQ3[%d]",t),sTerm[t].Data(),5*fnQ3bins,fnQ3min,fnQ3max); // TBI hardwired refined binning, remove eventually
  //fDistributionsVsQ3[t]->SetStats(kFALSE);
  fDistributionsVsQ3[t]->SetMarkerStyle(kFullSquare);
  fDistributionsVsQ3[t]->SetMarkerColor(kBlue);
  fDistributionsVsQ3[t]->SetLineColor(kBlue);
  fDistributionsVsQ3[t]->SetFillColor(kBlue-10);
  fDistributionsVsQ3[t]->GetXaxis()->SetTitle("Q_{3}");
  fDistributionsVsQ3[t]->GetYaxis()->SetTitle("counts");
  fHybridApproachList->Add(fDistributionsVsQ3[t]);
 }

 // b) Book fMixedEventsHA[3] and fGlobalTracksAODHA[3]:
 for(Int_t me=0;me<3;me++) // [0] is buffer for 1st event; [1] for 2nd, etc.
 {
  fMixedEventsHA[me] = new TClonesArray("AliAODTrack",10000);
  fGlobalTracksAODHA[me] = new TExMap();
 }

 // c) ...

} // void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForHybridApproach()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::BookAndNestAllLists()
{
 // Book and nest all lists nested in the base list fHistList.

 // a) Book and nest lists for control histograms;
 // b) Book and nest lists for eveny-by-event histograms;
 // c) Book and nest lists for correlation functions;
 // d) Book and nest lists for background;
 // e) Book and nest lists for buffers;
 // f) Book and nest lists for QA;
 // g) Book and nest lists for common global track cuts;
 // h) Book and nest lists for test correlation functions;
 // i) Book and nest lists for test background;
 // j) Book and nest lists for 'hybrid approach'.

 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::BookAndNestAllLists()";
 if(!fHistList){Fatal(sMethodName.Data(),"fHistList is NULL");}

 // a) Book and nest lists for control histograms:
 fControlHistogramsList = new TList();
 fControlHistogramsList->SetName("Control_histograms");
 fControlHistogramsList->SetOwner(kTRUE);
 fHistList->Add(fControlHistogramsList);
 fControlHistogramsEventList = new TList();
 fControlHistogramsEventList->SetName("Event");
 fControlHistogramsEventList->SetOwner(kTRUE);
 fControlHistogramsList->Add(fControlHistogramsEventList);
 fControlHistogramsNonIdentifiedParticlesList = new TList();
 fControlHistogramsNonIdentifiedParticlesList->SetName("Non-identified_particles");
 fControlHistogramsNonIdentifiedParticlesList->SetOwner(kTRUE);
 fControlHistogramsList->Add(fControlHistogramsNonIdentifiedParticlesList);
 fControlHistogramsNonIdentifiedParticlesFTSFList = new TList();
 fControlHistogramsNonIdentifiedParticlesFTSFList->SetName("Non-identified_particles (f.t.s.f.)"); // for the specified filterbit
 fControlHistogramsNonIdentifiedParticlesFTSFList->SetOwner(kTRUE);
 fControlHistogramsList->Add(fControlHistogramsNonIdentifiedParticlesFTSFList);
 fControlHistogramsIdentifiedParticlesList = new TList();
 fControlHistogramsIdentifiedParticlesList->SetName("Identified_particles");
 fControlHistogramsIdentifiedParticlesList->SetOwner(kTRUE);
 fControlHistogramsList->Add(fControlHistogramsIdentifiedParticlesList);
 fControlHistogramsV0sList = new TList();
 fControlHistogramsV0sList->SetName("V0s");
 fControlHistogramsV0sList->SetOwner(kTRUE);
 fControlHistogramsList->Add(fControlHistogramsV0sList);

 // b) Book and nest lists for eveny-by-event histograms:
 fEBEHistogramsList = new TList();
 fEBEHistogramsList->SetName("Event-by-event_histograms");
 fEBEHistogramsList->SetOwner(kTRUE);
 fHistList->Add(fEBEHistogramsList);

 // c) Book and nest lists for correlation functions:
 fCorrelationFunctionsList = new TList();
 fCorrelationFunctionsList->SetName("CorrelationFunctions");
 fCorrelationFunctionsList->SetOwner(kTRUE);
 fHistList->Add(fCorrelationFunctionsList);
 // Correlation functions sublists:
 for(Int_t cfs=0;cfs<3;cfs++)
 {
  fCorrelationFunctionsSublist[cfs] = new TList();
  fCorrelationFunctionsSublist[cfs]->SetName(Form("f%dpCorrelationFunctions",cfs+2));
  fCorrelationFunctionsSublist[cfs]->SetOwner(kTRUE);
  fCorrelationFunctionsList->Add(fCorrelationFunctionsSublist[cfs]);
 } // for(Int_t cfs=0;cfs<3;cfs++)

 // d) Book and nest lists for background:
 fBackgroundList = new TList();
 fBackgroundList->SetName("Background");
 fBackgroundList->SetOwner(kTRUE);
 fHistList->Add(fBackgroundList);
 // Background sublists:
 for(Int_t bs=0;bs<3;bs++)
 {
  fBackgroundSublist[bs] = new TList();
  fBackgroundSublist[bs]->SetName(Form("f%dpBackground",bs+2));
  fBackgroundSublist[bs]->SetOwner(kTRUE);
  fBackgroundList->Add(fBackgroundSublist[bs]);
 } // for(Int_t bs=0;bs<3;bs++)

 // e) Book and nest lists for buffers:
 fBuffersList = new TList();
 fBuffersList->SetName("Buffers");
 fBuffersList->SetOwner(kTRUE);
 fHistList->Add(fBuffersList);

 // f) Book and nest lists for QA:
 fQAList = new TList();
 fQAList->SetName("QA");
 fQAList->SetOwner(kTRUE);
 fHistList->Add(fQAList);
 if(fFillQAEvents)
 {
  fQAEventsList = new TList();
  fQAEventsList->SetName("Events");
  fQAEventsList->SetOwner(kTRUE);
  fQAList->Add(fQAEventsList);
 }
 if(fFillQAParticles)
 {
  fQAParticlesList = new TList();
  fQAParticlesList->SetName("Particles");
  fQAParticlesList->SetOwner(kTRUE);
  fQAList->Add(fQAParticlesList);
 }

 // g) Book and nest lists for common global track cuts:
 fGlobalTrackCutsList = new TList();
 fGlobalTrackCutsList->SetName("Global_track_cuts");
 fGlobalTrackCutsList->SetOwner(kTRUE);
 fHistList->Add(fGlobalTrackCutsList);

 // h) Book and nest lists for test correlation functions:
 fCorrelationFunctionsTESTList = new TList();
 fCorrelationFunctionsTESTList->SetName("CorrelationFunctionsTEST");
 fCorrelationFunctionsTESTList->SetOwner(kTRUE);
 fHistList->Add(fCorrelationFunctionsTESTList);
 const Int_t nTests = 10;
 for(Int_t t=0;t<nTests;t++)
 {
  if(!fFillCorrelationFunctionsTEST[t]){continue;}
  fCorrelationFunctionsTESTSublist[t] = new TList();
  fCorrelationFunctionsTESTSublist[t]->SetName(Form("TEST_%d",t));
  fCorrelationFunctionsTESTSublist[t]->SetOwner(kTRUE);
  fCorrelationFunctionsTESTList->Add(fCorrelationFunctionsTESTSublist[t]);
 } // for(Int_t t=0;t<nTests;t++)

 // i) Book and nest lists for test background:
 fBackgroundTESTList = new TList();
 fBackgroundTESTList->SetName("BackgroundTEST");
 fBackgroundTESTList->SetOwner(kTRUE);
 fHistList->Add(fBackgroundTESTList);
 for(Int_t t=0;t<nTests;t++)
 {
  if(!fFillBackgroundTEST[t]){continue;}
  fBackgroundTESTSublist[t] = new TList();
  fBackgroundTESTSublist[t]->SetName(Form("TEST_%d",t));
  fBackgroundTESTSublist[t]->SetOwner(kTRUE);
  fBackgroundTESTList->Add(fBackgroundTESTSublist[t]);
 } // for(Int_t t=0;t<nTests;t++)

 // j) Book and nest lists for 'hybrid approach':
 fHybridApproachList = new TList();
 fHybridApproachList->SetName("HybridApproach");
 fHybridApproachList->SetOwner(kTRUE);
 fHistList->Add(fHybridApproachList);

} // void AliAnalysisTaskMultiparticleFemtoscopy::BookAndNestAllLists()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::BookEverything()
{
 // Book all unclassified objects temporary here. TBI

 fAnalysisType = new TString();

 fPIDResponse = new AliPIDResponse();

 for(Int_t index=0;index<fMaxNoGlobalTracksAOD;index++)
 {
  fGlobalTracksAOD[index] = new TExMap();
 }

} // void AliAnalysisTaskMultiparticleFemtoscopy::BookEverything()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForControlHistograms()
{
 // Book all the stuff for control histograms.

 // a) Book the profile holding all the flags for control histograms;
 // b) Common variables;
 // c) Book all control histograms...
 //  c0) Event;
 //  c1) Non-identified particles (for AOD these are "normal global" tracks);
 //  c2) Non-identified particles for the specified filterbit (f.t.s.f.) (by default TPC-only);
 //  c3) Identified particles;
 //  c4) V0s.

 // a) Book the profile holding all the flags for control histograms: TBI stil incomplete
 fControlHistogramsFlagsPro = new TProfile("fControlHistogramsFlagsPro","Flags and settings for control histograms",1,0,1);
 fControlHistogramsFlagsPro->SetTickLength(-0.01,"Y");
 fControlHistogramsFlagsPro->SetMarkerStyle(25);
 fControlHistogramsFlagsPro->SetLabelSize(0.04);
 fControlHistogramsFlagsPro->SetLabelOffset(0.02,"Y");
 fControlHistogramsFlagsPro->SetStats(kFALSE);
 fControlHistogramsFlagsPro->SetFillColor(kGray);
 fControlHistogramsFlagsPro->SetLineColor(kBlack);
 fControlHistogramsFlagsPro->GetXaxis()->SetBinLabel(1,"fFillControlHistograms"); fControlHistogramsFlagsPro->Fill(0.5,fFillControlHistograms);
 fControlHistogramsList->Add(fControlHistogramsFlagsPro);

 if(!fFillControlHistograms){return;}

 // b) Common variables:
 TString sParticleLabel[5] = {"e","#mu","#pi","K","p"};
 Double_t dNominalMass[5] = {TDatabasePDG::Instance()->GetParticle(11)->Mass(),TDatabasePDG::Instance()->GetParticle(13)->Mass(),TDatabasePDG::Instance()->GetParticle(211)->Mass(),TDatabasePDG::Instance()->GetParticle(321)->Mass(),TDatabasePDG::Instance()->GetParticle(2212)->Mass()};
 // ...

 //  c0) Event:
 // Book the profile holding all the flags for control histograms for global event observables:
 fControlHistogramsEventFlagsPro = new TProfile("fControlHistogramsEventFlagsPro","Flags and settings for TBI",1,0,1);
 fControlHistogramsEventFlagsPro->SetTickLength(-0.01,"Y");
 fControlHistogramsEventFlagsPro->SetMarkerStyle(25);
 fControlHistogramsEventFlagsPro->SetLabelSize(0.04);
 fControlHistogramsEventFlagsPro->SetLabelOffset(0.02,"Y");
 fControlHistogramsEventFlagsPro->SetStats(kFALSE);
 fControlHistogramsEventFlagsPro->SetFillColor(kGray);
 fControlHistogramsEventFlagsPro->SetLineColor(kBlack);
 fControlHistogramsEventFlagsPro->GetXaxis()->SetBinLabel(1,"fFillControlHistogramsEvent"); fControlHistogramsEventFlagsPro->Fill(0.5,fFillControlHistogramsEvent);
 fControlHistogramsEventList->Add(fControlHistogramsEventFlagsPro);
 if(fFillControlHistogramsEvent)
 {
  fGetNumberOfTracksHist = new TH1I("fGetNumberOfTracksHist","aAOD->GetNumberOfTracks() (Remark: Not all of tracks are unique.)",10000,0,10000);
  //fGetNumberOfTracksHist->SetStats(kFALSE);
  fGetNumberOfTracksHist->SetFillColor(kBlue-10);
  fControlHistogramsEventList->Add(fGetNumberOfTracksHist);
  fGetNumberOfGlobalTracksHist = new TH1I("fGetNumberOfGlobalTracksHist","fGlobalTracksAOD[0]->GetSize()",10000,0,10000);
  //fGetNumberOfGlobalTracksHist->SetStats(kFALSE);
  fGetNumberOfGlobalTracksHist->SetFillColor(kBlue-10);
  fControlHistogramsEventList->Add(fGetNumberOfGlobalTracksHist);
  fGetNumberOfV0sHist = new TH1I("fGetNumberOfV0sHist","aAOD->GetNumberOfV0s() (Remark: Some V0s share the daughter.)",10000,0,10000);
  fGetNumberOfV0sHist->SetStats(kFALSE);
  fGetNumberOfV0sHist->SetFillColor(kBlue-10);
  fControlHistogramsEventList->Add(fGetNumberOfV0sHist);
  fGetNumberOfCascadesHist = new TH1I("fGetNumberOfCascadesHist","aAOD->GetNumberOfCascades() (TBI: Not validated.)",10000,0,10000);
  fGetNumberOfCascadesHist->SetStats(kFALSE);
  fGetNumberOfCascadesHist->SetFillColor(kBlue-10);
  fControlHistogramsEventList->Add(fGetNumberOfCascadesHist);
  fGetMagneticFieldHist = new TH1D("fGetMagneticFieldHist","aAOD->GetMagneticField()",200,-10.,10.);
  fGetMagneticFieldHist->SetFillColor(kBlue-10);
  fControlHistogramsEventList->Add(fGetMagneticFieldHist);
  fGetEventTypeHist = new TH1I("fGetEventTypeHist","aAOD->GetEventType()",1000,0,10);
  fGetEventTypeHist->SetFillColor(kBlue-10);
  fControlHistogramsEventList->Add(fGetEventTypeHist);
  fGetCentralityHist = new TH1D("fGetCentralityHist","aAOD->GetCentrality()",100,0.,100.);
  fGetCentralityHist->SetFillColor(kBlue-10);
  fControlHistogramsEventList->Add(fGetCentralityHist);
  TString sxyz[3] = {"X","Y","Z"};
  for(Int_t xyz=0;xyz<3;xyz++)
  {
   fVertexXYZ[xyz] = new TH1F(Form("fVertex%s",sxyz[xyz].Data()),Form("avtz->Get%s()",sxyz[xyz].Data()),100000,-50.,50);
   fVertexXYZ[xyz]->SetStats(kFALSE);
   fControlHistogramsEventList->Add(fVertexXYZ[xyz]);
  }
  fGetNContributorsHist = new TH1I("fGetNContributorsHist","avtx->GetNContributors()",10000,0,10000);
  fGetNContributorsHist->SetStats(kFALSE);
  fGetNContributorsHist->SetFillColor(kBlue-10);
  fControlHistogramsEventList->Add(fGetNContributorsHist);
  fGetChi2perNDFHist = new TH1F("fGetChi2perNDFHist","avtx->GetChi2perNDF()",5000,0.,50.);
  fGetChi2perNDFHist->SetStats(kFALSE);
  fGetChi2perNDFHist->SetFillColor(kBlue-10);
  fControlHistogramsEventList->Add(fGetChi2perNDFHist);
  fGetNDaughtersHist = new TH1I("GetNDaughtersHist","avtx->GetNDaughters()",10000,0,10000);
  fGetNDaughtersHist->SetStats(kFALSE);
  fGetNDaughtersHist->SetFillColor(kBlue-10);
  fControlHistogramsEventList->Add(fGetNDaughtersHist);
  // ...
 } // if(fFillControlHistogramsEvent)

 //  c1) Non-identified particles:
 // Book the profile holding all the flags for control histograms for non-identified particles:
 fControlHistogramsNonIdentifiedParticlesFlagsPro = new TProfile("fControlHistogramsNonIdentifiedParticlesFlagsPro","Flags and settings for TBI",1,0,1);
 fControlHistogramsNonIdentifiedParticlesFlagsPro->SetTickLength(-0.01,"Y");
 fControlHistogramsNonIdentifiedParticlesFlagsPro->SetMarkerStyle(25);
 fControlHistogramsNonIdentifiedParticlesFlagsPro->SetLabelSize(0.04);
 fControlHistogramsNonIdentifiedParticlesFlagsPro->SetLabelOffset(0.02,"Y");
 fControlHistogramsNonIdentifiedParticlesFlagsPro->SetStats(kFALSE);
 fControlHistogramsNonIdentifiedParticlesFlagsPro->SetFillColor(kGray);
 fControlHistogramsNonIdentifiedParticlesFlagsPro->SetLineColor(kBlack);
 fControlHistogramsNonIdentifiedParticlesFlagsPro->GetXaxis()->SetBinLabel(1,"fFillControlHistogramsNonIdentifiedParticles"); fControlHistogramsNonIdentifiedParticlesFlagsPro->Fill(0.5,fFillControlHistogramsNonIdentifiedParticles);
 fControlHistogramsNonIdentifiedParticlesList->Add(fControlHistogramsNonIdentifiedParticlesFlagsPro);
 if(fFillControlHistogramsNonIdentifiedParticles)
 {
  fChargeHist = new TH1I("fChargeHist","atrack->Charge()",5,-2,3);
  fChargeHist->SetStats(kFALSE);
  fChargeHist->SetFillColor(kBlue-10);
  fChargeHist->GetXaxis()->SetBinLabel(1,"-2");
  fChargeHist->GetXaxis()->SetBinLabel(2,"-1");
  fChargeHist->GetXaxis()->SetBinLabel(3,"0");
  fChargeHist->GetXaxis()->SetBinLabel(4,"1");
  fChargeHist->GetXaxis()->SetBinLabel(5,"2");
  fControlHistogramsNonIdentifiedParticlesList->Add(fChargeHist);
  fGetTPCNclsHist = new TH1I("fGetTPCNclsHist","atrack->fGetTPCNclsHist()",200,0,200);
  fGetTPCNclsHist->SetStats(kFALSE);
  fGetTPCNclsHist->SetFillColor(kBlue-10);
  fGetTPCNclsHist->GetXaxis()->SetTitle("TPCNcls");
  fControlHistogramsNonIdentifiedParticlesList->Add(fGetTPCNclsHist);
  fGetTPCsignalNHist = new TH1I("fGetTPCsignalNHist","atrack->fGetTPCsignalNHist()",200,0,200);
  fGetTPCsignalNHist->SetStats(kFALSE);
  fGetTPCsignalNHist->SetFillColor(kBlue-10);
  fGetTPCsignalNHist->GetXaxis()->SetTitle("TPCsignalN");
  fControlHistogramsNonIdentifiedParticlesList->Add(fGetTPCsignalNHist);
  fGetITSNclsHist = new TH1I("fGetITSNclsHist","atrack->fGetITSNclsHist()",200,0,200);
  fGetITSNclsHist->SetStats(kFALSE);
  fGetITSNclsHist->SetFillColor(kBlue-10);
  fGetITSNclsHist->GetXaxis()->SetTitle("ITSNcls");
  fControlHistogramsNonIdentifiedParticlesList->Add(fGetITSNclsHist);
  fdEdxVsPtHist = new TH2F("fdEdxVsPtHist","atrack->GetTPCmomentum(),atrack->GetTPCsignal()",1000,0.,20.,1000,-500.,500.);
  fdEdxVsPtHist->SetStats(kFALSE);
  fControlHistogramsNonIdentifiedParticlesList->Add(fdEdxVsPtHist);
  fPtHist = new TH1F("fPtHist","atrack->Pt()",1000,0.,20.);
  fPtHist->SetStats(kFALSE);
  fPtHist->SetFillColor(kBlue-10);
  fPtHist->SetMinimum(0.);
  fControlHistogramsNonIdentifiedParticlesList->Add(fPtHist);
  fEtaHist = new TH1F("fEtaHist","atrack->Eta()",200,-2.,2.);
  fEtaHist->SetStats(kFALSE);
  fEtaHist->SetFillColor(kBlue-10);
  fEtaHist->SetMinimum(0.);
  fControlHistogramsNonIdentifiedParticlesList->Add(fEtaHist);
  fPhiHist = new TH1F("fPhiHist","atrack->Phi()",360,0.,TMath::TwoPi());
  fPhiHist->SetStats(kFALSE);
  fPhiHist->SetFillColor(kBlue-10);
  fPhiHist->SetMinimum(0.);
  fControlHistogramsNonIdentifiedParticlesList->Add(fPhiHist);
  // TBI
  fMassHist = new TH1F("fMassHist","atrack->M()",10000,0.,10.);
  fMassHist->SetStats(kFALSE);
  fMassHist->SetFillColor(kBlue-10);
  fMassHist->SetMinimum(0.);
  for(Int_t nm=0;nm<5;nm++) // nominal masses
  {
   fMassHist->GetXaxis()->SetBinLabel(fMassHist->FindBin(dNominalMass[nm]),Form("m_{%s}",sParticleLabel[nm].Data()));
  }
  fControlHistogramsNonIdentifiedParticlesList->Add(fMassHist);
  fGetFilterMap = new TH1I("fGetFilterMap","atrack->fGetFilterMap()",10000,0,10000);
  fGetFilterMap->SetStats(kFALSE);
  fGetFilterMap->SetFillColor(kBlue-10);
  fGetFilterMap->SetMinimum(0.);
  fControlHistogramsNonIdentifiedParticlesList->Add(fGetFilterMap);
  fGetPdgCode = new TH1I("fGetPdgCode","atrack->fGetPdgCode()",20000,-10000,10000);
  fGetPdgCode->SetStats(kFALSE);
  fGetPdgCode->SetFillColor(kBlue-10);
  fGetPdgCode->SetMinimum(0.);
  fControlHistogramsNonIdentifiedParticlesList->Add(fGetPdgCode);
 } // if(fFillControlHistogramsNonIdentifiedParticles)

 //  c2) Non-identified particles for the specified filterbit (f.t.s.f.) (by default TPC-only);
 // Book the profile holding all the flags for control histograms for non-identified particles for the specified filterbit (f.t.s.f.) (by default TPC-only);
 fControlHistogramsNonIdentifiedParticlesFTSFFlagsPro = new TProfile("fControlHistogramsNonIdentifiedParticlesFTSFFlagsPro","Flags and settings for TBI",2,0,2);
 fControlHistogramsNonIdentifiedParticlesFTSFFlagsPro->SetTickLength(-0.01,"Y");
 fControlHistogramsNonIdentifiedParticlesFTSFFlagsPro->SetMarkerStyle(25);
 fControlHistogramsNonIdentifiedParticlesFTSFFlagsPro->SetLabelSize(0.04);
 fControlHistogramsNonIdentifiedParticlesFTSFFlagsPro->SetLabelOffset(0.02,"Y");
 fControlHistogramsNonIdentifiedParticlesFTSFFlagsPro->SetStats(kFALSE);
 fControlHistogramsNonIdentifiedParticlesFTSFFlagsPro->SetFillColor(kGray);
 fControlHistogramsNonIdentifiedParticlesFTSFFlagsPro->SetLineColor(kBlack);
 fControlHistogramsNonIdentifiedParticlesFTSFFlagsPro->GetXaxis()->SetBinLabel(1,"fFillControlHistogramsNonIdentifiedParticlesFTSF"); fControlHistogramsNonIdentifiedParticlesFTSFFlagsPro->Fill(0.5,fFillControlHistogramsNonIdentifiedParticlesFTSF);
 fControlHistogramsNonIdentifiedParticlesFTSFFlagsPro->GetXaxis()->SetBinLabel(2,"fFilterBitFTSF"); fControlHistogramsNonIdentifiedParticlesFTSFFlagsPro->Fill(1.5,fFilterBitFTSF);
 fControlHistogramsNonIdentifiedParticlesFTSFList->Add(fControlHistogramsNonIdentifiedParticlesFTSFFlagsPro);
 if(fFillControlHistogramsNonIdentifiedParticlesFTSF)
 {
  fChargeFTSFHist = new TH1I("fChargeFTSFHist",Form("atrack->Charge(), fb = %d",fFilterBitFTSF),5,-2,3);
  fChargeFTSFHist->SetStats(kFALSE);
  fChargeFTSFHist->SetFillColor(kBlue-10);
  fChargeFTSFHist->GetXaxis()->SetBinLabel(1,"-2");
  fChargeFTSFHist->GetXaxis()->SetBinLabel(2,"-1");
  fChargeFTSFHist->GetXaxis()->SetBinLabel(3,"0");
  fChargeFTSFHist->GetXaxis()->SetBinLabel(4,"1");
  fChargeFTSFHist->GetXaxis()->SetBinLabel(5,"2");
  fControlHistogramsNonIdentifiedParticlesFTSFList->Add(fChargeFTSFHist);
  fGetTPCNclsFTSFHist = new TH1I("fGetTPCNclsFTSFHist",Form("atrack->fGetTPCNclsFTSFHist(), fb = %d",fFilterBitFTSF),200,0,200);
  fGetTPCNclsFTSFHist->SetStats(kFALSE);
  fGetTPCNclsFTSFHist->SetFillColor(kBlue-10);
  fGetTPCNclsFTSFHist->GetXaxis()->SetTitle("TPCNcls");
  fControlHistogramsNonIdentifiedParticlesFTSFList->Add(fGetTPCNclsFTSFHist);
  fGetTPCsignalNFTSFHist = new TH1I("fGetTPCsignalNFTSFHist",Form("atrack->fGetTPCsignalNFTSFHist(), fb = %d",fFilterBitFTSF),200,0,200);
  fGetTPCsignalNFTSFHist->SetStats(kFALSE);
  fGetTPCsignalNFTSFHist->SetFillColor(kBlue-10);
  fGetTPCsignalNFTSFHist->GetXaxis()->SetTitle("TPCsignalN");
  fControlHistogramsNonIdentifiedParticlesFTSFList->Add(fGetTPCsignalNFTSFHist);
  fGetITSNclsFTSFHist = new TH1I("fGetITSNclsFTSFHist",Form("atrack->fGetITSNclsFTSFHist(), fb = %d",fFilterBitFTSF),200,0,200);
  fGetITSNclsFTSFHist->SetStats(kFALSE);
  fGetITSNclsFTSFHist->SetFillColor(kBlue-10);
  fGetITSNclsFTSFHist->GetXaxis()->SetTitle("ITSNcls");
  fControlHistogramsNonIdentifiedParticlesFTSFList->Add(fGetITSNclsFTSFHist);
  fdEdxVsPtFTSFHist = new TH2F("fdEdxVsPtFTSFHist",Form("atrack->GetTPCmomentum(),atrack->GetTPCsignal(), fb = %d",fFilterBitFTSF),1000,0.,20.,1000,-500.,500.);
  fdEdxVsPtFTSFHist->SetStats(kFALSE);
  fControlHistogramsNonIdentifiedParticlesFTSFList->Add(fdEdxVsPtFTSFHist);
  fPtFTSFHist = new TH1F("fPtFTSFHist",Form("atrack->Pt(), fb = %d",fFilterBitFTSF),1000,0.,20.);
  fPtFTSFHist->SetStats(kFALSE);
  fPtFTSFHist->SetFillColor(kBlue-10);
  fPtFTSFHist->SetMinimum(0.);
  fControlHistogramsNonIdentifiedParticlesFTSFList->Add(fPtFTSFHist);
  fEtaFTSFHist = new TH1F("fEtaFTSFHist",Form("atrack->Eta(), fb = %d",fFilterBitFTSF),200,-2.,2.);
  fEtaFTSFHist->SetStats(kFALSE);
  fEtaFTSFHist->SetFillColor(kBlue-10);
  fEtaFTSFHist->SetMinimum(0.);
  fControlHistogramsNonIdentifiedParticlesFTSFList->Add(fEtaFTSFHist);
  fPhiFTSFHist = new TH1F("fPhiFTSFHist",Form("atrack->Phi(), fb = %d",fFilterBitFTSF),360,0.,TMath::TwoPi());
  fPhiFTSFHist->SetStats(kFALSE);
  fPhiFTSFHist->SetFillColor(kBlue-10);
  fPhiFTSFHist->SetMinimum(0.);
  fControlHistogramsNonIdentifiedParticlesFTSFList->Add(fPhiFTSFHist);
  // TBI
  fMassFTSFHist = new TH1F("fMassFTSFHist",Form("atrack->M(), fb = %d",fFilterBitFTSF),10000,0.,10.);
  fMassFTSFHist->SetStats(kFALSE);
  fMassFTSFHist->SetFillColor(kBlue-10);
  fMassFTSFHist->SetMinimum(0.);
  for(Int_t nm=0;nm<5;nm++) // nominal masses
  {
   fMassFTSFHist->GetXaxis()->SetBinLabel(fMassFTSFHist->FindBin(dNominalMass[nm]),Form("m_{%s}",sParticleLabel[nm].Data()));
  }
  fControlHistogramsNonIdentifiedParticlesFTSFList->Add(fMassFTSFHist);
  fGetFilterMap = new TH1I("fGetFilterMap",Form("atrack->fGetFilterMap(), fb = %d",fFilterBitFTSF),10000,0,10000);
  fGetFilterMap->SetStats(kFALSE);
  fGetFilterMap->SetFillColor(kBlue-10);
  fGetFilterMap->SetMinimum(0.);
  fControlHistogramsNonIdentifiedParticlesFTSFList->Add(fGetFilterMap);
  fGetPdgCode = new TH1I("fGetPdgCode",Form("atrack->fGetPdgCode(), fb = %d",fFilterBitFTSF),20000,-10000,10000);
  fGetPdgCode->SetStats(kFALSE);
  fGetPdgCode->SetFillColor(kBlue-10);
  fGetPdgCode->SetMinimum(0.);
  fControlHistogramsNonIdentifiedParticlesFTSFList->Add(fGetPdgCode);
 } // if(fFillControlFTSFHistogramsNonIdentifiedParticlesFTSF)

 //  c3) Identified particles:
 fControlHistogramsIdentifiedParticlesFlagsPro = new TProfile("fControlHistogramsIdentifiedParticlesFlagsPro","Flags and settings for identified particles",4,0,4);
 fControlHistogramsIdentifiedParticlesFlagsPro->SetTickLength(-0.01,"Y");
 fControlHistogramsIdentifiedParticlesFlagsPro->SetMarkerStyle(25);
 fControlHistogramsIdentifiedParticlesFlagsPro->SetLabelSize(0.04);
 fControlHistogramsIdentifiedParticlesFlagsPro->SetLabelOffset(0.02,"Y");
 fControlHistogramsIdentifiedParticlesFlagsPro->SetStats(kFALSE);
 fControlHistogramsIdentifiedParticlesFlagsPro->SetFillColor(kGray);
 fControlHistogramsIdentifiedParticlesFlagsPro->SetLineColor(kBlack);
 fControlHistogramsIdentifiedParticlesFlagsPro->GetXaxis()->SetBinLabel(1,"fFillControlHistogramsIdentifiedParticles"); fControlHistogramsIdentifiedParticlesFlagsPro->Fill(0.5,fFillControlHistogramsIdentifiedParticles);
 fControlHistogramsIdentifiedParticlesFlagsPro->GetXaxis()->SetBinLabel(2,"fUseDefaultInclusiveSigmaCuts"); fControlHistogramsIdentifiedParticlesFlagsPro->Fill(1.5,(Int_t)fUseDefaultInclusiveSigmaCuts);
 fControlHistogramsIdentifiedParticlesFlagsPro->GetXaxis()->SetBinLabel(3,"fUseDefaultExclusiveSigmaCuts"); fControlHistogramsIdentifiedParticlesFlagsPro->Fill(2.5,(Int_t)fUseDefaultExclusiveSigmaCuts);
 fControlHistogramsIdentifiedParticlesFlagsPro->GetXaxis()->SetBinLabel(4,"fFillControlHistogramsWithGlobalTrackInfo"); fControlHistogramsIdentifiedParticlesFlagsPro->Fill(3.5,(Int_t)fFillControlHistogramsWithGlobalTrackInfo);
 fControlHistogramsIdentifiedParticlesList->Add(fControlHistogramsIdentifiedParticlesFlagsPro);

 if(fFillControlHistogramsIdentifiedParticles)
 {
  // Default inclusive sigma cuts:
  if(fUseDefaultInclusiveSigmaCuts)
  {
   fInclusiveSigmaCuts[2] = 3.; // i.e. in function Pion(...) the inclusive cut for pions is 2 sigma
   fInclusiveSigmaCuts[3] = 3.; // i.e. in function Kaon(...) the inclusive cut for kaons is 2 sigma
   fInclusiveSigmaCuts[4] = 3.; // i.e. in function Proton(...) the inclusive cut for protons is 2 sigma
  }
  const Int_t nPidFunctions = 5; //
  TString sPidFunctions[nPidFunctions] = {"Electron(...)","Muon(...)","Pion(...)","Kaon(...)","Proton(...)"};
  const Int_t nParticleSpecies = 5;
  TString sParticleSpecies[nParticleSpecies] = {"e","#mu","#pi","K","p"};
  fInclusiveSigmaCutsPro = new TProfile("fInclusiveSigmaCutsPro","Inclusive sigma cuts",nPidFunctions,0.,nPidFunctions);
  fInclusiveSigmaCutsPro->GetYaxis()->SetTitle("n#sigma");
  for(Int_t pidFunction=0;pidFunction<5;pidFunction++)
  {
   fInclusiveSigmaCutsPro->SetTickLength(-0.01,"Y");
   fInclusiveSigmaCutsPro->SetMarkerStyle(25);
   fInclusiveSigmaCutsPro->SetLabelSize(0.04);
   fInclusiveSigmaCutsPro->SetLabelOffset(0.02,"Y");
   fInclusiveSigmaCutsPro->SetStats(kFALSE);
   fInclusiveSigmaCutsPro->SetFillColor(kGray);
   fInclusiveSigmaCutsPro->SetLineColor(kBlack);
   fInclusiveSigmaCutsPro->GetXaxis()->SetBinLabel(pidFunction+1,sPidFunctions[pidFunction].Data());
   fInclusiveSigmaCutsPro->Fill(pidFunction+0.5,fInclusiveSigmaCuts[pidFunction]);
  } // for(Int_t pidFunction=0;pidFunction<5;pidFunction++)
  fControlHistogramsIdentifiedParticlesList->Add(fInclusiveSigmaCutsPro);

  // Default exclusive sigma cuts:
  if(fUseDefaultExclusiveSigmaCuts)
  {
   // Pion(...)
   fExclusiveSigmaCuts[2][3] = 4.; // i.e. in function Pion(...) the exclusive cut for kaons is 4 sigma
   fExclusiveSigmaCuts[2][4] = 4.; // i.e. in function Pion(...) the exclusive cut for protons is 4 sigma
   // Kaon(...)
   fExclusiveSigmaCuts[3][2] = 4.; // i.e. in function Kaon(...) the exclusive cut for pions is 4 sigma
   fExclusiveSigmaCuts[3][4] = 4.; // i.e. in function Kaon(...) the exclusive cut for protons is 4 sigma
   // Proton(...)
   fExclusiveSigmaCuts[4][2] = 4.; // i.e. in function Proton(...) the exclusive cut for pions is 4 sigma
   fExclusiveSigmaCuts[4][3] = 4.; // i.e. in function Proton(...) the exclusive cut for kaons is 4 sigma
  } // if(fUseDefaultExclusiveSigmaCuts)
  fExclusiveSigmaCutsPro = new TProfile2D("fExclusiveSigmaCutsPro","Exclusive sigma cuts",nPidFunctions,0.,nPidFunctions,nParticleSpecies,0.,nParticleSpecies);
  fExclusiveSigmaCutsPro->SetTickLength(-0.01,"Y");
  fExclusiveSigmaCutsPro->SetMarkerStyle(25);
  fExclusiveSigmaCutsPro->SetLabelSize(0.04);
  //fExclusiveSigmaCutsPro->SetLabelOffset(0.01,"X");
  //fExclusiveSigmaCutsPro->SetLabelOffset(0.01,"Y");
  fExclusiveSigmaCutsPro->SetTitleOffset(0.9,"Z");
  fExclusiveSigmaCutsPro->SetStats(kFALSE);
  fExclusiveSigmaCutsPro->SetFillColor(kGray);
  fExclusiveSigmaCutsPro->SetLineColor(kBlack);
  fExclusiveSigmaCutsPro->GetZaxis()->SetTitle("n#sigma");
  for(Int_t pidFunction=0;pidFunction<nPidFunctions;pidFunction++)
  {
   fExclusiveSigmaCutsPro->GetXaxis()->SetBinLabel(pidFunction+1,sPidFunctions[pidFunction].Data());
   for(Int_t particleSpecies=0;particleSpecies<nParticleSpecies;particleSpecies++)
   {
    if(0==pidFunction){fExclusiveSigmaCutsPro->GetYaxis()->SetBinLabel(particleSpecies+1,sParticleSpecies[particleSpecies].Data());}
    fExclusiveSigmaCutsPro->Fill(pidFunction+0.5,particleSpecies+0.5,fExclusiveSigmaCuts[pidFunction][particleSpecies]);
   } // for(Int_t pidFunction=0;pidFunction<nPidFunctions;pidFunction++)
  } // for(Int_t pidFunction=0;pidFunction<nPidFunctions;pidFunction++)
  fControlHistogramsIdentifiedParticlesList->Add(fExclusiveSigmaCutsPro);

  for(Int_t pid=0;pid<5;pid++) // [0=e,1=mu,2=pi,3=K,4=p]
  {
   for(Int_t pa=0;pa<2;pa++) // particle(+q)/antiparticle(-q)
   {
    for(Int_t ps=0;ps<2;ps++) // kPrimary/kFromDecayVtx
    {
     if(fFillControlHistogramsWithGlobalTrackInfo)
     {
      fMassPIDHist[pid][pa][ps] = new TH1F(Form("fMassPIDHist[%d][%d][%d]",pid,pa,ps),Form("fMassPIDHist[%d][%d][%d] (%s) ('gtrack' parameters)",pid,pa,ps,sParticleLabel[pid].Data()),10000,0.,10.);
     }
     else
     {
      fMassPIDHist[pid][pa][ps] = new TH1F(Form("fMassPIDHist[%d][%d][%d]",pid,pa,ps),Form("fMassPIDHist[%d][%d][%d] (%s) ('atrack' parameters)",pid,pa,ps,sParticleLabel[pid].Data()),10000,0.,10.);
     }
     fMassPIDHist[pid][pa][ps]->GetXaxis()->SetTitle("m [GeV/c^{2}]");
     for(Int_t nm=0;nm<5;nm++) // nominal masses
     {
      fMassPIDHist[pid][pa][ps]->GetXaxis()->SetBinLabel(fMassPIDHist[pid][pa][ps]->FindBin(dNominalMass[nm]),Form("m_{%s}",sParticleLabel[nm].Data()));
     }
     fMassPIDHist[pid][pa][ps]->SetFillColor(kBlue-10);
     fControlHistogramsIdentifiedParticlesList->Add(fMassPIDHist[pid][pa][ps]);
     if(fFillControlHistogramsWithGlobalTrackInfo)
     {
      fPtPIDHist[pid][pa][ps] = new TH1F(Form("fPtPIDHist[%d][%d][%d]",pid,pa,ps),Form("fPtPIDHist[%d][%d][%d] ('gtrack' parameters)",pid,pa,ps),1000,0.,10.);
      fPPIDHist[pid][pa][ps][0] = new TH1F(Form("fPPIDHist[%d][%d][%d][0]",pid,pa,ps),Form("fPPIDHist[%d][%d][%d][0] ('gtrack' parameters)",pid,pa,ps),2000,-10.,10.);
      fPPIDHist[pid][pa][ps][1] = new TH1F(Form("fPPIDHist[%d][%d][%d][1]",pid,pa,ps),Form("fPPIDHist[%d][%d][%d][1] ('gtrack' parameters)",pid,pa,ps),2000,-10.,10.);
      fPPIDHist[pid][pa][ps][2] = new TH1F(Form("fPPIDHist[%d][%d][%d][2]",pid,pa,ps),Form("fPPIDHist[%d][%d][%d][2] ('gtrack' parameters)",pid,pa,ps),2000,-10.,10.);
     }
     else
     {
      fPtPIDHist[pid][pa][ps] = new TH1F(Form("fPtPIDHist[%d][%d][%d]",pid,pa,ps),Form("fPtPIDHist[%d][%d][%d] ('atrack' parameters)",pid,pa,ps),1000,0.,10.);
      fPPIDHist[pid][pa][ps][0] = new TH1F(Form("fPPIDHist[%d][%d][%d][0]",pid,pa,ps),Form("fPPIDHist[%d][%d][%d][0] ('atrack' parameters)",pid,pa,ps),2000,-10.,10.);
      fPPIDHist[pid][pa][ps][1] = new TH1F(Form("fPPIDHist[%d][%d][%d][1]",pid,pa,ps),Form("fPPIDHist[%d][%d][%d][1] ('atrack' parameters)",pid,pa,ps),2000,-10.,10.);
      fPPIDHist[pid][pa][ps][2] = new TH1F(Form("fPPIDHist[%d][%d][%d][2]",pid,pa,ps),Form("fPPIDHist[%d][%d][%d][2] ('atrack' parameters)",pid,pa,ps),2000,-10.,10.);
     }
     fPtPIDHist[pid][pa][ps]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
     fPtPIDHist[pid][pa][ps]->SetFillColor(kBlue-10);
     fControlHistogramsIdentifiedParticlesList->Add(fPtPIDHist[pid][pa][ps]);
     fPPIDHist[pid][pa][ps][0]->GetXaxis()->SetTitle("p_{x} [GeV/c]");
     fPPIDHist[pid][pa][ps][0]->SetFillColor(kBlue-10);
     fControlHistogramsIdentifiedParticlesList->Add(fPPIDHist[pid][pa][ps][0]);
     fPPIDHist[pid][pa][ps][1]->GetXaxis()->SetTitle("p_{y} [GeV/c]");
     fPPIDHist[pid][pa][ps][1]->SetFillColor(kBlue-10);
     fControlHistogramsIdentifiedParticlesList->Add(fPPIDHist[pid][pa][ps][1]);
     fPPIDHist[pid][pa][ps][2]->GetXaxis()->SetTitle("p_{z} [GeV/c]");
     fPPIDHist[pid][pa][ps][2]->SetFillColor(kBlue-10);
     fControlHistogramsIdentifiedParticlesList->Add(fPPIDHist[pid][pa][ps][2]);

     if(fFillControlHistogramsWithGlobalTrackInfo)
     {
      fEtaPIDHist[pid][pa][ps] = new TH1F(Form("fEtaPIDHist[%d][%d][%d]",pid,pa,ps),Form("fEtaPIDHist[%d][%d][%d] ('gtrack' parameters)",pid,pa,ps),200000,-2.,2.);
     }
     else
     {
      fEtaPIDHist[pid][pa][ps] = new TH1F(Form("fEtaPIDHist[%d][%d][%d]",pid,pa,ps),Form("fEtaPIDHist[%d][%d][%d] ('atrack' parameters)",pid,pa,ps),200000,-2.,2.);
     }
     fEtaPIDHist[pid][pa][ps]->SetFillColor(kBlue-10);
     fControlHistogramsIdentifiedParticlesList->Add(fEtaPIDHist[pid][pa][ps]);
     if(fFillControlHistogramsWithGlobalTrackInfo)
     {
      fPhiPIDHist[pid][pa][ps] = new TH1F(Form("fPhiPIDHist[%d][%d][%d]",pid,pa,ps),Form("fPhiPIDHist[%d][%d][%d] ('gtrack' parameters)",pid,pa,ps),360,0.,TMath::TwoPi());
     }
     else
     {
      fPhiPIDHist[pid][pa][ps] = new TH1F(Form("fPhiPIDHist[%d][%d][%d]",pid,pa,ps),Form("fPhiPIDHist[%d][%d][%d] ('atrack' parameters)",pid,pa,ps),360,0.,TMath::TwoPi());
     }
     fPhiPIDHist[pid][pa][ps]->GetXaxis()->SetTitle("#phi");
     fPhiPIDHist[pid][pa][ps]->SetFillColor(kBlue-10);
     fControlHistogramsIdentifiedParticlesList->Add(fPhiPIDHist[pid][pa][ps]);
    } // for(Int_t ps=0;ps<2;ps++) // kPrimary/kFromDecayVtx
   } // for(Int_t pa=0;pa<2;pa++) // particle(+q)/antiparticle(-q)
  } // for(Int_t pid=0;pid<4;pid++) // [0=e,1=mu,2=pi,3=K,4=p]
 } // if(fFillControlHistogramsIdentifiedParticles)

 //  c4) V0s:
 // Book the profile holding all the flags for V0s:
 fControlHistogramsV0sFlagsPro = new TProfile("fControlHistogramsV0sFlagsPro","Flags and settings for V0s",1,0,1);
 fControlHistogramsV0sFlagsPro->SetTickLength(-0.01,"Y");
 fControlHistogramsV0sFlagsPro->SetMarkerStyle(25);
 fControlHistogramsV0sFlagsPro->SetLabelSize(0.04);
 fControlHistogramsV0sFlagsPro->SetLabelOffset(0.02,"Y");
 fControlHistogramsV0sFlagsPro->SetStats(kFALSE);
 fControlHistogramsV0sFlagsPro->SetFillColor(kGray);
 fControlHistogramsV0sFlagsPro->SetLineColor(kBlack);
 fControlHistogramsV0sFlagsPro->GetXaxis()->SetBinLabel(1,"fFillControlHistogramsV0s"); fControlHistogramsV0sFlagsPro->Fill(0.5,fFillControlHistogramsV0s);
 fControlHistogramsV0sList->Add(fControlHistogramsV0sFlagsPro);
 if(fFillControlHistogramsV0s)
 {
  fGetNProngsHist = new TH1I("fGetNProngsHist","aAODv0->GetNProngs()",10,0,10);
  fGetNProngsHist->SetStats(kFALSE);
  fGetNProngsHist->SetFillColor(kBlue-10);
  fControlHistogramsV0sList->Add(fGetNProngsHist);
  // TBI
  fMassK0ShortHist = new TH1F("fMassK0ShortHist","aAODv0->MassK0Short()",1000000,0.,100.);
  //fMassK0ShortHist->SetStats(kFALSE);
  fMassK0ShortHist->SetFillColor(kBlue-10);
  Double_t dMassK0Short = TDatabasePDG::Instance()->GetParticle(310)->Mass(); // nominal mass
  //fMassK0ShortHist->GetXaxis()->SetBinLabel(fMassK0ShortHist->FindBin(dMassK0Short),Form("m_{K_{S}^{0}} = %f",dMassK0Short));
  fMassK0ShortHist->SetBinContent(fMassK0ShortHist->FindBin(dMassK0Short),1e6);
  fControlHistogramsV0sList->Add(fMassK0ShortHist);
  // TBI
  fMassLambdaHist = new TH1F("fMassLambdaHist","aAODv0->MassLambda()",1000000,0.,100.);
  //fMassLambdaHist->SetStats(kFALSE);
  fMassLambdaHist->SetFillColor(kBlue-10);
  Double_t dMassLambda = TDatabasePDG::Instance()->GetParticle(3122)->Mass(); // nominal mass
  //fMassLambdaHist->GetXaxis()->SetBinLabel(fMassLambdaHist->FindBin(dMassLambda),Form("m_{#Lambda^{0}} = %f",dMassLambda));
  fMassLambdaHist->SetBinContent(fMassLambdaHist->FindBin(dMassLambda),1e6);
  fControlHistogramsV0sList->Add(fMassLambdaHist);
  // TBI
  fMassAntiLambdaHist = new TH1F("fMassAntiLambdaHist","aAODv0->MassAntiLambda()",1000000,0.,100.);
  //fMassAntiLambdaHist->SetStats(kFALSE);
  fMassAntiLambdaHist->SetFillColor(kBlue-10);
  Double_t dMassAntiLambda = TDatabasePDG::Instance()->GetParticle(-3122)->Mass(); // nominal mass
  //fMassAntiLambdaHist->GetXaxis()->SetBinLabel(fMassAntiLambdaHist->FindBin(dMassAntiLambda),Form("m_{#bar{Lambda}^{0}} = %f",dMassAntiLambda));
  fMassAntiLambdaHist->SetBinContent(fMassAntiLambdaHist->FindBin(dMassAntiLambda),1e6);
  fControlHistogramsV0sList->Add(fMassAntiLambdaHist);
  // TBI
  fOpenAngleV0Hist = new TH1F("fOpenAngleV0Hist","aAODv0->fOpenAngleV0()",10000,-0.044,TMath::Pi()+0.044);
  fOpenAngleV0Hist->SetStats(kFALSE);
  fOpenAngleV0Hist->SetFillColor(kBlue-10);
  fControlHistogramsV0sList->Add(fOpenAngleV0Hist);
  // TBI
  fRadiusV0Hist = new TH1F("fRadiusV0Hist","aAODv0->fRadiusV0()",10000,0.,1000.);
  fRadiusV0Hist->SetStats(kFALSE);
  fRadiusV0Hist->SetFillColor(kBlue-10);
  fControlHistogramsV0sList->Add(fRadiusV0Hist);
  // TBI
  fDcaV0ToPrimVertexHist = new TH1F("fDcaV0ToPrimVertexHist","aAODv0->fDcaV0ToPrimVertex()",10000,0.,1000.);
  fDcaV0ToPrimVertexHist->SetStats(kFALSE);
  fDcaV0ToPrimVertexHist->SetFillColor(kBlue-10);
  fControlHistogramsV0sList->Add(fDcaV0ToPrimVertexHist);
  // TBI
  fMomV0XHist = new TH1F("fMomV0XHist","aAODv0->fMomV0X() = px(+) + px(-)",10000,-1000.,1000.);
  fMomV0XHist->SetStats(kFALSE);
  fMomV0XHist->SetFillColor(kBlue-10);
  fControlHistogramsV0sList->Add(fMomV0XHist);
  // TBI
  fMomV0YHist = new TH1F("fMomV0YHist","aAODv0->fMomV0Y() = py(+) + py(-)",10000,-1000.,1000.);
  fMomV0YHist->SetStats(kFALSE);
  fMomV0YHist->SetFillColor(kBlue-10);
  fControlHistogramsV0sList->Add(fMomV0YHist);
  // TBI
  fMomV0ZHist = new TH1F("fMomV0ZHist","aAODv0->fMomV0Z() = pz(+) + pz(-)",10000,-1000.,1000.);
  fMomV0ZHist->SetStats(kFALSE);
  fMomV0ZHist->SetFillColor(kBlue-10);
  fControlHistogramsV0sList->Add(fMomV0ZHist);
  // TBI
  fPtV0Hist = new TH1F("fPtV0Hist","pow(aAODv0->fPt2V0(),0.5)",10000,0.,100.);
  fPtV0Hist->SetStats(kFALSE);
  fPtV0Hist->SetFillColor(kBlue-10);
  fControlHistogramsV0sList->Add(fPtV0Hist);
  // TBI
  fPseudoRapV0Hist = new TH1F("fPseudoRapV0Hist","aAODv0->PseudoRapV0()",1000,-10.,10.);
  fPseudoRapV0Hist->SetStats(kFALSE);
  fPseudoRapV0Hist->SetFillColor(kBlue-10);
  fControlHistogramsV0sList->Add(fPseudoRapV0Hist);
  // TBI
  fPAHist = new TH2F("fPAHist","TBI",100,-2.,2.,100,0.,1.);
  fPAHist->SetStats(kFALSE);
  fPAHist->GetXaxis()->SetTitle("#alpha");
  fPAHist->GetYaxis()->SetTitle("p_{T}");
  fControlHistogramsV0sList->Add(fPAHist);
 } // if(fFillControlHistogramsV0s)

} // void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForControlHistograms()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForEBEObjects()
{
 // Book all the stuff for event-by-event objects.

 // a) Book the profile holding all the flags for event-by-event objects;
 // b) Book all event-by-event objects.

 // a) Book the profile holding all the flags for EBE objects:
 fEBEObjectsFlagsPro = new TProfile("fEBEObjectsFlagsPro","Flags and settings for event-by-event histograms",1,0,1);
 fEBEObjectsFlagsPro->SetTickLength(-0.01,"Y");
 fEBEObjectsFlagsPro->SetMarkerStyle(25);
 fEBEObjectsFlagsPro->SetLabelSize(0.04);
 fEBEObjectsFlagsPro->SetLabelOffset(0.02,"Y");
 fEBEObjectsFlagsPro->SetStats(kFALSE);
 fEBEObjectsFlagsPro->SetFillColor(kGray);
 fEBEObjectsFlagsPro->SetLineColor(kBlack);
 //fEBEObjectsFlagsPro->GetXaxis()->SetBinLabel(1,"fFillEBEHistograms"); fEBEObjectsFlagsPro->Fill(0.5,fFillEBEHistograms);
 fEBEHistogramsList->Add(fEBEObjectsFlagsPro);

 //if(!fFillEBEHistograms){return;} // TBI rethink

 // TBI 
 fUniqueIDHistEBE = new TH1I("fUniqueIDHistEBE","TBI",40000,-20000,20000);
 fUniqueIDHistEBE->SetStats(kFALSE);
 fUniqueIDHistEBE->SetFillColor(kBlue-10);
 // TBI I do not want to store this histogram, right?

 // TBI add comment
 for(Int_t pid=0;pid<5;pid++) // [0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pa=0;pa<2;pa++) // particle(+q)/antiparticle(-q)
  {
   for(Int_t ps=0;ps<2;ps++) // kPrimary/kFromDecayVtx
   {
    fPIDCA[pid][pa][ps] = new TClonesArray("AliAODTrack",10000);
   }
  }
 }

 // TBI add comment
 for(Int_t pid=0;pid<1;pid++) // [0=Lambda,1=...] 
 {
  fPIDV0sCA[pid] = new TClonesArray("AliAODv0",10000);
 }

} // void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForEBEObjects()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForCorrelationFunctions()
{
 // Book all the stuff for correlation functions objects.

 // a) Book the profile holding all the flags for correlation functions objects;
 // b) Book all 2p correlation functions;
 // c) Book TExMap *fCorrelationFunctionsIndices;
 // d) Book all 3p correlation functions;
 // e) Book all 4p correlation functions.

 // a) Book the profile holding all the flags for correlation functions objects:
 fCorrelationFunctionsFlagsPro = new TProfile("fCorrelationFunctionsFlagsPro","Flags and settings for correlation functions histograms",8,0,8);
 fCorrelationFunctionsFlagsPro->SetTickLength(-0.01,"Y");
 fCorrelationFunctionsFlagsPro->SetMarkerStyle(25);
 fCorrelationFunctionsFlagsPro->SetLabelSize(0.04);
 fCorrelationFunctionsFlagsPro->SetLabelOffset(0.02,"Y");
 fCorrelationFunctionsFlagsPro->SetStats(kFALSE);
 fCorrelationFunctionsFlagsPro->SetFillColor(kGray);
 fCorrelationFunctionsFlagsPro->SetLineColor(kBlack);
 fCorrelationFunctionsFlagsPro->GetXaxis()->SetBinLabel(1,"fFillCorrelationFunctions"); fCorrelationFunctionsFlagsPro->Fill(0.5,fFillCorrelationFunctions);
 fCorrelationFunctionsFlagsPro->GetXaxis()->SetBinLabel(2,"fNormalizeCorrelationFunctions"); fCorrelationFunctionsFlagsPro->Fill(1.5,fNormalizeCorrelationFunctions);
 fCorrelationFunctionsFlagsPro->GetXaxis()->SetBinLabel(3,"fFill3pCorrelationFunctions"); fCorrelationFunctionsFlagsPro->Fill(2.5,fFill3pCorrelationFunctions);
 fCorrelationFunctionsFlagsPro->GetXaxis()->SetBinLabel(4,"fFill4pCorrelationFunctions"); fCorrelationFunctionsFlagsPro->Fill(3.5,fFill4pCorrelationFunctions);
 fCorrelationFunctionsFlagsPro->GetXaxis()->SetBinLabel(5,"fNormalizationOption"); fCorrelationFunctionsFlagsPro->Fill(4.5,fNormalizationOption);
 fCorrelationFunctionsFlagsPro->GetXaxis()->SetBinLabel(6,"fNormalizationInterval[0]"); fCorrelationFunctionsFlagsPro->Fill(5.5,fNormalizationInterval[0]);
 fCorrelationFunctionsFlagsPro->GetXaxis()->SetBinLabel(7,"fNormalizationInterval[1]"); fCorrelationFunctionsFlagsPro->Fill(6.5,fNormalizationInterval[1]);
 fCorrelationFunctionsFlagsPro->GetXaxis()->SetBinLabel(8,"fnMergedBins"); fCorrelationFunctionsFlagsPro->Fill(7.5,fnMergedBins);
 fCorrelationFunctionsList->Add(fCorrelationFunctionsFlagsPro);
 // 2p:
 f2pCorrelationFunctionsFlagsPro = new TProfile("f2pCorrelationFunctionsFlagsPro","Flags and settings for correlation functions histograms",1,0,1);
 f2pCorrelationFunctionsFlagsPro->SetTickLength(-0.01,"Y");
 f2pCorrelationFunctionsFlagsPro->SetMarkerStyle(25);
 f2pCorrelationFunctionsFlagsPro->SetLabelSize(0.04);
 f2pCorrelationFunctionsFlagsPro->SetLabelOffset(0.02,"Y");
 f2pCorrelationFunctionsFlagsPro->SetStats(kFALSE);
 f2pCorrelationFunctionsFlagsPro->SetFillColor(kGray);
 f2pCorrelationFunctionsFlagsPro->SetLineColor(kBlack);
 f2pCorrelationFunctionsFlagsPro->GetXaxis()->SetBinLabel(1,"TBI"); f2pCorrelationFunctionsFlagsPro->Fill(0.5,0.); // TBI dummy at the moment
 fCorrelationFunctionsSublist[0]->Add(f2pCorrelationFunctionsFlagsPro);
 // 3p:
 f3pCorrelationFunctionsFlagsPro = new TProfile("f3pCorrelationFunctionsFlagsPro","Flags and settings for correlation functions histograms",1,0,1);
 f3pCorrelationFunctionsFlagsPro->SetTickLength(-0.01,"Y");
 f3pCorrelationFunctionsFlagsPro->SetMarkerStyle(25);
 f3pCorrelationFunctionsFlagsPro->SetLabelSize(0.04);
 f3pCorrelationFunctionsFlagsPro->SetLabelOffset(0.02,"Y");
 f3pCorrelationFunctionsFlagsPro->SetStats(kFALSE);
 f3pCorrelationFunctionsFlagsPro->SetFillColor(kGray);
 f3pCorrelationFunctionsFlagsPro->SetLineColor(kBlack);
 f3pCorrelationFunctionsFlagsPro->GetXaxis()->SetBinLabel(1,"TBI"); f3pCorrelationFunctionsFlagsPro->Fill(0.5,0.); // TBI dummy at the moment
 fCorrelationFunctionsSublist[1]->Add(f3pCorrelationFunctionsFlagsPro);
 // 4p:
 f4pCorrelationFunctionsFlagsPro = new TProfile("f4pCorrelationFunctionsFlagsPro","Flags and settings for correlation functions histograms",1,0,1);
 f4pCorrelationFunctionsFlagsPro->SetTickLength(-0.01,"Y");
 f4pCorrelationFunctionsFlagsPro->SetMarkerStyle(25);
 f4pCorrelationFunctionsFlagsPro->SetLabelSize(0.04);
 f4pCorrelationFunctionsFlagsPro->SetLabelOffset(0.02,"Y");
 f4pCorrelationFunctionsFlagsPro->SetStats(kFALSE);
 f4pCorrelationFunctionsFlagsPro->SetFillColor(kGray);
 f4pCorrelationFunctionsFlagsPro->SetLineColor(kBlack);
 f4pCorrelationFunctionsFlagsPro->GetXaxis()->SetBinLabel(1,"TBI"); f4pCorrelationFunctionsFlagsPro->Fill(0.5,0.); // TBI dummy at the moment
 fCorrelationFunctionsSublist[2]->Add(f4pCorrelationFunctionsFlagsPro);

 if(!fFillCorrelationFunctions){return;} // TBI is this safe? It is not, because now if I want to fill 3-p, I always have to fill also 2-p

 const Int_t nParticleSpecies = 5; // Supported at the moment: 0=e,1=mu,2=pi,3=K,4=p
 TString sParticles[2*nParticleSpecies] = {"e^{+}","#mu^{+}","#pi^{+}","K^{+}","p^{+}","e^{-}","#mu^{-}","#pi^{-}","K^{-}","p^{-}"};

 // b) Book all 2p correlation functions:
 // Remark 0: First particle in the pair is always the one with positive charge.
 // Remark 1: Diagonal elements are particle/antiparticle pairs.
 for(Int_t pid1=0;pid1<2*nParticleSpecies;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pid2=pid1;pid2<2*nParticleSpecies;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   // Correlation functions:
   fCorrelationFunctions[pid1][pid2] = new TH1F(Form("fCorrelationFunctions[%d][%d]",pid1,pid2),Form("fCorrelationFunctions[%d][%d] = (%s,%s)",pid1,pid2,sParticles[pid1].Data(),sParticles[pid2].Data()),10000,0.,10.); // TBI rethink the boundaries and nbins
   fCorrelationFunctions[pid1][pid2]->SetStats(kTRUE);
   fCorrelationFunctions[pid1][pid2]->SetFillColor(kBlue-10);
   fCorrelationFunctions[pid1][pid2]->SetXTitle("Q_{2}");
   fCorrelationFunctions[pid1][pid2]->SetYTitle("C(Q_{2})");
   fCorrelationFunctionsSublist[0]->Add(fCorrelationFunctions[pid1][pid2]);
  } // for(Int_t pid2=0;pid2<5;pid2++) // anti-particle [0=e,1=mu,2=pi,3=K,4=p]
 } // for(Int_t pid=0;pid<5;pid++) // particle [0=e,1=mu,2=pi,3=K,4=p]

 // c) Book TExMap *fCorrelationFunctionsIndices:
 //    Chosen convention: [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 5=e,6=mu,7=pi,8=K,9=p]
 fCorrelationFunctionsIndices = new TExMap();
 fCorrelationFunctionsIndices->Add(11,0);
 fCorrelationFunctionsIndices->Add(13,1);
 fCorrelationFunctionsIndices->Add(211,2);
 fCorrelationFunctionsIndices->Add(321,3);
 fCorrelationFunctionsIndices->Add(2212,4);
 fCorrelationFunctionsIndices->Add(-11,5);
 fCorrelationFunctionsIndices->Add(-13,6);
 fCorrelationFunctionsIndices->Add(-211,7);
 fCorrelationFunctionsIndices->Add(-321,8);
 fCorrelationFunctionsIndices->Add(-2212,9);

 // d) Book all 3p correlation functions:
 if(fFill3pCorrelationFunctions)
 {
  for(Int_t pid1=0;pid1<2*nParticleSpecies;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   for(Int_t pid2=0;pid2<2*nParticleSpecies;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
   {
    for(Int_t pid3=0;pid3<2*nParticleSpecies;pid3++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
    {
     // 3-p correlation functions:
     f3pCorrelationFunctions[pid1][pid2][pid3] = new TH1F(Form("f3pCorrelationFunctions[%d][%d][%d]",pid1,pid2,pid3),Form("f3pCorrelationFunctions[%d][%d][%d] = (%s,%s,%s)",pid1,pid2,pid3,sParticles[pid1].Data(),sParticles[pid2].Data(),sParticles[pid3].Data()),10000,0.,10.); // TBI rethink the boundaries and nbins
     f3pCorrelationFunctions[pid1][pid2][pid3]->SetStats(kTRUE);
     f3pCorrelationFunctions[pid1][pid2][pid3]->SetFillColor(kBlue-10);
     f3pCorrelationFunctions[pid1][pid2][pid3]->SetXTitle("Q_{3}");
     f3pCorrelationFunctions[pid1][pid2][pid3]->SetYTitle("C(Q_{3})");
     fCorrelationFunctionsSublist[1]->Add(f3pCorrelationFunctions[pid1][pid2][pid3]);
    } // for(Int_t pid2=0;pid2<5;pid2++) // anti-particle [0=e,1=mu,2=pi,3=K,4=p]
   } // for(Int_t pid2=0;pid2<2*nParticleSpecies;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  } // for(Int_t pid=0;pid<5;pid++) // particle [0=e,1=mu,2=pi,3=K,4=p]
 } // if(fFill3pCorrelationFunctions)

 // e) Book all 4p correlation functions:
 if(fFill4pCorrelationFunctions)
 {
  for(Int_t pid1=0;pid1<2*nParticleSpecies;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   for(Int_t pid2=0;pid2<2*nParticleSpecies;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
   {
    for(Int_t pid3=0;pid3<2*nParticleSpecies;pid3++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
    {
     for(Int_t pid4=0;pid4<2*nParticleSpecies;pid4++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
     {
      // 4-p correlation functions:
      f4pCorrelationFunctions[pid1][pid2][pid3][pid4] = new TH1F(Form("f4pCorrelationFunctions[%d][%d][%d][%d]",pid1,pid2,pid3,pid4),Form("f4pCorrelationFunctions[%d][%d][%d][%d] = (%s,%s,%s,%s)",pid1,pid2,pid3,pid4,sParticles[pid1].Data(),sParticles[pid2].Data(),sParticles[pid3].Data(),sParticles[pid4].Data()),10000,0.,10.); // TBI rethink the boundaries and nbins
      f4pCorrelationFunctions[pid1][pid2][pid3][pid4]->SetStats(kTRUE);
      f4pCorrelationFunctions[pid1][pid2][pid3][pid4]->SetFillColor(kBlue-10);
      f4pCorrelationFunctions[pid1][pid2][pid3][pid4]->SetXTitle("Q_{4}");
      f4pCorrelationFunctions[pid1][pid2][pid3][pid4]->SetYTitle("C(Q_{4})");
      fCorrelationFunctionsSublist[2]->Add(f4pCorrelationFunctions[pid1][pid2][pid3][pid4]);
     } // for(Int_t pid4=0;pid4<2*nParticleSpecies;pid4++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
    } // for(Int_t pid3=0;pid3<2*nParticleSpecies;pid3++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
   } // for(Int_t pid2=0;pid2<2*nParticleSpecies;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  } // for(Int_t pid1=0;pid1<2*nParticleSpecies;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 } // if(fFill4pCorrelationFunctions)

} // void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForCorrelationFunctions()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForBackground()
{
 // Book all the stuff for background objects.

 // a) Book the profile holding all the flags for background objects;
 // b) Book all histograms for 2p background;
 // c) Book all histograms for 3p background;
 // d) Book all histograms for 4p background;
 // e) Book buffer objects.

 // a) Book the profile holding all the flags for background objects:
 fBackgroundFlagsPro = new TProfile("fBackgroundFlagsPro","Flags and settings for background histograms",4,0,4);
 fBackgroundFlagsPro->SetTickLength(-0.01,"Y");
 fBackgroundFlagsPro->SetMarkerStyle(25);
 fBackgroundFlagsPro->SetLabelSize(0.04);
 fBackgroundFlagsPro->SetLabelOffset(0.02,"Y");
 fBackgroundFlagsPro->SetStats(kFALSE);
 fBackgroundFlagsPro->SetFillColor(kGray);
 fBackgroundFlagsPro->SetLineColor(kBlack);
 fBackgroundFlagsPro->GetXaxis()->SetBinLabel(1,"fBackgroundOption"); fBackgroundFlagsPro->Fill(0.5,fBackgroundOption);
 fBackgroundFlagsPro->GetXaxis()->SetBinLabel(2,"fEstimate2pBackground"); fBackgroundFlagsPro->Fill(1.5,fEstimate2pBackground);
 fBackgroundFlagsPro->GetXaxis()->SetBinLabel(3,"fEstimate3pBackground"); fBackgroundFlagsPro->Fill(2.5,fEstimate3pBackground);
 fBackgroundFlagsPro->GetXaxis()->SetBinLabel(4,"fEstimate4pBackground"); fBackgroundFlagsPro->Fill(3.5,fEstimate4pBackground);
 fBackgroundList->Add(fBackgroundFlagsPro);
 // 2p:
 f2pBackgroundFlagsPro = new TProfile("f2pBackgroundFlagsPro","Flags and settings for 2p background histograms",1,0,1);
 f2pBackgroundFlagsPro->SetTickLength(-0.01,"Y");
 f2pBackgroundFlagsPro->SetMarkerStyle(25);
 f2pBackgroundFlagsPro->SetLabelSize(0.04);
 f2pBackgroundFlagsPro->SetLabelOffset(0.02,"Y");
 f2pBackgroundFlagsPro->SetStats(kFALSE);
 f2pBackgroundFlagsPro->SetFillColor(kGray);
 f2pBackgroundFlagsPro->SetLineColor(kBlack);
 f2pBackgroundFlagsPro->GetXaxis()->SetBinLabel(1,"TBI"); f2pBackgroundFlagsPro->Fill(0.5,0); // TBI dummy at the moment
 fBackgroundSublist[0]->Add(f2pBackgroundFlagsPro);
 // 3p:
 f3pBackgroundFlagsPro = new TProfile("f3pBackgroundFlagsPro","Flags and settings for 3p background histograms",1,0,1);
 f3pBackgroundFlagsPro->SetTickLength(-0.01,"Y");
 f3pBackgroundFlagsPro->SetMarkerStyle(25);
 f3pBackgroundFlagsPro->SetLabelSize(0.04);
 f3pBackgroundFlagsPro->SetLabelOffset(0.02,"Y");
 f3pBackgroundFlagsPro->SetStats(kFALSE);
 f3pBackgroundFlagsPro->SetFillColor(kGray);
 f3pBackgroundFlagsPro->SetLineColor(kBlack);
 f3pBackgroundFlagsPro->GetXaxis()->SetBinLabel(1,"TBI"); f3pBackgroundFlagsPro->Fill(0.5,0); // TBI dummy at the moment
 fBackgroundSublist[1]->Add(f3pBackgroundFlagsPro);
 // 4p:
 f4pBackgroundFlagsPro = new TProfile("f4pBackgroundFlagsPro","Flags and settings for 4p background histograms",1,0,1);
 f4pBackgroundFlagsPro->SetTickLength(-0.01,"Y");
 f4pBackgroundFlagsPro->SetMarkerStyle(25);
 f4pBackgroundFlagsPro->SetLabelSize(0.04);
 f4pBackgroundFlagsPro->SetLabelOffset(0.02,"Y");
 f4pBackgroundFlagsPro->SetStats(kFALSE);
 f4pBackgroundFlagsPro->SetFillColor(kGray);
 f4pBackgroundFlagsPro->SetLineColor(kBlack);
 f4pBackgroundFlagsPro->GetXaxis()->SetBinLabel(1,"TBI"); f4pBackgroundFlagsPro->Fill(0.5,0); // TBI dummy at the moment
 fBackgroundSublist[2]->Add(f4pBackgroundFlagsPro);

 //if(!fEstimateBackground){return;} // TBI is this safe?

 const Int_t nParticleSpecies = 5; // Supported at the moment: 0=e,1=mu,2=pi,3=K,4=p
 TString sParticles[2*nParticleSpecies] = {"e^{+}","#mu^{+}","#pi^{+}","K^{+}","p^{+}","e^{-}","#mu^{-}","#pi^{-}","K^{-}","p^{-}"};

 // b) Book all histograms for 2p background:
 if(fEstimate2pBackground)
 {
  // Remark 0: First particle in the pair is always the one with positive charge.
  // Remark 1: Diagonal elements are particle/antiparticle pairs.
  for(Int_t pid1=0;pid1<2*nParticleSpecies;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   for(Int_t pid2=pid1;pid2<2*nParticleSpecies;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
   {
    // Background:
    f2pBackground[pid1][pid2] = new TH1F(Form("f2pBackground[%d][%d]",pid1,pid2),Form("f2pBackground[%d][%d] = (%s,%s)",pid1,pid2,sParticles[pid1].Data(),sParticles[pid2].Data()),10000,0.,10.); // TBI rethink the boundaries and nbins
    f2pBackground[pid1][pid2]->SetStats(kTRUE);
    f2pBackground[pid1][pid2]->SetFillColor(kRed-10);
    f2pBackground[pid1][pid2]->SetLineColor(kRed);
    f2pBackground[pid1][pid2]->SetXTitle("k");
    f2pBackground[pid1][pid2]->SetYTitle("B(k)");
    fBackgroundSublist[0]->Add(f2pBackground[pid1][pid2]);
   } // for(Int_t pid2=0;pid2<5;pid2++) // anti-particle [0=e,1=mu,2=pi,3=K,4=p]
  } // for(Int_t pid=0;pid<5;pid++) // particle [0=e,1=mu,2=pi,3=K,4=p]
 } // if(fEstimate2pBackground)

 // c) Book all histograms for 3p background:
 if(fEstimate3pBackground)
 {
  for(Int_t pid1=0;pid1<2*nParticleSpecies;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   for(Int_t pid2=0;pid2<2*nParticleSpecies;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
   {
    for(Int_t pid3=0;pid3<2*nParticleSpecies;pid3++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
    {
     // Background:
     f3pBackground[pid1][pid2][pid3] = new TH1F(Form("f3pBackground[%d][%d][%d]",pid1,pid2,pid3),Form("f3pBackground[%d][%d][%d] = (%s,%s,%s)",pid1,pid2,pid3,sParticles[pid1].Data(),sParticles[pid2].Data(),sParticles[pid3].Data()),10000,0.,10.); // TBI rethink the boundaries and nbins
     f3pBackground[pid1][pid2][pid3]->SetStats(kTRUE);
     f3pBackground[pid1][pid2][pid3]->SetFillColor(kRed-10);
     f3pBackground[pid1][pid2][pid3]->SetLineColor(kRed);
     f3pBackground[pid1][pid2][pid3]->SetXTitle("Q_{3}");
     f3pBackground[pid1][pid2][pid3]->SetYTitle("B(Q_{3})");
     fBackgroundSublist[1]->Add(f3pBackground[pid1][pid2][pid3]);
    } // for(Int_t pid3=0;pid3<2*nParticleSpecies;pid3++) // anti-particle [0=e,1=mu,2=pi,3=K,4=p]
   } // for(Int_t pid2=0;pid2<2*nParticleSpecies;pid2++) // anti-particle [0=e,1=mu,2=pi,3=K,4=p]
  } // for(Int_t pid=0;pid<2*nParticleSpecies;pid++) // particle [0=e,1=mu,2=pi,3=K,4=p]
 } // if(fEstimate2pBackground)

 // d) Book all histograms for 4p background:
 if(fEstimate4pBackground)
 {
  for(Int_t pid1=0;pid1<2*nParticleSpecies;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   for(Int_t pid2=0;pid2<2*nParticleSpecies;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
   {
    for(Int_t pid3=0;pid3<2*nParticleSpecies;pid3++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
    {
     for(Int_t pid4=0;pid4<2*nParticleSpecies;pid4++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
     {
      // Background:
      f4pBackground[pid1][pid2][pid3][pid4] = new TH1F(Form("f4pBackground[%d][%d][%d][%d]",pid1,pid2,pid3,pid4),Form("f4pBackground[%d][%d][%d][%d] = (%s,%s,%s,%s)",pid1,pid2,pid3,pid4,sParticles[pid1].Data(),sParticles[pid2].Data(),sParticles[pid3].Data(),sParticles[pid4].Data()),10000,0.,10.); // TBI rethink the boundaries and nbins
      f4pBackground[pid1][pid2][pid3][pid4]->SetStats(kTRUE);
      f4pBackground[pid1][pid2][pid3][pid4]->SetFillColor(kRed-10);
      f4pBackground[pid1][pid2][pid3][pid4]->SetLineColor(kRed);
      f4pBackground[pid1][pid2][pid3][pid4]->SetXTitle("Q_{4}");
      f4pBackground[pid1][pid2][pid3][pid4]->SetYTitle("B(Q_{4})");
      fBackgroundSublist[2]->Add(f4pBackground[pid1][pid2][pid3][pid4]);
     } // for(Int_t pid4=0;pid4<2*nParticleSpecies;pid4++) // anti-particle [0=e,1=mu,2=pi,3=K,4=p]
    } // for(Int_t pid3=0;pid3<2*nParticleSpecies;pid3++) // anti-particle [0=e,1=mu,2=pi,3=K,4=p]
   } // for(Int_t pid2=0;pid2<2*nParticleSpecies;pid2++) // anti-particle [0=e,1=mu,2=pi,3=K,4=p]
  } // for(Int_t pid=0;pid<2*nParticleSpecies;pid++) // particle [0=e,1=mu,2=pi,3=K,4=p]
 } // if(fEstimate2pBackground)

 // e) Book buffer objects:
 switch(fBackgroundOption)
 {
  case 0: // shifting (see documentation in EstimateBackground(...))
   for(Int_t me=0;me<3;me++) // [0] is buffer for 1st event; [1] for 2nd, etc.
   {
    fMixedEvents0[me] = new TClonesArray("AliAODTrack",10000);
   }
  break; // case 0: // shifting

  case 1: // permutations + binning in vertex z-position
   for(Int_t vzr=0;vzr<10;vzr++) // vertez-z range
   {
    for(Int_t me=0;me<fMaxBufferSize1;me++) // buffer for mixed events
    {
     fMixedEvents1[vzr][me] = new TClonesArray("AliAODTrack",10000);
     fGlobalTracksAOD1[vzr][me] = new TExMap();
    }
   }
  break; // permutations + binning in vertex z-position

  /*
  default:
   cout<<Form("And the fatal 'fBackgroundOption' value is... %d. Congratulations!!",fBackgroundOption)<<endl;
   Fatal(sMethodName.Data(),"switch(fBackgroundOption)");
  */
 } // switch(fBackgroundOption)

} // void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForBackground()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForBuffers()
{
 // Book all the stuff for buffers.

 // a) Book the profile holding all the flags for buffers;
 // b) Book TClonesArray *fChargedParticlesCA[2][10][10000];
 // c) Book TExMap *fChargedParticlesEM[10];

 // a) Book the profile holding all the flags for buffers:
 fBuffersFlagsPro = new TProfile("fBuffersFlagsPro","Flags and settings for buffers",2,0,2);
 fBuffersFlagsPro->SetTickLength(-0.01,"Y");
 fBuffersFlagsPro->SetMarkerStyle(25);
 fBuffersFlagsPro->SetLabelSize(0.04);
 fBuffersFlagsPro->SetLabelOffset(0.02,"Y");
 fBuffersFlagsPro->SetStats(kFALSE);
 fBuffersFlagsPro->SetFillColor(kGray);
 fBuffersFlagsPro->SetLineColor(kBlack);
 fBuffersFlagsPro->GetXaxis()->SetBinLabel(1,"fFillBuffers"); fBuffersFlagsPro->Fill(0.5,(Int_t)fFillBuffers);
 fBuffersFlagsPro->GetXaxis()->SetBinLabel(2,"fMaxBuffer"); fBuffersFlagsPro->Fill(1.5,fMaxBuffer);
 fBuffersList->Add(fBuffersFlagsPro);

 if(!fFillBuffers){return;} // TBI is this safe? Well, let's hope so for the time being...

 // b) Book TClonesArray *fChargedParticlesCA[2][10][10000]:
 for(Int_t am=0;am<1+(Int_t)fProcessBothKineAndReco;am++) // analysis method [0=AOD||ESD,1=MC]
 {
  for(Int_t e=0;e<fMaxBuffer;e++) // number of events to buffer, max = apparently 10
  {
   for(Int_t p=0;p<10000;p++) // charged particles, max = apparently 10000
   {
    fChargedParticlesCA[am][e][p] = NULL; // TBI not sure if I need this really here, this way...
    // Example: fChargedParticlesCA[0][4][6] is 6th charged particles in 5th buffered event, in AOD||ESD analysis
   } // for(Int_t p=0;p<10000;p++) // charged particles, max = apparently 10000
  } // for(Int_t e=0;e<10;e++) // number of events to buffer, max = apparently 10
 } // for(Int_t am=0;am<2;am++) // analysis method [0=AOD||ESD,1=MC]

 // c) Book TExMap *fChargedParticlesEM[10]:
 for(Int_t e=0;e<fMaxBuffer;e++) // number of events to buffer, max = apparently 10
 {
  fChargedParticlesEM[e] = new TExMap();
 } // for(Int_t e=0;e<fMaxBuffer;e++) // number of events to buffer, max = apparently 10

} // void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForBuffers()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForQA()
{
 // Book all the stuff for QA.

 // a) Book the profile holding all the flags for QA;
 // b) Book all objects for "QA events";
 // c) Book all objects for "QA particles".

 // a) Book the profile holding all the flags for QA:
 fQAFlagsPro = new TProfile("fQAFlagsPro","Flags and settings for QA",4,0.,4.);
 fQAFlagsPro->SetTickLength(-0.01,"Y");
 fQAFlagsPro->SetMarkerStyle(25);
 fQAFlagsPro->SetLabelSize(0.04);
 fQAFlagsPro->SetLabelOffset(0.02,"Y");
 fQAFlagsPro->SetStats(kFALSE);
 fQAFlagsPro->SetFillColor(kGray);
 fQAFlagsPro->SetLineColor(kBlack);
 fQAFlagsPro->GetXaxis()->SetBinLabel(1,"fFillQA"); fQAFlagsPro->Fill(0.5,(Int_t)fFillQA);
 fQAFlagsPro->GetXaxis()->SetBinLabel(2,"fBailOutAfterQA"); fQAFlagsPro->Fill(1.5,(Int_t)fBailOutAfterQA);
 fQAFlagsPro->GetXaxis()->SetBinLabel(3,"fFillQAEvents"); fQAFlagsPro->Fill(2.5,(Int_t)fFillQAEvents);
 fQAFlagsPro->GetXaxis()->SetBinLabel(4,"fFillQAParticles"); fQAFlagsPro->Fill(3.5,(Int_t)fFillQAParticles);
 fQAList->Add(fQAFlagsPro);

 if(!fFillQA){return;}

 // b) Book all objects for "QA events":
 if(fFillQAEvents)
 {
  // ...
 } // if(fFillQAEvents)

 // c) Book all objects for "QA particles":
 if(fFillQAParticles)
 {
  // c0) Book fQAFilterBitScan:
  Int_t nFilterBits = 14;
  fQAFilterBitScan = new TH1I("fQAFilterBitScan","fQAFilterBitScan",nFilterBits,0,nFilterBits);
  fQAFilterBitScan->SetStats(kFALSE);
  fQAFilterBitScan->SetFillColor(kBlue-10);
  fQAFilterBitScan->SetXTitle("FilterBit");
  fQAFilterBitScan->SetYTitle("# particles");
  for(Int_t fb=0;fb<nFilterBits;fb++) // 'fb' is a 'left shifting opetator', i.e. 1<<fb = filterbit
  {
   fQAFilterBitScan->GetXaxis()->SetBinLabel(fb+1,Form("%d",1<<fb));
  }
  fQAParticlesList->Add(fQAFilterBitScan);

  // c1) Book fQAIDvsFilterBit:
  Int_t nIDsMax = 10000;
  fQAIDvsFilterBit = new TH2I("fQAIDvsFilterBit","fQAIDvsFilterBit",nFilterBits,0,nFilterBits,2*nIDsMax,-nIDsMax,nIDsMax);
  //fQAIDvsFilterBit->SetStats(kFALSE);
  //fQAIDvsFilterBit->SetFillColor(kBlue-10);
  fQAIDvsFilterBit->SetXTitle("FilterBit");
  fQAIDvsFilterBit->SetYTitle("atrack->GetID()");
  for(Int_t fb=0;fb<nFilterBits;fb++) // 'fb' is a 'left shifting opetator', i.e. 1<<fb = filterbit
  {
   fQAIDvsFilterBit->GetXaxis()->SetBinLabel(fb+1,Form("%d",1<<fb));
  }
  fQAParticlesList->Add(fQAIDvsFilterBit);

  // c2) Book fQAParticleHist[2][10][10];
  TString sBeforeAfter[2] = {"before","after"};
  TString sDistributions[10] = {"p_{T}","#eta","","","","","","","",""}; // Only the ones with non-empty entries are booked and filled, see below if(sDistributions[di].EqualTo("")){continue;}
  TString sCuts[10] = {"PID","","","","","","","","",""}; // Only the ones with non-empty entries are booked and filled, see below if(sCuts[ci].EqualTo("")){continue;}
  for(Int_t ba=0;ba<2;ba++) // 0="before rain",1="after rain"
  {
   for(Int_t di=0;di<10;di++) // number of distinct distributions, maximum is 10 (at the moment, as it seems...)
   {
    if(sDistributions[di].EqualTo("")){continue;}
    for(Int_t ci=0;ci<10;ci++) // number of distinct track cuts, maximum is 10 (at the moment, as it seems...)
    {
     if(sCuts[ci].EqualTo("")){continue;}
      fQAParticleHist[ba][di][ci] = new TH1F(Form("fQAParticleHist[%d][%d][%d]",ba,di,ci),Form("fQAParticleHist[%d][%d][%d] = (%s,%s,%s)",ba,di,ci,sBeforeAfter[ba].Data(),sDistributions[di].Data(),sCuts[ci].Data() ),10000,0.,10.); // TBI rethink the boundaries and nbins
      fQAParticleHist[ba][di][ci]->SetStats(kFALSE);
      fQAParticleHist[ba][di][ci]->SetFillColor(kBlue-10);
      fQAParticleHist[ba][di][ci]->SetXTitle("TBI");
      fQAParticleHist[ba][di][ci]->SetYTitle("TBI");
      fQAParticlesList->Add(fQAParticleHist[ba][di][ci]);
    } // for(Int_t ci=0;ci<10;ci++) // number of distinct track cuts, maximum is 10 (at the moment)
   } // for(Int_t di=0;di<10;di++) // number of distinct distributions, maximum is 10 (at the moment)
  } // for(Int_t ba=0;ba<2;ba++) // 0="before rain",1="after rain"
 } // if(fFillQAParticles)

} // void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForQA()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForGlobalTrackCuts()
{
 // Book all objects for global track cuts (applied only on "normal" global tracks in AOD).

 // a) Book the profile holding all the flags;
 // b) ...

 // a) Book the profile holding all the flags:
 fGlobalTrackCutsFlagsPro = new TProfile("fGlobalTrackCutsFlagsPro","Flags and settings for global track cuts",7,0.,7.);
 fGlobalTrackCutsFlagsPro->SetTickLength(-0.01,"Y");
 fGlobalTrackCutsFlagsPro->SetMarkerStyle(25);
 fGlobalTrackCutsFlagsPro->SetLabelSize(0.04);
 fGlobalTrackCutsFlagsPro->SetLabelOffset(0.02,"Y");
 fGlobalTrackCutsFlagsPro->SetStats(kFALSE);
 fGlobalTrackCutsFlagsPro->SetFillColor(kGray);
 fGlobalTrackCutsFlagsPro->SetLineColor(kBlack);
 fGlobalTrackCutsFlagsPro->GetXaxis()->SetBinLabel(1,"fApplyGlobalTrackCuts"); fGlobalTrackCutsFlagsPro->Fill(0.5,(Int_t)fApplyGlobalTrackCuts);
 fGlobalTrackCutsFlagsPro->GetXaxis()->SetBinLabel(2,"fPtRange[0]"); fGlobalTrackCutsFlagsPro->Fill(1.5,fPtRange[0]);
 fGlobalTrackCutsFlagsPro->GetXaxis()->SetBinLabel(3,"fPtRange[1]"); fGlobalTrackCutsFlagsPro->Fill(2.5,fPtRange[1]);
 fGlobalTrackCutsFlagsPro->GetXaxis()->SetBinLabel(4,"fEtaRange[0]"); fGlobalTrackCutsFlagsPro->Fill(3.5,fEtaRange[0]);
 fGlobalTrackCutsFlagsPro->GetXaxis()->SetBinLabel(5,"fEtaRange[1]"); fGlobalTrackCutsFlagsPro->Fill(4.5,fEtaRange[1]);
 fGlobalTrackCutsFlagsPro->GetXaxis()->SetBinLabel(6,"fPhiRange[0]"); fGlobalTrackCutsFlagsPro->Fill(5.5,fPhiRange[0]);
 fGlobalTrackCutsFlagsPro->GetXaxis()->SetBinLabel(7,"fPhiRange[1]"); fGlobalTrackCutsFlagsPro->Fill(6.5,fPhiRange[1]);
 fGlobalTrackCutsList->Add(fGlobalTrackCutsFlagsPro);

} // void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForGlobalTrackCuts()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::GlobalTracksAOD(AliAODEvent *aAOD, Int_t index)
{
 // Filter out unique global tracks in AOD and store them in fGlobalTracksAOD[index].

 // Remark 0: All global tracks have positive ID, the duplicated TPC-only tracks have -(ID+1);
 // Remark 1: The issue here is that there are apparently two sets of global tracks: a) "normal" and b) constrained to primary vertex.
 //           However, only the "normal" global tracks come with positive ID, additionally they can be discriminated simply via: aodTrack->IsGlobalConstrained()
 //           Global constrained tracks have the same negative ID as the TPC-only tracks, both associated with the same "normal global" tracks. E.g. we can have
 //            iTrack: atrack->GetID(): atrack->Pt() atrack->Eta() atrack->Phi()
 //                 1:               0:     2.073798     -0.503640      2.935432
 //                19:              -1:     2.075537     -0.495988      2.935377 => this is TPC-only
 //                35:              -1:     2.073740     -0.493576      2.935515 => this is IsGlobalConstrained()
 //           In fact, this is important, otherwise there is double or even triple counting in some cases.
 // Remark 2: There are tracks for which: 0 == aodTrack->GetFilterMap()
 //           a) Basically all of them pass: atrack->GetType() == AliAODTrack::kFromDecayVtx , but few exceptions also pass atrack->GetType() == AliAODTrack::kPrimary
 //           b) All of them apparently have positive ID, i.e. these are global tracks
 //           c) Clearly, we cannot use TestFilterBit() on them
 //           d) None of them apparently satisfies: atrack->IsGlobalConstrained()
 // Remark 3: There is a performance penalty when fGlobalTracksAOD[1] and fGlobalTracksAOD[2] needed for mixed events are calculated.
 //           Yes, I can get them directly from fGlobalTracksAOD[0], without calling this method for them again. TBI today

 // a) Insanity checks;
 // b) Determine the map.

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::GlobalTracksAOD(AliAODEvent *aAOD, Int_t index)";
 if(!fGlobalTracksAOD[index]){Fatal(sMethodName.Data(),"fGlobalTracksAOD[%d]",index);}
 if(0 != fGlobalTracksAOD[index]->GetSize()){fGlobalTracksAOD[index]->Delete();} // yes, this method determines mapping from scratch each time

 // b) Determine the map:
 for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++)
 {
  AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack));
  if(aodTrack)
  {
   Int_t id = aodTrack->GetID();
   //if(id>=0 && aodTrack->GetFilterMap()>0 && !aodTrack->IsGlobalConstrained()) // TBI rethink this
   if(id>=0 && !aodTrack->IsGlobalConstrained()) // TBI rethink this, it seems that id>=0 is just enough, the second constraint is most likely just an overkill
   {
    fGlobalTracksAOD[index]->Add(id,iTrack); // "key" = id, "value" = iTrack
   } // if(id>=0 && !aodTrack->IsGlobalConstrained())
  } // if(aodTrack)
 } // for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++)

} // void AliAnalysisTaskMultiparticleFemtoscopy::GlobalTracksAOD(AliAODEvent *aAOD)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::GlobalTracksAODTEST(AliAODEvent *aAOD, Int_t index)
{
 // Filter out unique global tracks in AOD and store them in fGlobalTracksAODTEST[index].

 // Remark 0: All global tracks have positive ID, the duplicated TPC-only tracks have -(ID+1);
 // Remark 1: The issue here is that there are apparently two sets of global tracks: a) "normal" and b) constrained to primary vertex.
 //           However, only the "normal" global tracks come with positive ID, additionally they can be discriminated simply via: aodTrack->IsGlobalConstrained()
 //           Global constrained tracks have the same negative ID as the TPC-only tracks, both associated with the same "normal global" tracks. E.g. we can have
 //            iTrack: atrack->GetID(): atrack->Pt() atrack->Eta() atrack->Phi()
 //                 1:               0:     2.073798     -0.503640      2.935432
 //                19:              -1:     2.075537     -0.495988      2.935377 => this is TPC-only
 //                35:              -1:     2.073740     -0.493576      2.935515 => this is IsGlobalConstrained()
 //           In fact, this is important, otherwise there is double or even triple counting in some cases.
 // Remark 2: There are tracks for which: 0 == aodTrack->GetFilterMap()
 //           a) Basically all of them pass: atrack->GetType() == AliAODTrack::kFromDecayVtx , but few exceptions also pass atrack->GetType() == AliAODTrack::kPrimary
 //           b) All of them apparently have positive ID, i.e. these are global tracks
 //           c) Clearly, we cannot use TestFilterBit() on them
 //           d) None of them apparently satisfies: atrack->IsGlobalConstrained()
 // Remark 3: There is a performance penalty when fGlobalTracksAODTEST[1] and fGlobalTracksAODTEST[2] needed for mixed events are calculated.
 //           Yes, I can get them directly from fGlobalTracksAODTEST[0], without calling this method for them again. TBI today

 // a) Insanity checks;
 // b) Determine the map.

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::GlobalTracksAODTEST(AliAODEvent *aAOD, Int_t index)";
 if(!fGlobalTracksAODTEST[index]){Fatal(sMethodName.Data(),"fGlobalTracksAODTEST[%d]",index);}
 if(0 != fGlobalTracksAODTEST[index]->GetSize()){fGlobalTracksAODTEST[index]->Delete();} // yes, this method determines mapping from scratch each time

 // b) Determine the map:
 for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++)
 {
  AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack));
  if(aodTrack)
  {
   Int_t id = aodTrack->GetID();
   //if(id>=0 && aodTrack->GetFilterMap()>0 && !aodTrack->IsGlobalConstrained()) // TBI rethink this
   if(id>=0 && !aodTrack->IsGlobalConstrained()) // TBI rethink this, it seems that id>=0 is just enough, the second constraint is most likely just an overkill
   {
    fGlobalTracksAODTEST[index]->Add(id,iTrack); // "key" = id, "value" = iTrack
   } // if(id>=0 && !aodTrack->IsGlobalConstrained())
  } // if(aodTrack)
 } // for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++)

} // void AliAnalysisTaskMultiparticleFemtoscopy::GlobalTracksAODTEST(AliAODEvent *aAOD)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::GlobalTracksAODHA(AliAODEvent *aAOD, Int_t index)
{
 // Filter out unique global tracks in AOD and store them in fGlobalTracksAODHA[index].

 // Remark 0: All global tracks have positive ID, the duplicated TPC-only tracks have -(ID+1);
 // Remark 1: The issue here is that there are apparently two sets of global tracks: a) "normal" and b) constrained to primary vertex.
 //           However, only the "normal" global tracks come with positive ID, additionally they can be discriminated simply via: aodTrack->IsGlobalConstrained()
 //           Global constrained tracks have the same negative ID as the TPC-only tracks, both associated with the same "normal global" tracks. E.g. we can have
 //            iTrack: atrack->GetID(): atrack->Pt() atrack->Eta() atrack->Phi()
 //                 1:               0:     2.073798     -0.503640      2.935432
 //                19:              -1:     2.075537     -0.495988      2.935377 => this is TPC-only
 //                35:              -1:     2.073740     -0.493576      2.935515 => this is IsGlobalConstrained()
 //           In fact, this is important, otherwise there is double or even triple counting in some cases.
 // Remark 2: There are tracks for which: 0 == aodTrack->GetFilterMap()
 //           a) Basically all of them pass: atrack->GetType() == AliAODTrack::kFromDecayVtx , but few exceptions also pass atrack->GetType() == AliAODTrack::kPrimary
 //           b) All of them apparently have positive ID, i.e. these are global tracks
 //           c) Clearly, we cannot use TestFilterBit() on them
 //           d) None of them apparently satisfies: atrack->IsGlobalConstrained()
 // Remark 3: There is a performance penalty when fGlobalTracksAODHA[1] and fGlobalTracksAODHA[2] needed for mixed events are calculated.
 //           Yes, I can get them directly from fGlobalTracksAODHA[0], without calling this method for them again. TBI today

 // a) Insanity checks;
 // b) Determine the map.

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::GlobalTracksAODHA(AliAODEvent *aAOD, Int_t index)";
 if(!fGlobalTracksAODHA[index]){Fatal(sMethodName.Data(),"fGlobalTracksAODHA[%d]",index);}
 if(0 != fGlobalTracksAODHA[index]->GetSize()){fGlobalTracksAODHA[index]->Delete();} // yes, this method determines mapping from scratch each time

 // b) Determine the map:
 for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++)
 {
  AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack));
  if(aodTrack)
  {
   Int_t id = aodTrack->GetID();
   //if(id>=0 && aodTrack->GetFilterMap()>0 && !aodTrack->IsGlobalConstrained()) // TBI rethink this
   if(id>=0 && !aodTrack->IsGlobalConstrained()) // TBI rethink this, it seems that id>=0 is just enough, the second constraint is most likely just an overkill
   {
    fGlobalTracksAODHA[index]->Add(id,iTrack); // "key" = id, "value" = iTrack
   } // if(id>=0 && !aodTrack->IsGlobalConstrained())
  } // if(aodTrack)
 } // for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++)

} // void AliAnalysisTaskMultiparticleFemtoscopy::GlobalTracksAODHA(AliAODEvent *aAOD, Int_t index)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::GlobalTracksAOD(AliAODEvent *aAOD, Int_t indexX, Int_t indexY)
{
 // Filter out unique global tracks in AOD and store them in fGlobalTracksAOD1[10][5]. TBI hardwired "1" in an array, not in the function name, so this is a bit incosistent...

 // Remark 0: All global tracks have positive ID, the duplicated TPC-only tracks have -(ID+1);
 // Remark 1: The issue here is that there are apparently two sets of global tracks: a) "normal" and b) constrained to primary vertex.
 //           However, only the "normal" global tracks come with positive ID, additionally they can be discriminated simply via: aodTrack->IsGlobalConstrained()
 //           Global constrained tracks have the same negative ID as the TPC-only tracks, both associated with the same "normal global" tracks. E.g. we can have
 //            iTrack: atrack->GetID(): atrack->Pt() atrack->Eta() atrack->Phi()
 //                 1:               0:     2.073798     -0.503640      2.935432
 //                19:              -1:     2.075537     -0.495988      2.935377 => this is TPC-only
 //                35:              -1:     2.073740     -0.493576      2.935515 => this is IsGlobalConstrained()
 //           In fact, this is important, otherwise there is double or even triple counting in some cases.
 // Remark 2: There are tracks for which: 0 == aodTrack->GetFilterMap()
 //           a) Basically all of them pass: atrack->GetType() == AliAODTrack::kFromDecayVtx , but few exceptions also pass atrack->GetType() == AliAODTrack::kPrimary
 //           b) All of them apparently have positive ID, i.e. these are global tracks
 //           c) Clearly, we cannot use TestFilterBit() on them
 //           d) None of them apparently satisfies: atrack->IsGlobalConstrained()
 // Remark 3: Used only for fBackgroundOption = 1, i.e. permutations + binning in vertex z-position

 // a) Insanity checks;
 // b) Determine the map.

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::GlobalTracksAOD(AliAODEvent *aAOD, Int_t indexX, Int_t indexY)";
 if(!fGlobalTracksAOD1[indexX][indexY]){Fatal(sMethodName.Data(),"fGlobalTracksAOD[%d][%d]",indexX,indexY);}
 if(0 != fGlobalTracksAOD1[indexX][indexY]->GetSize()){fGlobalTracksAOD1[indexX][indexY]->Delete();} // yes, this method determines mapping from scratch each time

 // b) Determine the map:
 for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++)
 {
  AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack));
  if(aodTrack)
  {
   Int_t id = aodTrack->GetID();
   //if(id>=0 && aodTrack->GetFilterMap()>0 && !aodTrack->IsGlobalConstrained()) // TBI rethink this
   if(id>=0 && !aodTrack->IsGlobalConstrained()) // TBI rethink this, it seems that id>=0 is just enough, the second constraint is most likely just an overkill
   {
    fGlobalTracksAOD1[indexX][indexY]->Add(id,iTrack); // "key" = id, "value" = iTrack
   } // if(id>=0 && !aodTrack->IsGlobalConstrained())
  } // if(aodTrack)
 } // for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++)

} // void AliAnalysisTaskMultiparticleFemtoscopy::GlobalTracksAOD(AliAODEvent *aAOD, Int_t indexX, Int_t indexY)

//=======================================================================================================================

Bool_t AliAnalysisTaskMultiparticleFemtoscopy::SpecifiedEvent(UInt_t run, UShort_t bunchCross, UInt_t orbit, UInt_t period)
{
 // Check if this is event specified in a steering macro via the setter void SetWaitForSpecifiedEvent(UInt_t run, UShort_t bunchCross, UInt_t orbit, UInt_t period).

 if(run != fRun) return kFALSE;
 else if(bunchCross != fBunchCross) return kFALSE;
 else if(orbit != fOrbit) return kFALSE;
 else if(period != fPeriod) return kFALSE;

 return kTRUE;

} // void AliAnalysisTaskMultiparticleFemtoscopy::SpecifiedEvent(UInt_t run, UShort_t bunchCross, UInt_t orbit, UInt_t period)

//=======================================================================================================================

Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesGlobalTrackCuts(AliAODTrack *gtrack)
{
 // Check if the track passes common global track cuts (irrespectively of PID).

 TString sMethodName = "Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesGlobalTrackCuts(AliAODTrack *gtrack)";
 if(!gtrack){Fatal(sMethodName.Data(),"!gtrack");}

 // To do: add data members and corresponding setters:
 // fTPCNclsMin, fTPCNclsMax
 // fTPCsignalNMin, fTPCsignalNMax
 if(gtrack->Pt()<fPtRange[0]) return kFALSE;
 if(gtrack->Pt()>=fPtRange[1]) return kFALSE;
 if(gtrack->Eta()<fEtaRange[0]) return kFALSE;
 if(gtrack->Eta()>=fEtaRange[1]) return kFALSE;
 if(gtrack->Phi()<fPhiRange[0]) return kFALSE;
 if(gtrack->Phi()>=fPhiRange[1]) return kFALSE;
 //if(gtrack->GetTPCNcls()<70) return kFALSE;
 //if(gtrack->GetTPCsignalN()<70) return kFALSE;

 if(fRejectFakeTracks && gtrack->GetLabel()<0) return kFALSE; // TBI is it more precise <=0 ?

 return kTRUE;

} // Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesGlobalTrackCuts(AliAODTrack *gtrack)

//=======================================================================================================================

Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesCommonTrackCuts(AliAODTrack *atrack)
{
 // Check if the track passes common analysis specific track (e.g. TPC-only) cuts.
 // Therefore we can cut independetly on global track parameters, and on TPC-only cut parameters.

 TString sMethodName = "Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesCommonTrackCuts(AliAODTrack *atrack)";
 if(!atrack){Fatal(sMethodName.Data(),"!atrack");}

 if(!atrack->TestFilterBit(128)) return kFALSE; // TPC-only TBI setter TBI#2 There might be some conflict with FillControlHistogramsNonIdentifiedParticlesFTSF

 // TBI: implement eventually separate cuts for 'atracks' and 'gtracks'
 if(atrack->Pt()<fPtRange[0]) return kFALSE;
 if(atrack->Pt()>=fPtRange[1]) return kFALSE;
 if(atrack->Eta()<fEtaRange[0]) return kFALSE;
 if(atrack->Eta()>=fEtaRange[1]) return kFALSE;
 if(atrack->Phi()<fPhiRange[0]) return kFALSE;
 if(atrack->Phi()>=fPhiRange[1]) return kFALSE;
 //if(gtrack->GetTPCNcls()<70) return kFALSE;
 //if(gtrack->GetTPCsignalN()<70) return kFALSE;

 if(fRejectFakeTracks && atrack->GetLabel()<0) return kFALSE; // TBI is it more precise <=0 ?

 return kTRUE;

} // Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesCommonTrackCuts(AliAODTrack *atrack)

//=======================================================================================================================

Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesCommonTrackCuts(AliAODMCParticle *amcparticle)
{
 // TBI not validated. this method applies only to MC, make it uniform with AOD.

 TString sMethodName = "Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesCommonTrackCuts(AliAODMCParticle *amcparticle)";
 if(!amcparticle){Fatal(sMethodName.Data(),"!amcparticle");}

 if(TMath::Abs(amcparticle->Charge())<1.e-10) return kFALSE; // skipping neutral particles TBI
 if(!amcparticle->IsPhysicalPrimary()) return kFALSE; // skipping secondaries TBI
 if(amcparticle->Pt()<0.2) return kFALSE;
 if(amcparticle->Pt()>=5.0) return kFALSE;
 if(amcparticle->Eta()<-0.8) return kFALSE;
 if(amcparticle->Eta()>=0.8) return kFALSE;

 return kTRUE;

} // Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesCommonTrackCuts(AliAODMCParticle *amcparticle)

//=======================================================================================================================

Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesCommonEventCuts(AliVEvent *ave)
{
 // Check if the event passes common event cuts.

 // TBI: add support for fProcessBothKineAndReco

 // a) Determine Ali{MC,ESD,AOD}Event;

 // a) Determine Ali{MC,ESD,AOD}Event:
 AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 if(aMC && fProcessOnlyKine)
 {
  // TBI
 }
 else if(aESD)
 {
  // TBI
 }
 else if(aAOD) // shall work both for fProcessBothKineAndReco = kTRUE and fProcessOnlyReco = kTRUE
 {
  // a) Cuts on AliAODEvent:
  if(fRejectEventsWithoutPrimaryVertex && !aAOD->GetPrimaryVertex()) return kFALSE;
  if(TMath::Abs(aAOD->GetMagneticField())<fMinMagneticField) return kFALSE;
  if(fCutOnNumberOfTracks)
  {
   if(aAOD->GetNumberOfTracks() < fMinNumberOfTracks) return kFALSE;
   if(aAOD->GetNumberOfTracks() > fMaxNumberOfTracks) return kFALSE;
  }
  if(fCutOnNumberOfV0s)
  {
   if(aAOD->GetNumberOfV0s() < fMinNumberOfV0s) return kFALSE;
   if(aAOD->GetNumberOfV0s() > fMaxNumberOfV0s) return kFALSE;
  }
  if(fCutOnNumberOfCascades)
  {
   if(aAOD->GetNumberOfCascades() < fMinNumberOfCascades) return kFALSE;
   if(aAOD->GetNumberOfCascades() > fMaxNumberOfCascades) return kFALSE;
  }

  // b) Cuts on AliAODVertex:
  AliAODVertex *avtx = (AliAODVertex*)aAOD->GetPrimaryVertex();
  if(fCutOnVertexX)
  {
   if(avtx->GetX() < fMinVertexX) return kFALSE;
   if(avtx->GetX() > fMaxVertexX) return kFALSE;
  }
  if(fCutOnVertexY)
  {
   if(avtx->GetY() < fMinVertexY) return kFALSE;
   if(avtx->GetY() > fMaxVertexY) return kFALSE;
  }
  if(fCutOnVertexZ)
  {
   if(avtx->GetZ() < fMinVertexZ) return kFALSE;
   if(avtx->GetZ() > fMaxVertexZ) return kFALSE;
  }
  if(fCutOnNContributors)
  {
   if(avtx->GetNContributors() < fMinNContributors) return kFALSE;
   if(avtx->GetNContributors() > fMaxNContributors) return kFALSE;
  }
 } // else if(aAOD)

 return kTRUE;

} // Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesCommonEventCuts(AliVEvent *ave)

//=======================================================================================================================

Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesMixedEventCuts(AliVEvent *ave)
{
 // Check if the event passes mixed event cuts.

 // Remark: PassesCommonEventCuts() is already checked beforehand in UserExec(), therefore here only some new cuts characteristic only for mixed events shall be implemented.

 // a) Determine Ali{MC,ESD,AOD}Event;

 // a) Determine Ali{MC,ESD,AOD}Event:
 AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 if(aMC)
 {
  // TBI
 }
 else if(aESD)
 {
  // TBI
 }
 else if(aAOD)
 {

  // TBI
  /* Check whether the cuts below are already being applied in PassesCommonEventCuts()
  if(!aAOD->GetPrimaryVertex()) return kFALSE;
  if(TMath::Abs(aAOD->GetMagneticField())<0.001) return kFALSE;
  AliAODVertex *avtx = (AliAODVertex*)aAOD->GetPrimaryVertex();
  //if(TMath::Abs(avtx->GetX())>10.0) return kFALSE;
  //if(TMath::Abs(avtx->GetY())>10.0) return kFALSE;
  if(TMath::Abs(avtx->GetZ())>10.0) return kFALSE; // TBI setter
  if(avtx->GetNContributors()<=2) return kFALSE; // TBI setter
  */
 }

 return kTRUE;

} // Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesMixedEventCuts(AliVEvent *ave)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::CalculateCorrelationFunctions(AliAODEvent *aAOD)
{
 // Calculate correlation functions.

 // a) Insanity checks;
 // b) Two nested loops to calculate C(k), just an example; TBI
 // c) Three nested loops to calculate C(Q3), just an example; TBI
 // d) Four nested loops to calculate C(Q4), just an example. TBI

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::CalculateCorrelationFunctions(AliAODEvent *aAOD)";
 if(!aAOD){Fatal(sMethodName.Data(),"!aAOD");}
 if(0 == fGlobalTracksAOD[0]->GetSize()){Fatal(sMethodName.Data(),"0 == fGlobalTracksAOD[0]->GetSize()");} // this case shall be already treated in UserExec

 // b) Two nested loops to calculate C(k), just an example:
 Int_t nTracks = aAOD->GetNumberOfTracks();
 for(Int_t iTrack1=0;iTrack1<nTracks;iTrack1++)
 {
  AliAODTrack *atrack1 = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack1));
  // TBI Temporary track insanity checks:
  if(!atrack1){Fatal(sMethodName.Data(),"!atrack1");} // TBI keep this for some time, eventually just continue
  if(atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(atrack1->TestFilterBit(128) && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->TestFiletrBit(128) && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(!PassesCommonTrackCuts(atrack1)){continue;} // TBI re-think
  // Corresponding AOD global track:
  Int_t id1 = atrack1->GetID();
  AliAODTrack *gtrack1 = dynamic_cast<AliAODTrack*>(id1>=0 ? aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(id1)) : aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(-(id1+1))));
  if(!gtrack1){Fatal(sMethodName.Data(),"!gtrack1");} // TBI keep this for some time, eventually just continue
  Int_t gid1 = (id1>=0 ? id1 : -(id1+1)); // ID of corresponding global track
  // Common track selection criteria for all "normal" global tracks:
  if(!PassesGlobalTrackCuts(gtrack1)){continue;}

  // Loop over the 2nd particle:
  for(Int_t iTrack2=iTrack1+1;iTrack2<nTracks;iTrack2++)
  {
   AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack2));
   // TBI Temporary track insanity checks:
   if(!atrack2){Fatal(sMethodName.Data(),"!atrack2");} // TBI keep this for some time, eventually just continue
   if(atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(atrack2->TestFilterBit(128) && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->TestFiletrBit(128) && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(!PassesCommonTrackCuts(atrack2)){continue;} // TBI re-think
   // Corresponding AOD global track:
   Int_t id2 = atrack2->GetID();
   AliAODTrack *gtrack2 = dynamic_cast<AliAODTrack*>(id2>=0 ? aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(id2)) : aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(-(id2+1))));
   if(!gtrack2){Fatal(sMethodName.Data(),"!gtrack2");} // TBI keep this for some time, eventually just continue
   Int_t gid2 = (id2>=0 ? id2 : -(id2+1)); // ID of corresponding global track
   if(gid1==gid2){continue;} // Eliminate self-correlations:

   // Common track selection criteria for all "normal" global tracks:
   if(!PassesGlobalTrackCuts(gtrack2)){continue;}

   // Okay, so we have two tracks, let's decide whether kinemetics from 'atrack' or 'gtrack' is taken, let's check PID, and fill the correlation functions:
   AliAODTrack *agtrack1 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
   AliAODTrack *agtrack2 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
   if(fFillControlHistogramsWithGlobalTrackInfo)
   {
    agtrack1 = gtrack1;
    agtrack2 = gtrack2;
   }
   else
   {
    agtrack1 = atrack1;
    agtrack2 = atrack2;
   } // if(fFillControlHistogramsWithGlobalTrackInfo)
   if(!agtrack1){Fatal(sMethodName.Data(),"!agtrack1");}
   if(!agtrack2){Fatal(sMethodName.Data(),"!agtrack2");}

   // 1.) Same particle species:

   // a) pion-pion:
   //  a1) pi+pi+ [2][2]:
   if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE))
   {
    fCorrelationFunctions[2][2]->Fill(Q2(agtrack1,agtrack2));
   }
   //  a2) pi-pi- [7][7]:
   if(Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[7][7]->Fill(Q2(agtrack1,agtrack2));
   }
   //  a3) pi+pi- || pi-pi+ [2][7]:
   if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,1,kTRUE)))
   {
    fCorrelationFunctions[2][7]->Fill(Q2(agtrack1,agtrack2));
   }

   // b) kaon-kaon:
   //  b1) K+K+ [3][3]:
   if(Kaon(gtrack1,1,kTRUE) && Kaon(gtrack2,1,kTRUE))
   {
    fCorrelationFunctions[3][3]->Fill(Q2(agtrack1,agtrack2));
   }
   //  b2) K-K- [8][8]:
   if(Kaon(gtrack1,-1,kTRUE) && Kaon(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[8][8]->Fill(Q2(agtrack1,agtrack2));
   }
   //  b3) K+K- || K-K+ [3][8]:
   if((Kaon(gtrack1,1,kTRUE) && Kaon(gtrack2,-1,kTRUE)) || (Kaon(gtrack1,-1,kTRUE) && Kaon(gtrack2,1,kTRUE)))
   {
    fCorrelationFunctions[3][8]->Fill(Q2(agtrack1,agtrack2));
   }

   // c) proton-proton:
   //  c1) p+p+ [4][4]:
   if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE))
   {
    fCorrelationFunctions[4][4]->Fill(Q2(agtrack1,agtrack2));
   }
   //  c2) p-p- [9][9]:
   if(Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[9][9]->Fill(Q2(agtrack1,agtrack2));
   }
   //  c3) p+p- || p-p+ [4][9]:
   if((Proton(gtrack1,1,kTRUE) && Proton(gtrack2,-1,kTRUE)) || (Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,1,kTRUE)))
   {
    fCorrelationFunctions[4][9]->Fill(Q2(agtrack1,agtrack2));
   }

   // 2.) Mixed particle species:
   // a) pion-kaon
   //  a1) pi+K+ [2][3]:
   if(Pion(gtrack1,1,kTRUE) && Kaon(gtrack2,1,kTRUE))
   {
    fCorrelationFunctions[2][3]->Fill(Q2(agtrack1,agtrack2));
   }
   //  a2) pi+K- [2][8]
   if(Pion(gtrack1,1,kTRUE) && Kaon(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[2][8]->Fill(Q2(agtrack1,agtrack2));
   }
   //  a3) K+pi- [3][7]
   if(Kaon(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[3][7]->Fill(Q2(agtrack1,agtrack2));
   }
   //  a4) pi-K- [7][8]
   if(Pion(gtrack1,-1,kTRUE) && Kaon(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[7][8]->Fill(Q2(agtrack1,agtrack2));
   }
   // b) pion-proton
   //  b1) pi+p+ [2][4]:
   if(Pion(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE))
   {
    fCorrelationFunctions[2][4]->Fill(Q2(agtrack1,agtrack2));
   }
   //  b2) pi+p- [2][9]
   if(Pion(gtrack1,1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[2][9]->Fill(Q2(agtrack1,agtrack2));
   }
   //  b3) p+pi- [4][7]
   if(Proton(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[4][7]->Fill(Q2(agtrack1,agtrack2));
   }
   //  b4) pi-p- [7][9]
   if(Pion(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[7][9]->Fill(Q2(agtrack1,agtrack2));
   }
   // c) kaon-proton
   //  c1) K+p+ [3][4]:
   if(Kaon(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE))
   {
    fCorrelationFunctions[3][4]->Fill(Q2(agtrack1,agtrack2));
   }
   //  c2) K+p- [3][9]
   if(Kaon(gtrack1,1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[3][9]->Fill(Q2(agtrack1,agtrack2));
   }
   //  c3) p+K- [4][8]
   if(Proton(gtrack1,1,kTRUE) && Kaon(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[4][8]->Fill(Q2(agtrack1,agtrack2));
   }
   //  c4) K-p- [8][9]
   if(Kaon(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[8][9]->Fill(Q2(agtrack1,agtrack2));
   }
  } // for(Int_t iTrack2=0;iTrack2<nTracks;iTrack2++)
 } // for(Int_t iTrack1=0;iTrack1<nTracks;iTrack1++)

 // c) Three nested loops to calculate C(Q3), just an example; TBI
 if(fFill3pCorrelationFunctions) this->Calculate3pCorrelationFunctions(aAOD);

 // d) Four nested loops to calculate C(Q4), just an example. TBI
 if(fFill4pCorrelationFunctions) this->Calculate4pCorrelationFunctions(aAOD);

} // void AliAnalysisTaskMultiparticleFemtoscopy::CalculateCorrelationFunctions(AliAODEvent *aAOD)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::Calculate3pCorrelationFunctions(AliAODEvent *aAOD)
{
 // Calculate 3-particle correlation functions.

 // a) Insanity checks;
 // b) Three nested loops to calculate C(Q3), just an example.

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::Calculate3pCorrelationFunctions(AliAODEvent *aAOD)";
 if(!aAOD){Fatal(sMethodName.Data(),"!aAOD");}
 if(0 == fGlobalTracksAOD[0]->GetSize()){Fatal(sMethodName.Data(),"0 == fGlobalTracksAOD[0]->GetSize()");} // this case shall be already treated in UserExec

 // b) Three nested loops to calculate C(Q3), just an example:
 Int_t nTracks = aAOD->GetNumberOfTracks();
 for(Int_t iTrack1=0;iTrack1<nTracks;iTrack1++)
 {
  AliAODTrack *atrack1 = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack1));
  // TBI Temporary track insanity checks:
  if(!atrack1){Fatal(sMethodName.Data(),"!atrack1");} // TBI keep this for some time, eventually just continue
  if(atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(atrack1->TestFilterBit(128) && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->TestFiletrBit(128) && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(!PassesCommonTrackCuts(atrack1)){continue;} // TBI re-think
  // Corresponding AOD global track:
  Int_t id1 = atrack1->GetID();
  AliAODTrack *gtrack1 = dynamic_cast<AliAODTrack*>(id1>=0 ? aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(id1)) : aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(-(id1+1))));
  if(!gtrack1){Fatal(sMethodName.Data(),"!gtrack1");} // TBI keep this for some time, eventually just continue
  Int_t gid1 = (id1>=0 ? id1 : -(id1+1)); // ID of corresponding global track
  // Common track selection criteria for all "normal" global tracks:
  if(!PassesGlobalTrackCuts(gtrack1)){continue;}

  // Loop over the 2nd particle:
  for(Int_t iTrack2=0;iTrack2<nTracks;iTrack2++)
  {
   if(iTrack2<=iTrack1){continue;} // Eliminate self-evident self-correlations, and permutations as well
   AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack2));
   // TBI Temporary track insanity checks:
   if(!atrack2){Fatal(sMethodName.Data(),"!atrack2");} // TBI keep this for some time, eventually just continue
   if(atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(atrack2->TestFilterBit(128) && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->TestFiletrBit(128) && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(!PassesCommonTrackCuts(atrack2)){continue;} // TBI re-think
   // Corresponding AOD global track:
   Int_t id2 = atrack2->GetID();
   AliAODTrack *gtrack2 = dynamic_cast<AliAODTrack*>(id2>=0 ? aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(id2)) : aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(-(id2+1))));
   if(!gtrack2){Fatal(sMethodName.Data(),"!gtrack2");} // TBI keep this for some time, eventually just continue
   Int_t gid2 = (id2>=0 ? id2 : -(id2+1)); // ID of corresponding global track
   if(gid1==gid2){continue;} // Eliminate not-so-evident self-correlations:
   // Common track selection criteria for all "normal" global tracks:
   if(!PassesGlobalTrackCuts(gtrack2)){continue;}

   // Loop over the 3rd particle:
   for(Int_t iTrack3=0;iTrack3<nTracks;iTrack3++)
   {
    if(iTrack3<=iTrack2 || iTrack3<=iTrack1){continue;} // Eliminate self-evident self-correlations, and permutations as well
    AliAODTrack *atrack3 = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack3));
    // TBI Temporary track insanity checks:
    if(!atrack3){Fatal(sMethodName.Data(),"!atrack3");} // TBI keep this for some time, eventually just continue
    if(atrack3->GetID()>=0 && atrack3->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack3->GetID()>=0 && atrack3->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
    if(atrack3->TestFilterBit(128) && atrack3->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack3->TestFiletrBit(128) && atrack3->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
    if(!PassesCommonTrackCuts(atrack3)){continue;} // TBI re-think
    // Corresponding AOD global track:
    Int_t id3 = atrack3->GetID();
    AliAODTrack *gtrack3 = dynamic_cast<AliAODTrack*>(id3>=0 ? aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(id3)) : aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(-(id3+1))));
    if(!gtrack3){Fatal(sMethodName.Data(),"!gtrack3");} // TBI keep this for some time, eventually just continue
    Int_t gid3 = (id3>=0 ? id3 : -(id3+1)); // ID of corresponding global track
    if(gid3==gid2 || gid3==gid1){continue;} // Eliminate not-so-evident self-correlations
    // Common track selection criteria for all "normal" global tracks:
    if(!PassesGlobalTrackCuts(gtrack3)){continue;}

    // Okay, so we have three tracks, let's decide whether kinemetics from 'atrack' or 'gtrack' is taken, let's check PID, and fill the correlation functions:
    AliAODTrack *agtrack1 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
    AliAODTrack *agtrack2 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
    AliAODTrack *agtrack3 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
    if(fFillControlHistogramsWithGlobalTrackInfo)
    {
     agtrack1 = gtrack1;
     agtrack2 = gtrack2;
     agtrack3 = gtrack3;
    }
    else
    {
     agtrack1 = atrack1;
     agtrack2 = atrack2;
     agtrack3 = atrack3;
    } // if(fFillControlHistogramsWithGlobalTrackInfo)
    if(!agtrack1){Fatal(sMethodName.Data(),"!agtrack1");}
    if(!agtrack2){Fatal(sMethodName.Data(),"!agtrack2");}
    if(!agtrack3){Fatal(sMethodName.Data(),"!agtrack3");}

    // Cases of interest:
    // a) Same species and same charge;
    // b) Same species but different charge combinations (modulo permutations);
    // c) Two pions + something else (modulo permutations);
    // d) Two nucleons + something else (modulo permutations).

    // a) Same species:
    // pi+pi+pi+
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE))
    {
     f3pCorrelationFunctions[2][2][2]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi-pi-pi-
    if(Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE))
    {
     f3pCorrelationFunctions[7][7][7]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // K+K+K+
    if(Kaon(gtrack1,1,kTRUE) && Kaon(gtrack2,1,kTRUE) && Kaon(gtrack3,1,kTRUE))
    {
     f3pCorrelationFunctions[3][3][3]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // K-K-K-
    if(Kaon(gtrack1,-1,kTRUE) && Kaon(gtrack2,-1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pCorrelationFunctions[8][8][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p+p+p+
    if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE) && Proton(gtrack3,1,kTRUE))
    {
     f3pCorrelationFunctions[4][4][4]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p-p-p-
    if(Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE) && Proton(gtrack3,-1,kTRUE))
    {
     f3pCorrelationFunctions[9][9][9]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }

    // b) Same species but different charge combinations (modulo permutations):
    // pi+pi+pi-
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,-1,kTRUE))
    {
     f3pCorrelationFunctions[2][2][7]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi-pi-
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE))
    {
     f3pCorrelationFunctions[2][7][7]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // K+K+K-
    if(Kaon(gtrack1,1,kTRUE) && Kaon(gtrack2,1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pCorrelationFunctions[3][3][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // K+K-K-
    if(Kaon(gtrack1,1,kTRUE) && Kaon(gtrack2,-1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pCorrelationFunctions[3][8][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p+p+p-
    if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE) && Proton(gtrack3,-1,kTRUE))
    {
     f3pCorrelationFunctions[4][4][9]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p+p-p-
    if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,-1,kTRUE) && Proton(gtrack3,-1,kTRUE))
    {
     f3pCorrelationFunctions[4][9][9]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }

    // c) Two pions + something else (modulo permutations):
    // pi+pi+K+
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Kaon(gtrack3,1,kTRUE))
    {
     f3pCorrelationFunctions[2][2][3]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi+K-
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pCorrelationFunctions[2][2][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi-pi-K+
    if(Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Kaon(gtrack3,1,kTRUE))
    {
     f3pCorrelationFunctions[7][7][3]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi-pi-K-
    if(Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pCorrelationFunctions[7][7][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi-K+
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Kaon(gtrack3,1,kTRUE))
    {
     f3pCorrelationFunctions[2][7][3]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi-K-
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pCorrelationFunctions[2][7][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi+p+
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Proton(gtrack3,1,kTRUE))
    {
     f3pCorrelationFunctions[2][2][4]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi+p-
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Proton(gtrack3,-1,kTRUE))
    {
     f3pCorrelationFunctions[2][2][9]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi-pi-p+
    if(Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Proton(gtrack3,1,kTRUE))
    {
     f3pCorrelationFunctions[7][7][4]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi-pi-p-
    if(Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Proton(gtrack3,-1,kTRUE))
    {
     f3pCorrelationFunctions[7][7][9]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi-p+
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Proton(gtrack3,1,kTRUE))
    {
     f3pCorrelationFunctions[2][7][4]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi-p-
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Proton(gtrack3,-1,kTRUE))
    {
     f3pCorrelationFunctions[2][7][9]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }

    // d) Two nucleons + something else (modulo permutations):
    // p+p+pi+
    if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE))
    {
     f3pCorrelationFunctions[4][4][2]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p+p+pi-
    if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE) && Pion(gtrack3,-1,kTRUE))
    {
     f3pCorrelationFunctions[4][4][7]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p+p+K+
    if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE) && Kaon(gtrack3,1,kTRUE))
    {
     f3pCorrelationFunctions[4][4][3]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p+p+K-
    if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pCorrelationFunctions[4][4][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p-p-pi+
    if(Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE) && Pion(gtrack3,1,kTRUE))
    {
     f3pCorrelationFunctions[9][9][2]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p-p-pi-
    if(Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE))
    {
     f3pCorrelationFunctions[9][9][7]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p-p-K+
    if(Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE) && Kaon(gtrack3,1,kTRUE))
    {
     f3pCorrelationFunctions[9][9][3]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p-p-K-
    if(Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pCorrelationFunctions[9][9][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }

   } // for(Int_t iTrack3=0;iTrack3<nTracks;iTrack3++)
  } // for(Int_t iTrack3=0;iTrack3<nTracks;iTrack3++)
 } // for(Int_t iTrack3=0;iTrack3<nTracks;iTrack3++)

} // void AliAnalysisTaskMultiparticleFemtoscopy::Calculate3pCorrelationFunctions()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::Calculate4pCorrelationFunctions(AliAODEvent *aAOD)
{
 // Calculate 4-particle correlation functions.

 // a) Insanity checks;
 // b) Four nested loops to calculate C(Q4), just an example.

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::Calculate4pCorrelationFunctions(AliAODEvent *aAOD)";
 if(!aAOD){Fatal(sMethodName.Data(),"!aAOD");}
 if(0 == fGlobalTracksAOD[0]->GetSize()){Fatal(sMethodName.Data(),"0 == fGlobalTracksAOD[0]->GetSize()");} // this case shall be already treated in UserExec

 // b) Four nested loops to calculate C(Q4), just an example:
 Int_t nTracks = aAOD->GetNumberOfTracks();
 for(Int_t iTrack1=0;iTrack1<nTracks;iTrack1++)
 {
  AliAODTrack *atrack1 = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack1));
  // TBI Temporary track insanity checks:
  if(!atrack1){Fatal(sMethodName.Data(),"!atrack1");} // TBI keep this for some time, eventually just continue
  if(atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(atrack1->TestFilterBit(128) && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->TestFiletrBit(128) && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(!PassesCommonTrackCuts(atrack1)){continue;} // TBI re-think
  // Corresponding AOD global track:
  Int_t id1 = atrack1->GetID();
  AliAODTrack *gtrack1 = dynamic_cast<AliAODTrack*>(id1>=0 ? aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(id1)) : aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(-(id1+1))));
  if(!gtrack1){Fatal(sMethodName.Data(),"!gtrack1");} // TBI keep this for some time, eventually just continue
  Int_t gid1 = (id1>=0 ? id1 : -(id1+1)); // ID of corresponding global track
  // Common track selection criteria for all "normal" global tracks:
  if(!PassesGlobalTrackCuts(gtrack1)){continue;}

  // Loop over the 2nd particle:
  for(Int_t iTrack2=0;iTrack2<nTracks;iTrack2++)
  {
   if(iTrack1==iTrack2){continue;} // Eliminate self-evident self-correlations
   AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack2));
   // TBI Temporary track insanity checks:
   if(!atrack2){Fatal(sMethodName.Data(),"!atrack2");} // TBI keep this for some time, eventually just continue
   if(atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(atrack2->TestFilterBit(128) && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->TestFiletrBit(128) && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(!PassesCommonTrackCuts(atrack2)){continue;} // TBI re-think
   // Corresponding AOD global track:
   Int_t id2 = atrack2->GetID();
   AliAODTrack *gtrack2 = dynamic_cast<AliAODTrack*>(id2>=0 ? aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(id2)) : aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(-(id2+1))));
   if(!gtrack2){Fatal(sMethodName.Data(),"!gtrack2");} // TBI keep this for some time, eventually just continue
   Int_t gid2 = (id2>=0 ? id2 : -(id2+1)); // ID of corresponding global track
   if(gid1==gid2){continue;} // Eliminate not-so-evident self-correlations:
   // Common track selection criteria for all "normal" global tracks:
   if(!PassesGlobalTrackCuts(gtrack2)){continue;}

   // Loop over the 3rd particle:
   for(Int_t iTrack3=0;iTrack3<nTracks;iTrack3++)
   {
    if(iTrack3==iTrack2 || iTrack3==iTrack1){continue;} // Eliminate self-evident self-correlations
    AliAODTrack *atrack3 = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack3));
    // TBI Temporary track insanity checks:
    if(!atrack3){Fatal(sMethodName.Data(),"!atrack3");} // TBI keep this for some time, eventually just continue
    if(atrack3->GetID()>=0 && atrack3->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack3->GetID()>=0 && atrack3->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
    if(atrack3->TestFilterBit(128) && atrack3->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack3->TestFiletrBit(128) && atrack3->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
    if(!PassesCommonTrackCuts(atrack3)){continue;} // TBI re-think
    // Corresponding AOD global track:
    Int_t id3 = atrack3->GetID();
    AliAODTrack *gtrack3 = dynamic_cast<AliAODTrack*>(id3>=0 ? aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(id3)) : aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(-(id3+1))));
    if(!gtrack3){Fatal(sMethodName.Data(),"!gtrack3");} // TBI keep this for some time, eventually just continue
    Int_t gid3 = (id3>=0 ? id3 : -(id3+1)); // ID of corresponding global track
    if(gid3==gid2 || gid3==gid1){continue;} // Eliminate not-so-evident self-correlations
    // Common track selection criteria for all "normal" global tracks:
    if(!PassesGlobalTrackCuts(gtrack3)){continue;}

    // Loop over the 4th particle:
    for(Int_t iTrack4=0;iTrack4<nTracks;iTrack4++)
    {
     if(iTrack4==iTrack3 || iTrack4==iTrack2 || iTrack4==iTrack1){continue;} // Eliminate self-evident self-correlations
     AliAODTrack *atrack4 = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack4));
     // TBI Temporary track insanity checks:
     if(!atrack4){Fatal(sMethodName.Data(),"!atrack4");} // TBI keep this for some time, eventually just continue
     if(atrack4->GetID()>=0 && atrack4->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack4->GetID()>=0 && atrack4->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
     if(atrack4->TestFilterBit(128) && atrack4->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack4->TestFiletrBit(128) && atrack4->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
     if(!PassesCommonTrackCuts(atrack4)){continue;} // TBI re-think
     // Corresponding AOD global track:
     Int_t id4 = atrack4->GetID();
     AliAODTrack *gtrack4 = dynamic_cast<AliAODTrack*>(id4>=0 ? aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(id4)) : aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(-(id4+1))));
     if(!gtrack4){Fatal(sMethodName.Data(),"!gtrack4");} // TBI keep this for some time, eventually just continue
     Int_t gid4 = (id4>=0 ? id4 : -(id4+1)); // ID of corresponding global track
     if(gid4==gid3 || gid4==gid2 || gid4==gid1){continue;} // Eliminate not-so-evident self-correlations
     // Common track selection criteria for all "normal" global tracks:
     if(!PassesGlobalTrackCuts(gtrack4)){continue;}

     // Okay, so we have four tracks, let's decide whether kinemetics from 'atrack' or 'gtrack' is taken, let's check PID, and fill the correlation functions:
     AliAODTrack *agtrack1 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
     AliAODTrack *agtrack2 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
     AliAODTrack *agtrack3 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
     AliAODTrack *agtrack4 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
     if(fFillControlHistogramsWithGlobalTrackInfo)
     {
      agtrack1 = gtrack1;
      agtrack2 = gtrack2;
      agtrack3 = gtrack3;
      agtrack4 = gtrack4;
     }
     else
     {
      agtrack1 = atrack1;
      agtrack2 = atrack2;
      agtrack3 = atrack3;
      agtrack4 = atrack4;
     } // if(fFillControlHistogramsWithGlobalTrackInfo)
     if(!agtrack1){Fatal(sMethodName.Data(),"!agtrack1");}
     if(!agtrack2){Fatal(sMethodName.Data(),"!agtrack2");}
     if(!agtrack3){Fatal(sMethodName.Data(),"!agtrack3");}
     if(!agtrack4){Fatal(sMethodName.Data(),"!agtrack4");}


     // TBI
     // First test example: pi+pi+pi+pi+
     if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE) && Pion(gtrack4,1,kTRUE))
     {
      f4pCorrelationFunctions[2][2][2][2]->Fill(Q4(agtrack1,agtrack2,agtrack3,agtrack4)); // Lorentz invariant Q4
     }

     // Second test example: pi-pi-pi-pi-
     if(Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE) && Pion(gtrack4,-1,kTRUE))
     {
      f4pCorrelationFunctions[7][7][7][7]->Fill(Q4(agtrack1,agtrack2,agtrack3,agtrack4)); // Lorentz invariant Q4
     }


    } // for(Int_t iTrack4=0;iTrack4<nTracks;iTrack4++)
   } // for(Int_t iTrack3=0;iTrack3<nTracks;iTrack3++)
  } // for(Int_t iTrack2=0;iTrack2<nTracks;iTrack2++)
 } // for(Int_t iTrack1=0;iTrack1<nTracks;iTrack1++)

} // void AliAnalysisTaskMultiparticleFemtoscopy::Calculate4pCorrelationFunctions()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::CalculateCorrelationFunctions(AliMCEvent *aMC)
{
 // Calculate correlation functions for Monte Carlo.

 // a) Insanity checks;
 // b) Pattern from supported PDG codes;
 // c) Two nested loops to calculate C(k), just an example.

 // b) Pattern from supported PDG codes:
 TString pattern = ".11.13.211.321.2212."; // TBI this can be done programatically

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::CalculateCorrelationFunctions(AliMCEvent *aMC)";
 if(0 == aMC->GetNumberOfTracks()){return;} // TBI re-think

 // c) Two nested loops to calculate C(k), just an example:
 Int_t nTracks = aMC->GetNumberOfTracks();
 for(Int_t iTrack1=0;iTrack1<nTracks;iTrack1++)
 {
  AliAODMCParticle *amcparticle1 = dynamic_cast<AliAODMCParticle*>(aMC->GetTrack(iTrack1));
  if(!amcparticle1){Fatal(sMethodName.Data(),"!amcparticle1");}

  if(!PassesCommonTrackCuts(amcparticle1)){continue;} // TBI re-think, see implemntation of this method

  // Loop over the 2nd particle:
  for(Int_t iTrack2=iTrack1+1;iTrack2<nTracks;iTrack2++)
  {
   AliAODMCParticle *amcparticle2 = dynamic_cast<AliAODMCParticle*>(aMC->GetTrack(iTrack2));
   if(!amcparticle2){Fatal(sMethodName.Data(),"!amcparticle2");}

   if(!PassesCommonTrackCuts(amcparticle2)){continue;} // TBI re-think, see implemntation of this method

   // Okay, so we have two tracks, let's check PID, and fill the correlation functions:
   // Check if this PID is supported in current implementation:
   if(!(pattern.Contains(Form(".%d.",TMath::Abs(amcparticle1->GetPdgCode()))) && pattern.Contains(Form(".%d.",TMath::Abs(amcparticle2->GetPdgCode()))))){continue;}

   // Determine the indices of correlation function to be filled:
   Int_t index1 = -44;
   Int_t index2 = -44;
   if(fCorrelationFunctionsIndices->GetValue(amcparticle1->GetPdgCode())<=fCorrelationFunctionsIndices->GetValue(amcparticle2->GetPdgCode()))
   {
    index1 = fCorrelationFunctionsIndices->GetValue(amcparticle1->GetPdgCode());
    index2 = fCorrelationFunctionsIndices->GetValue(amcparticle2->GetPdgCode());
   }
   else
   {
    index1 = fCorrelationFunctionsIndices->GetValue(amcparticle2->GetPdgCode());
    index2 = fCorrelationFunctionsIndices->GetValue(amcparticle1->GetPdgCode());
   }
   if(-44 == index1){Fatal(sMethodName.Data(),"-44 == index1");}
   if(-44 == index2){Fatal(sMethodName.Data(),"-44 == index2");}

   // Fill or die:
   fCorrelationFunctions[index1][index2]->Fill(RelativeMomenta(amcparticle1,amcparticle2)); // for the relative momenta ordering doesn't matter

  } // for(Int_t iTrack2=iTrack1+1;iTrack2<nTracks;iTrack2++)
 } // for(Int_t iTrack1=0;iTrack1<nTracks;iTrack1++)

} // void AliAnalysisTaskMultiparticleFemtoscopy::CalculateCorrelationFunctions(AliMCEvent *aMC)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::CalculateCorrelationFunctionsTEST(AliAODEvent *aAOD)
{
 // Calculate correlation functions.

 // a) Insanity checks;
 // b) Book local TProfile's to hold single-event averages vs. Q2;
 // c) Book local TProfile's to hold single-event averages vs. Q3;
 // d) Nested loops to calculate single-event averages; // Remark: To evaluate 3p loops, enable the flag fFill3pCorrelationFunctions;
 // e) Build all-event correlations and cumulants from single-event averages;
 // f) Delete local TProfile's to hold single-event averages.

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::CalculateCorrelationFunctionsTEST(AliAODEvent *aAOD)";
 if(!aAOD){Fatal(sMethodName.Data(),"!aAOD");}
 if(0 == fGlobalTracksAOD[0]->GetSize()){Fatal(sMethodName.Data(),"0 == fGlobalTracksAOD[0]->GetSize()");} // this case shall be already treated in UserExec

 // b) Book local TProfile's to hold single-event averages vs. k:
 const Int_t nTestsMax = 10; // see what is hardwired in .h file
 TString sXYZ[3] = {"x","y","z"};
 const Int_t n2pCumulantTerms = 3;
 TString s2pCumulantTerms[n2pCumulantTerms] = {"#LTX_{1}#GT","#LTX_{2}#GT","#LTX_{1}X_{2}#GT"};
 TProfile *singleEventAverageCorrelationsVsQ2[nTestsMax][n2pCumulantTerms][3] = {{{NULL}}}; // [3 = x,y,z components]
 for(Int_t t=0;t<nTestsMax;t++) // test No
 {
  if(!fFillCorrelationFunctionsTEST[t]){continue;}
  for(Int_t ct=0;ct<n2pCumulantTerms;ct++)
  {
   for(Int_t xyz=0;xyz<3;xyz++)
   {
    singleEventAverageCorrelationsVsQ2[t][ct][xyz] = new TProfile(Form("singleEventAverageCorrelationsVsQ2[%d][%d][%d]",t,ct,xyz),"single-event averages",fnQ2bins,fnQ2min,fnQ2max);
    singleEventAverageCorrelationsVsQ2[t][ct][xyz]->GetXaxis()->SetTitle("Q_{2}");
    singleEventAverageCorrelationsVsQ2[t][ct][xyz]->GetYaxis()->SetTitle(Form("%s_{%s}",s2pCumulantTerms[ct].Data(),sXYZ[xyz].Data()));
   } // for(Int_t xyz=0;xyz<3;xyz++)
  } // for(Int_t ct=0;ct<n2pCumulantTerms;ct++)
 } // for(Int_t t=0;t<nTestsMax;t++) // test No

 // c) Book local TProfile's to hold single-event averages vs. Q3:
 const Int_t n3pCumulantTerms = 7;
 TString s3pCumulantTerms[n3pCumulantTerms] = {"#LTX_{1}#GT","#LTX_{2}#GT","#LTX_{3}#GT","#LTX_{1}X_{2}#GT","#LTX_{1}X_{3}#GT","#LTX_{2}X_{3}#GT","#LTX_{1}X_{2}X_{3}#GT"};
 TProfile *singleEventAverageCorrelationsVsQ3[nTestsMax][n3pCumulantTerms][3] = {{{NULL}}}; // [3 = x,y,z components]
 for(Int_t t=0;t<nTestsMax;t++) // test No
 {
  if(!fFillCorrelationFunctionsTEST[t]){continue;}
  for(Int_t ct=0;ct<n3pCumulantTerms;ct++)
  {
   if(!fFill3pCorrelationFunctions){break;}
   for(Int_t xyz=0;xyz<3;xyz++)
   {
    singleEventAverageCorrelationsVsQ3[t][ct][xyz] = new TProfile(Form("singleEventAverageCorrelationsVsQ3[%d][%d][%d]",t,ct,xyz),"single-event averages",fnQ3bins,fnQ3min,fnQ3max);
    singleEventAverageCorrelationsVsQ3[t][ct][xyz]->GetXaxis()->SetTitle("Q_{3}");
    singleEventAverageCorrelationsVsQ3[t][ct][xyz]->GetYaxis()->SetTitle(Form("%s_{%s}",s3pCumulantTerms[ct].Data(),sXYZ[xyz].Data()));
   } // for(Int_t xyz=0;xyz<3;xyz++)
  } // for(Int_t ct=0;ct<n3pCumulantTerms;ct++)
 } // for(Int_t t=0;t<nTestsMax;t++) // test No

 // ) Differential yield counters (needed only for "Test 2" at the moment):
 const Int_t nLoopsVsQ2 = 2;
 const Int_t nLoopsVsQ3 = 3;
 const Int_t nMaxEntries = 10000;
 TString *dycVsQ2[nLoopsVsQ2][nMaxEntries] = {{NULL}};
 TString *dycVsQ3[nLoopsVsQ3][nMaxEntries] = {{NULL}};
 for(Int_t nl=0;nl<nLoopsVsQ2;nl++)
 {
  for(Int_t nme=1;nme<=nMaxEntries;nme++) // yes, from 1, since I do not have later bin 0
  {
   if(nme > fnQ2bins){break;}
   dycVsQ2[nl][nme] = new TString();
  }
 }
 for(Int_t nl=0;nl<nLoopsVsQ3;nl++)
 {
  for(Int_t nme=1;nme<=nMaxEntries;nme++) // yes, from 1, since I do not have later bin 0
  {
   if(nme > fnQ3bins){break;}
   dycVsQ3[nl][nme] = new TString();
  }
 }
 
 // Counters for N_\pi+ and N_\pi^-, for each Q2 bin (needed only for "Test 7" at the moment):
 TString *test7Counters[10000][2] = {{NULL}}; // N[44][0] is then N_\pi+ in the 44th Q2 bin in this event, while N[44][1] is then N_\pi- in the 44th Q2 bin in this event TBI hardwired 10000
 for(Int_t q2=0;q2<fnQ2bins;q2++)
 {
  test7Counters[q2][0] = new TString();
  *test7Counters[q2][0] = "";
  
  test7Counters[q2][1] = new TString();
  *test7Counters[q2][1] = "";
 }

 // d) Nested loops to calculate single-event averages:
 Int_t nTracks = aAOD->GetNumberOfTracks();
 for(Int_t iTrack1=0;iTrack1<nTracks;iTrack1++)
 {
  AliAODTrack *atrack1 = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack1));
  // TBI Temporary track insanity checks:
  if(!atrack1){Fatal(sMethodName.Data(),"!atrack1");} // TBI keep this for some time, eventually just continue
  if(atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(atrack1->TestFilterBit(128) && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->TestFiletrBit(128) && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(!PassesCommonTrackCuts(atrack1)){continue;} // TBI re-think
  // Corresponding AOD global track:
  Int_t id1 = atrack1->GetID();
  AliAODTrack *gtrack1 = dynamic_cast<AliAODTrack*>(id1>=0 ? aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(id1)) : aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(-(id1+1))));
  if(!gtrack1){Fatal(sMethodName.Data(),"!gtrack1");} // TBI keep this for some time, eventually just continue
  Int_t gid1 = (id1>=0 ? id1 : -(id1+1)); // ID of corresponding global track
  // Common track selection criteria for all "normal" global tracks:
  if(!PassesGlobalTrackCuts(gtrack1)){continue;}

  // Loop over the 2nd particle:
  for(Int_t iTrack2=0;iTrack2<nTracks;iTrack2++) // TBI: I can gain in performance if I do not start from 0
  {
   //if(iTrack2<=iTrack1){continue;} // Eliminate self-evident self-correlations, and permutations as well
   if(iTrack2==iTrack1){continue;} // Eliminate self-evident self-correlations
   AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack2));
   // TBI Temporary track insanity checks:
   if(!atrack2){Fatal(sMethodName.Data(),"!atrack2");} // TBI keep this for some time, eventually just continue
   if(atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(atrack2->TestFilterBit(128) && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->TestFiletrBit(128) && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(!PassesCommonTrackCuts(atrack2)){continue;} // TBI re-think
   // Corresponding AOD global track:
   Int_t id2 = atrack2->GetID();
   AliAODTrack *gtrack2 = dynamic_cast<AliAODTrack*>(id2>=0 ? aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(id2)) : aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(-(id2+1))));
   if(!gtrack2){Fatal(sMethodName.Data(),"!gtrack2");} // TBI keep this for some time, eventually just continue
   Int_t gid2 = (id2>=0 ? id2 : -(id2+1)); // ID of corresponding global track
   if(gid1==gid2){continue;} // Eliminate self-correlations:
   // Common track selection criteria for all "normal" global tracks:
   if(!PassesGlobalTrackCuts(gtrack2)){continue;}

   // Okay, so we have two tracks, let's decide whether kinemetics from 'atrack' or 'gtrack' is taken, let's check PID, and fill the single event correlation functions projected onto Q2:
   AliAODTrack *agtrack1 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
   AliAODTrack *agtrack2 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
   if(fFillControlHistogramsWithGlobalTrackInfo)
   {
    agtrack1 = gtrack1;
    agtrack2 = gtrack2;
   }
   else
   {
    agtrack1 = atrack1;
    agtrack2 = atrack2;
   } // if(fFillControlHistogramsWithGlobalTrackInfo)
   if(!agtrack1){Fatal(sMethodName.Data(),"!agtrack1");}
   if(!agtrack2){Fatal(sMethodName.Data(),"!agtrack2");}

   // Test 0: "Same charge pions, 2p correlations and cumulants projected onto k, for x, y and z components separately", not Lorentz invariant
   if(fFillCorrelationFunctionsTEST[0])
   {
    if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE)))
    {
     Double_t dk = RelativeMomenta(agtrack1,agtrack2); // relative momenta k = \frac{1}{2}|\vec{p_1}-\vec{p_2}|

     fSignalYieldTEST[0]->Fill(dk); // TBI temporarily hardwired here

     // x-components:
     singleEventAverageCorrelationsVsQ2[0][0][0]->Fill(dk,agtrack1->Px()); // <X1>_x
     singleEventAverageCorrelationsVsQ2[0][1][0]->Fill(dk,agtrack2->Px()); // <X2>_x
     singleEventAverageCorrelationsVsQ2[0][2][0]->Fill(dk,agtrack1->Px()*agtrack2->Px()); // <X1X2>_x

     // y-components:
     singleEventAverageCorrelationsVsQ2[0][0][1]->Fill(dk,agtrack1->Py()); // <X1>_y
     singleEventAverageCorrelationsVsQ2[0][1][1]->Fill(dk,agtrack2->Py()); // <X2>_y
     singleEventAverageCorrelationsVsQ2[0][2][1]->Fill(dk,agtrack1->Py()*agtrack2->Py()); // <X1X2>_y

     // z-components:
     singleEventAverageCorrelationsVsQ2[0][0][2]->Fill(dk,agtrack1->Pz()); // <X1>_z
     singleEventAverageCorrelationsVsQ2[0][1][2]->Fill(dk,agtrack2->Pz()); // <X2>_z
     singleEventAverageCorrelationsVsQ2[0][2][2]->Fill(dk,agtrack1->Pz()*agtrack2->Pz()); // <X1X2>_z

    } // if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE)))
   } // if(fFillCorrelationFunctionsTEST[0])

   // Test 1: "Same charge pions, 2p correlations and cumulants projected onto Lorentz invariant Q2"
   if(fFillCorrelationFunctionsTEST[1])
   {
    if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE)))
    {

     Double_t dQ2 = Q2(agtrack1,agtrack2); // Lorentz invariant Q2

     singleEventAverageCorrelationsVsQ2[1][0][0]->Fill(dQ2,agtrack1->Px()); // <X1>_x
     singleEventAverageCorrelationsVsQ2[1][1][0]->Fill(dQ2,agtrack2->Px()); // <X2>_x
     singleEventAverageCorrelationsVsQ2[1][2][0]->Fill(dQ2,agtrack1->Px()*agtrack2->Px()); // <X1X2>_x

     // y-components:
     singleEventAverageCorrelationsVsQ2[1][0][1]->Fill(dQ2,agtrack1->Py()); // <X1>_y
     singleEventAverageCorrelationsVsQ2[1][1][1]->Fill(dQ2,agtrack2->Py()); // <X2>_y
     singleEventAverageCorrelationsVsQ2[1][2][1]->Fill(dQ2,agtrack1->Py()*agtrack2->Py()); // <X1X2>_y

     // z-components:
     singleEventAverageCorrelationsVsQ2[1][0][2]->Fill(dQ2,agtrack1->Pz()); // <X1>_z
     singleEventAverageCorrelationsVsQ2[1][1][2]->Fill(dQ2,agtrack2->Pz()); // <X2>_z
     singleEventAverageCorrelationsVsQ2[1][2][2]->Fill(dQ2,agtrack1->Pz()*agtrack2->Pz()); // <X1X2>_z

    } // if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE)))
   } // if(fFillCorrelationFunctionsTEST[1])

   // Test 4: "Same charge pions, 2p correlations and cumulants projected onto Lorentz invariant Q2, where energy prefactors have been taken into account, to make correlations and cumulants manifestly Lorenz invariant, even before projection"
   if(fFillCorrelationFunctionsTEST[4])
   {
    if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE)))
    {

     Double_t E1 = agtrack1->E(); // energy of 1st track TBI check the mass hypothesis
     Double_t E2 = agtrack2->E(); // energy of 2nd track TBI check the mass hypothesis
     Double_t dQ2 = Q2(agtrack1,agtrack2); // Lorentz invariant Q2

     singleEventAverageCorrelationsVsQ2[4][0][0]->Fill(dQ2,E1); // <X1>
     singleEventAverageCorrelationsVsQ2[4][1][0]->Fill(dQ2,E2); // <X2>
     singleEventAverageCorrelationsVsQ2[4][2][0]->Fill(dQ2,E1*E2); // <X1X2>

    } // if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE)))
   } // if(fFillCorrelationFunctionsTEST[4])

   // Test 5: "Same charge pions, 2p correlations and cumulants projected onto Lorentz invariant Q2, where energy prefactors have been taken into account, to make correlations and cumulants manifestly Lorenz invariant, even before projection. Here I am boosting all particles to a new frame, to check whether the final results are Lorentz invariant:"
   if(fFillCorrelationFunctionsTEST[5])
   {
    if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE)))
    {
     // p_1:
     Double_t p1x = agtrack1->Px();
     Double_t p1y = agtrack1->Py();
     Double_t p1z = agtrack1->Pz();
     Double_t e1  = agtrack1->E();
     // p_2:
     Double_t p2x = agtrack2->Px();
     Double_t p2y = agtrack2->Py();
     Double_t p2z = agtrack2->Pz();
     Double_t e2  = agtrack2->E();  

     // Corresponding energy-momentum four-vectors:
     TLorentzVector track1(p1x,p1y,p1z,e1);
     TLorentzVector track2(p2x,p2y,p2z,e2);

     // Boosting both energy-momentum four-vectors:
     // N.B. I have a flag fBoost, for future cases...
     track1.Boost(fBoostVelocity);
     track2.Boost(fBoostVelocity);

     Double_t E1 = track1.E(); // energy of 1st track TBI check the mass hypothesis
     Double_t E2 = track2.E(); // energy of 2nd track TBI check the mass hypothesis
     Double_t dQ2 = Q2(agtrack1,agtrack2); // Lorentz invariant Q2

     singleEventAverageCorrelationsVsQ2[5][0][0]->Fill(dQ2,E1); // <X1>
     singleEventAverageCorrelationsVsQ2[5][1][0]->Fill(dQ2,E2); // <X2>
     singleEventAverageCorrelationsVsQ2[5][2][0]->Fill(dQ2,E1*E2); // <X1X2>

    } // if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE)))
   } // if(fFillCorrelationFunctionsTEST[5])

   // Test 6: "Same charge pions, E_{a}^{(b)} projected onto Lorentz invariant Q2 (see the write-up). This distribution shall be manifestly Lorentz invariant!"
   if(fFillCorrelationFunctionsTEST[6])
   {
    if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE)))
    {
     // p_1:
     Double_t p1x = agtrack1->Px();
     Double_t p1y = agtrack1->Py();
     Double_t p1z = agtrack1->Pz();
     Double_t e1  = agtrack1->E();
     // p_2:
     Double_t p2x = agtrack2->Px();
     Double_t p2y = agtrack2->Py();
     Double_t p2z = agtrack2->Pz();
     Double_t e2  = agtrack2->E();  

     // Corresponding energy-momentum four-vectors:
     TLorentzVector track1(p1x,p1y,p1z,e1);
     TLorentzVector track2(p2x,p2y,p2z,e2);

     // Calculating L.I. E_{a}^{(b)} (i.e. energy of the 1st particle in the (rest) frame of the 2nd):
     TVector3 vb = TVector3(p2x/e2,p2y/e2,p2z/e2); // velocity of 2nd particle
     track1.Boost(-vb); // boosting the 1st track in the rest frame of 2nd
     //track2.Boost(-vb); // boosting the 2nd track in its rest frame 
     Double_t Eab = track1.E(); // energy of 1st track in the rest frame of 2nd. Yes, this is L.I. quantity

     // Fill Eab vs. dQ2 distribution, which is then manifestly L.I.
     Double_t dQ2 = Q2(agtrack1,agtrack2); // Lorentz invariant Q2
     singleEventAverageCorrelationsVsQ2[6][0][0]->Fill(dQ2,Eab); // <E_{a}^{(b)}> vs. Q2

     fEab_TEST6[0]->Fill(dQ2,Eab);

    } // if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE)))
   } // if(fFillCorrelationFunctionsTEST[6])

   // Test 7: "Differential yield, see the write-up. Instead of counting # of pairs in Q2, I am counting # of particle "a" and # of particle "b" separately. 
   //          Then, the two independent observables are N_a and N_b, and they fluctuate e-b-e. That's how I get the notion of cumulants, and the projection 
   //          onto Q2 shall render it manifestly Lorentz invariant!"
   // To do:
   // a) Condidering in this test \pi+ \pi- only, later this will be generalized
   // b) Check explicitly the Lorentz invariance 
   // c) This clearly will work only for non-identical particles :'( Which is perhaps not than bad, given the fact that I can start to do my thing already with Lp+p-      
   // d) I cannot do anything for a given event, as it seems... :'( 
   if(fFillCorrelationFunctionsTEST[7])
   {
    //if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE)) // TBI 20170724
    if( (Pion(gtrack1,1,kTRUE) || Pion(gtrack1,-1,kTRUE)) && (Proton(gtrack2,1,kTRUE) || Proton(gtrack2,-1,kTRUE)) )
    {
     // p_1:
     Double_t p1x = agtrack1->Px();
     Double_t p1y = agtrack1->Py();
     Double_t p1z = agtrack1->Pz();
     Double_t e1  = agtrack1->E();
     // p_2:
     Double_t p2x = agtrack2->Px();
     Double_t p2y = agtrack2->Py();
     Double_t p2z = agtrack2->Pz();
     Double_t e2  = agtrack2->E();  

     // Corresponding energy-momentum four-vectors:
     TLorentzVector track1(p1x,p1y,p1z,e1);
     TLorentzVector track2(p2x,p2y,p2z,e2);

     // Pseudo-code: 
     // 0) for each Q2 bin, formed from p_1 and p_2, count number of unique \pi+ and \pi-
     // 1) Then, treat N_\pi+ and N_\pi- and two independent observables, and after I run over all pairs, fill TProfile to get <N_\pi+ N_\pi-> - <N_\pi+> <N_\pi-> for each Q2 bin
     // 2) See if there is something non-trivial happening only at low Q2
     // 3) Do the same thing for background
     // 4) Generalize for 3-particles

     Double_t dQ2 = Q2(agtrack1,agtrack2); // Lorentz invariant Q2

     Int_t nBinQ2Number = this->BinNoForSpecifiedValue(fCorrelationFunctionsTEST[7][0][0][0],dQ2);
     //cout<<"BIN_NO: "<<nBinQ2Number<<endl;
     //cout<<"ENERGIES: "<<agtrack1->E()<<" "<<agtrack2->E()<<endl;
     if(!test7Counters[nBinQ2Number][0]->Contains(Form("_%.6f_",agtrack1->E())))
     {
      *test7Counters[nBinQ2Number][0]+=Form("_%.6f_",agtrack1->E());
     }
     //cout<<"STRING1: "<<test7Counters[nBinQ2Number][0]->Data()<<endl;

     if(!test7Counters[nBinQ2Number][1]->Contains(Form("_%.6f_",agtrack2->E())))
     {
      *test7Counters[nBinQ2Number][1]+=Form("_%.6f_",agtrack2->E());
     }
     //cout<<"STRING2: "<<test7Counters[nBinQ2Number][1]->Data()<<endl;


//  TProfile *fCorrelationFunctionsTEST[10][2][7][10]; //! [testNo][0=vs Q2, 1=vs Q3][example [0=<x1>][1=<x2>], ...,[6=<x1x2x3>]][differential index, e.g. for test 0 [0=Cx][1=Cy][2=Cz]]


/*

 TString *test7Counters[10000][2] = {{NULL}}; // N[44][0] is then N_\pi+ in the 44th Q2 bin in this event, while N[44][1] is then N_\pi- in the 44th Q2 bin in this event TBI hardwired 10000
 for(Int_t q2=0;q2<fnQ2bins;q2++)
 {
  test7Counters[q2][0] = new TString();
  *test7Counters[q2][0] = "";
  
  test7Counters[q2][1] = new TString();
  *test7Counters[q2][1] = "";
 }

*/

/*
     // Calculating L.I. E_{a}^{(b)} (i.e. energy of the 1st particle in the (rest) frame of the 2nd):
     TVector3 vb = TVector3(p2x/e2,p2y/e2,p2z/e2); // velocity of 2nd particle
     track1.Boost(-vb); // boosting the 1st track in the rest frame of 2nd
     //track2.Boost(-vb); // boosting the 2nd track in its rest frame 
     Double_t Eab = track1.E(); // energy of 1st track in the rest frame of 2nd. Yes, this is L.I. quantity

     // Fill Eab vs. dQ2 distribution, which is then manifestly L.I.
     singleEventAverageCorrelationsVsQ2[6][0][0]->Fill(dQ2,Eab); // <E_{a}^{(b)}> vs. Q2
*/
 

    } // if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE))
   } // if(fFillCorrelationFunctionsTEST[7])


   // Loop over the 3rd particle:
   for(Int_t iTrack3=0;iTrack3<nTracks;iTrack3++)
   {
    if(!fFill3pCorrelationFunctions){break;}
    //if(iTrack3<=iTrack2 || iTrack3<=iTrack1){continue;} // Eliminate self-evident self-correlations, and permutations as well
    if(iTrack3==iTrack2 || iTrack3==iTrack1){continue;} // Eliminate self-evident self-correlations
    AliAODTrack *atrack3 = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack3));
    if(!atrack3){Fatal(sMethodName.Data(),"!atrack3");} // TBI keep this for some time, eventually just continue
    if(atrack3->GetID()>=0 && atrack3->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack3->GetID()>=0 && atrack3->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
    if(atrack3->TestFilterBit(128) && atrack3->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack3->TestFiletrBit(128) && atrack3->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
    if(!PassesCommonTrackCuts(atrack3)){continue;} // TBI re-think
    // Corresponding AOD global track:
    Int_t id3 = atrack3->GetID();
    AliAODTrack *gtrack3 = dynamic_cast<AliAODTrack*>(id3>=0 ? aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(id3)) : aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(-(id3+1))));
    if(!gtrack3){Fatal(sMethodName.Data(),"!gtrack3");} // TBI keep this for some time, eventually just continue
    Int_t gid3 = (id3>=0 ? id3 : -(id3+1)); // ID of corresponding global track
    if(gid3==gid2 || gid3==gid1){continue;} // Eliminate not-so-evident self-correlations
    // Common track selection criteria for all "normal" global tracks:
    if(!PassesGlobalTrackCuts(gtrack3)){continue;}

    // Okay, so we have the third track, let's decide whether kinemetics from 'atrack' or 'gtrack' is taken, let's check PID, and fill the correlation functions:
    AliAODTrack *agtrack3 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
    if(fFillControlHistogramsWithGlobalTrackInfo)
    {
     agtrack3 = gtrack3;
    }
    else
    {
     agtrack3 = atrack3;
    } // if(fFillControlHistogramsWithGlobalTrackInfo)
    if(!agtrack3){Fatal(sMethodName.Data(),"!agtrack3");}

    // Test 0: "Same charge pions, 3p correlations and cumulants projected onto Q3, for x, y and z components separately", Lorentz invariant, as of 20170319 TBI
    if(fFillCorrelationFunctionsTEST[0])
    {
     if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE)))
     {
      // Q3:
      Double_t dQ3 = Q3(agtrack1,agtrack2,agtrack3); // TBI this is NOT Lorentz invariant, use Q3_NEW instead

      fSignalYieldTEST[1]->Fill(dQ3); // TBI temporarily hardwired here

      // x-components:
      singleEventAverageCorrelationsVsQ3[0][0][0]->Fill(dQ3,agtrack1->Px()); // <X1>_x
      singleEventAverageCorrelationsVsQ3[0][1][0]->Fill(dQ3,agtrack2->Px()); // <X2>_x
      singleEventAverageCorrelationsVsQ3[0][2][0]->Fill(dQ3,agtrack3->Px()); // <X3>_x
      singleEventAverageCorrelationsVsQ3[0][3][0]->Fill(dQ3,agtrack1->Px()*agtrack2->Px()); // <X1X2>_x
      singleEventAverageCorrelationsVsQ3[0][4][0]->Fill(dQ3,agtrack1->Px()*agtrack3->Px()); // <X1X3>_x
      singleEventAverageCorrelationsVsQ3[0][5][0]->Fill(dQ3,agtrack2->Px()*agtrack3->Px()); // <X2X3>_x
      singleEventAverageCorrelationsVsQ3[0][6][0]->Fill(dQ3,agtrack1->Px()*agtrack2->Px()*agtrack3->Px()); // <X1X2X3>_x

      // y-components:
      singleEventAverageCorrelationsVsQ3[0][0][1]->Fill(dQ3,agtrack1->Py()); // <X1>_y
      singleEventAverageCorrelationsVsQ3[0][1][1]->Fill(dQ3,agtrack2->Py()); // <X2>_y
      singleEventAverageCorrelationsVsQ3[0][2][1]->Fill(dQ3,agtrack3->Py()); // <X3>_y
      singleEventAverageCorrelationsVsQ3[0][3][1]->Fill(dQ3,agtrack1->Py()*agtrack2->Py()); // <X1X2>_y
      singleEventAverageCorrelationsVsQ3[0][4][1]->Fill(dQ3,agtrack1->Py()*agtrack3->Py()); // <X1X3>_y
      singleEventAverageCorrelationsVsQ3[0][5][1]->Fill(dQ3,agtrack2->Py()*agtrack3->Py()); // <X2X3>_y
      singleEventAverageCorrelationsVsQ3[0][6][1]->Fill(dQ3,agtrack1->Py()*agtrack2->Py()*agtrack3->Py()); // <X1X2X3>_y

      // z-components:
      singleEventAverageCorrelationsVsQ3[0][0][2]->Fill(dQ3,agtrack1->Pz()); // <X1>_z
      singleEventAverageCorrelationsVsQ3[0][1][2]->Fill(dQ3,agtrack2->Pz()); // <X2>_z
      singleEventAverageCorrelationsVsQ3[0][2][2]->Fill(dQ3,agtrack3->Pz()); // <X3>_z
      singleEventAverageCorrelationsVsQ3[0][3][2]->Fill(dQ3,agtrack1->Pz()*agtrack2->Pz()); // <X1X2>_z
      singleEventAverageCorrelationsVsQ3[0][4][2]->Fill(dQ3,agtrack1->Pz()*agtrack3->Pz()); // <X1X3>_z
      singleEventAverageCorrelationsVsQ3[0][5][2]->Fill(dQ3,agtrack2->Pz()*agtrack3->Pz()); // <X2X3>_z
      singleEventAverageCorrelationsVsQ3[0][6][2]->Fill(dQ3,agtrack1->Pz()*agtrack2->Pz()*agtrack3->Pz()); // <X1X2X3>_z

     } // if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE)))
    } // if(fFillCorrelationFunctionsTEST[0])

    // Test 1: "Same charge pions, 2p correlations and cumulants projected onto Lorentz invariant Q3""
    if(fFillCorrelationFunctionsTEST[1])
    {
     if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE)))
     {

      // Q3:
      Double_t dQ3 = Q3(agtrack1,agtrack2,agtrack3); // Lorentz invariant Q3

      // x-components:
      singleEventAverageCorrelationsVsQ3[1][0][0]->Fill(dQ3,agtrack1->Px()); // <X1>_x
      singleEventAverageCorrelationsVsQ3[1][1][0]->Fill(dQ3,agtrack2->Px()); // <X2>_x
      singleEventAverageCorrelationsVsQ3[1][2][0]->Fill(dQ3,agtrack3->Px()); // <X3>_x
      singleEventAverageCorrelationsVsQ3[1][3][0]->Fill(dQ3,agtrack1->Px()*agtrack2->Px()); // <X1X2>_x
      singleEventAverageCorrelationsVsQ3[1][4][0]->Fill(dQ3,agtrack1->Px()*agtrack3->Px()); // <X1X3>_x
      singleEventAverageCorrelationsVsQ3[1][5][0]->Fill(dQ3,agtrack2->Px()*agtrack3->Px()); // <X2X3>_x
      singleEventAverageCorrelationsVsQ3[1][6][0]->Fill(dQ3,agtrack1->Px()*agtrack2->Px()*agtrack3->Px()); // <X1X2X3>_x

      // y-components:
      singleEventAverageCorrelationsVsQ3[1][0][1]->Fill(dQ3,agtrack1->Py()); // <X1>_y
      singleEventAverageCorrelationsVsQ3[1][1][1]->Fill(dQ3,agtrack2->Py()); // <X2>_y
      singleEventAverageCorrelationsVsQ3[1][2][1]->Fill(dQ3,agtrack3->Py()); // <X3>_y
      singleEventAverageCorrelationsVsQ3[1][3][1]->Fill(dQ3,agtrack1->Py()*agtrack2->Py()); // <X1X2>_y
      singleEventAverageCorrelationsVsQ3[1][4][1]->Fill(dQ3,agtrack1->Py()*agtrack3->Py()); // <X1X3>_y
      singleEventAverageCorrelationsVsQ3[1][5][1]->Fill(dQ3,agtrack2->Py()*agtrack3->Py()); // <X2X3>_y
      singleEventAverageCorrelationsVsQ3[1][6][1]->Fill(dQ3,agtrack1->Py()*agtrack2->Py()*agtrack3->Py()); // <X1X2X3>_y

      // z-components:
      singleEventAverageCorrelationsVsQ3[1][0][2]->Fill(dQ3,agtrack1->Pz()); // <X1>_z
      singleEventAverageCorrelationsVsQ3[1][1][2]->Fill(dQ3,agtrack2->Pz()); // <X2>_z
      singleEventAverageCorrelationsVsQ3[1][2][2]->Fill(dQ3,agtrack3->Pz()); // <X3>_z
      singleEventAverageCorrelationsVsQ3[1][3][2]->Fill(dQ3,agtrack1->Pz()*agtrack2->Pz()); // <X1X2>_z
      singleEventAverageCorrelationsVsQ3[1][4][2]->Fill(dQ3,agtrack1->Pz()*agtrack3->Pz()); // <X1X3>_z
      singleEventAverageCorrelationsVsQ3[1][5][2]->Fill(dQ3,agtrack2->Pz()*agtrack3->Pz()); // <X2X3>_z
      singleEventAverageCorrelationsVsQ3[1][6][2]->Fill(dQ3,agtrack1->Pz()*agtrack2->Pz()*agtrack3->Pz()); // <X1X2X3>_z

     } // if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE)))
    } // if(fFillCorrelationFunctionsTEST[1])


    // Test 2: "Same charge pions-kaons-protons, 2p correlations and cumulants projected onto Lorentz invariant Q3, using yield, and not momenta like in 'Test 0' and 'Test 1'"
    //         the first particle in the loop is always pion, the second is kaon, the third is proton. [xyz] is always se to 0
    if(fFillCorrelationFunctionsTEST[2])
    {
     if((Pion(gtrack1,1,kTRUE) && Kaon(gtrack2,1,kTRUE) && Proton(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Kaon(gtrack2,-1,kTRUE) && Proton(gtrack3,-1,kTRUE)))
     //if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE)))
     {
      // Count the yield differentially:
      Double_t dQ3 = Q3(agtrack1,agtrack2,agtrack3); // This is Lorentz invariant
      Int_t binNoQ3 = this->BinNoForSpecifiedValue(singleEventAverageCorrelationsVsQ3[2][0][0],dQ3); // TBI cross-check with some other profile
      if(binNoQ3>0) // otherwise it's either underflow or overflow
      {
       TString pattern1 = Form(":track%d:",iTrack1);
       TString pattern2 = Form(":track%d:",iTrack2);
       TString pattern3 = Form(":track%d:",iTrack3);
       if(!dycVsQ3[0][binNoQ3]->Contains(pattern1.Data()))
       {
        *dycVsQ3[0][binNoQ3] = Form("%s%s",dycVsQ3[0][binNoQ3]->Data(),pattern1.Data());
       }
       if(!dycVsQ3[1][binNoQ3]->Contains(pattern2.Data()))
       {
        *dycVsQ3[1][binNoQ3] = Form("%s%s",dycVsQ3[1][binNoQ3]->Data(),pattern2.Data());
       }
       if(!dycVsQ3[2][binNoQ3]->Contains(pattern3.Data()))
       {
        *dycVsQ3[2][binNoQ3] = Form("%s%s",dycVsQ3[2][binNoQ3]->Data(),pattern3.Data());
       }
      } // if(binNoQ3>0)

     } // if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE)))
    } // if(fFillCorrelationFunctionsTEST[2])

    // Test 4: "Same charge pions, 3p correlations and cumulants projected onto Lorentz invariant Q3, where energy prefactors have been taken into account, to make correlations and cumulants manifestly Lorenz invariant, even before projection"
    if(fFillCorrelationFunctionsTEST[4])
    {
     if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE)))
     {
      Double_t E1 = agtrack1->E(); // energy of 1st track TBI check the mass hypothesis
      Double_t E2 = agtrack2->E(); // energy of 2nd track TBI check the mass hypothesis
      Double_t E3 = agtrack3->E(); // energy of 2nd track TBI check the mass hypothesis
      Double_t dQ3 = Q3(agtrack1,agtrack2,agtrack3); // Lorentz invariant Q3

      singleEventAverageCorrelationsVsQ3[4][0][0]->Fill(dQ3,E1); // <X1>_x
      singleEventAverageCorrelationsVsQ3[4][1][0]->Fill(dQ3,E2); // <X2>_x
      singleEventAverageCorrelationsVsQ3[4][2][0]->Fill(dQ3,E3); // <X3>_x
      singleEventAverageCorrelationsVsQ3[4][3][0]->Fill(dQ3,E1*E2); // <X1X2>_x
      singleEventAverageCorrelationsVsQ3[4][4][0]->Fill(dQ3,E1*E3); // <X1X3>_x
      singleEventAverageCorrelationsVsQ3[4][5][0]->Fill(dQ3,E2*E3); // <X2X3>_x
      singleEventAverageCorrelationsVsQ3[4][6][0]->Fill(dQ3,E1*E2*E3); // <X1X2X3>_x
     } // if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE)))

    } // if(fFillCorrelationFunctionsTEST[4]) 

    // Test 5: "Same charge pions, 2p correlations and cumulants projected onto Lorentz invariant Q2, where energy prefactors have been taken into account, to make correlations and cumulants manifestly Lorenz invariant, even before projection. Here I am boosting all particles to a new frame, to check whether the final results are Lorentz invariant:"
    if(fFillCorrelationFunctionsTEST[5])
    {
     // p_1:
     Double_t p1x = agtrack1->Px();
     Double_t p1y = agtrack1->Py();
     Double_t p1z = agtrack1->Pz();
     Double_t e1  = agtrack1->E();
     // p_2:
     Double_t p2x = agtrack2->Px();
     Double_t p2y = agtrack2->Py();
     Double_t p2z = agtrack2->Pz();
     Double_t e2  = agtrack2->E();  
     // p_3:
     Double_t p3x = agtrack3->Px();
     Double_t p3y = agtrack3->Py();
     Double_t p3z = agtrack3->Pz();
     Double_t e3  = agtrack3->E();  

     // Corresponding energy-momentum four-vectors:
     TLorentzVector track1(p1x,p1y,p1z,e1);
     TLorentzVector track2(p2x,p2y,p2z,e2);
     TLorentzVector track3(p3x,p3y,p3z,e3);

     // Boosting both energy-momentum four-vectors:
     // N.B. I have a flag fBoost, for future cases...
     track1.Boost(fBoostVelocity);
     track2.Boost(fBoostVelocity);
     track3.Boost(fBoostVelocity);

     Double_t E1 = track1.E(); // energy of 1st track TBI check the mass hypothesis
     Double_t E2 = track2.E(); // energy of 2nd track TBI check the mass hypothesis
     Double_t E3 = track3.E(); // energy of 3rd track TBI check the mass hypothesis

     Double_t dQ3 = Q3(agtrack1,agtrack2,agtrack3); // Lorentz invariant Q3

     singleEventAverageCorrelationsVsQ3[5][0][0]->Fill(dQ3,E1); // <X1>_x
     singleEventAverageCorrelationsVsQ3[5][1][0]->Fill(dQ3,E2); // <X2>_x
     singleEventAverageCorrelationsVsQ3[5][2][0]->Fill(dQ3,E3); // <X3>_x
     singleEventAverageCorrelationsVsQ3[5][3][0]->Fill(dQ3,E1*E2); // <X1X2>_x
     singleEventAverageCorrelationsVsQ3[5][4][0]->Fill(dQ3,E1*E3); // <X1X3>_x
     singleEventAverageCorrelationsVsQ3[5][5][0]->Fill(dQ3,E2*E3); // <X2X3>_x
     singleEventAverageCorrelationsVsQ3[5][6][0]->Fill(dQ3,E1*E2*E3); // <X1X2X3>_x

    } // if(fFillCorrelationFunctionsTEST[5])

   } // for(Int_t iTrack3=0;iTrack3<nTracks;iTrack3++)
  } // for(Int_t iTrack2=0;iTrack2<nTracks;iTrack2++)
 } // for(Int_t iTrack1=0;iTrack1<nTracks;iTrack1++)


 // e) Build all-event correlations and cumulants from single-event averages:
 for(Int_t t=0;t<nTestsMax;t++) // test No
 {
  if(!fFillCorrelationFunctionsTEST[t]){continue;}
  Int_t nBins = singleEventAverageCorrelationsVsQ2[t][0][0]->GetXaxis()->GetNbins();
  for(Int_t b=0;b<nBins;b++) // TBI at the moment, I use the same binning for Q2 and Q3
  {
   for(Int_t xyz=0;xyz<3;xyz++)
   {
    Double_t dX1 = 0., dX1Err = 0.; // <X1>
    Double_t dX2 = 0., dX2Err = 0.; // <X2>
    Double_t dX3 = 0., dX3Err = 0.; // <X3>
    Double_t dX1X2 = 0., dX1X2Err = 0.; // <X1X2>
    Double_t dX1X3 = 0., dX1X3Err = 0.; // <X1X3>
    Double_t dX2X3 = 0., dX2X3Err = 0.; // <X2X3>
    Double_t dX1X2X3 = 0., dX1X2X3Err = 0.; // <X1X2X3>
    Double_t dC2 = 0., dC2Err = 0.; // C2 = <X1X2> - <X1><X2> vs. Q2
    Double_t dC12 = 0., dC12Err = 0.; // C12 = <X1X2> - <X1><X2> vs. Q3
    Double_t dC13 = 0., dC13Err = 0.; // C13 = <X1X3> - <X1><X3> vs. Q3
    Double_t dC23 = 0., dC23Err = 0.; // C23 = <X2X3> - <X2><X3> vs. Q3
    Double_t dC3 = 0., dC3Err = 0.; // C3 = <X1X2X3> - <X1X2><X3> - <X1X3><X2> - <X2X3><X1> + 2<X1><X2><X3> vs. Q3
    Double_t dBinCenter = 0.; // TBI assuming binning is everywhere the same. TBI it is used also below when filling 3p
    // 2p correlation terms:
    dX1 = singleEventAverageCorrelationsVsQ2[t][0][xyz]->GetBinContent(b+1); // <X1> vs. Q2
    dX1Err = singleEventAverageCorrelationsVsQ2[t][0][xyz]->GetBinError(b+1);
    dX2 = singleEventAverageCorrelationsVsQ2[t][1][xyz]->GetBinContent(b+1); // <X2> vs. Q2
    dX2Err = singleEventAverageCorrelationsVsQ2[t][1][xyz]->GetBinError(b+1);
    dX1X2 = singleEventAverageCorrelationsVsQ2[t][2][xyz]->GetBinContent(b+1); // <X1X2> vs. Q2
    dX1X2Err = singleEventAverageCorrelationsVsQ2[t][2][xyz]->GetBinError(b+1);
    if(TMath::Abs(dX1)>1.e-14 && TMath::Abs(dX2)>1.e-14 && TMath::Abs(dX1X2)>1.e-14) // basically, do calculations only if all terms in the 2p cumulant definition are non-zero
    {
     // 2p cumulant:
     dC2 = dX1X2 - dX1*dX2;
     dC2Err = 0.; // TBI propagate an error one day
     // fill 2p:
     dBinCenter = fCorrelationFunctionsTEST[t][0][0][xyz]->GetBinCenter(b+1); // TBI assuming binning is everywhere the same. TBI it is used also below when filling 3p
     fCorrelationFunctionsTEST[t][0][0][xyz]->Fill(dBinCenter,dX1); // <X1>
     fCorrelationFunctionsTEST[t][0][1][xyz]->Fill(dBinCenter,dX2); // <X2>
     fCorrelationFunctionsTEST[t][0][2][xyz]->Fill(dBinCenter,dX1X2); // <X1X2>
     fSignalCumulantsTEST[t][0][0][xyz]->Fill(dBinCenter,dC2); // C2 = <X1X2> - <X1><X2>
    }

    // TBI this was needed for 'Test 6', commenting out now temporarily:
    //dBinCenter = fCorrelationFunctionsTEST[t][0][0][xyz]->GetBinCenter(b+1); // TBI assuming binning is everywhere the same. TBI it is used also below when filling 3p
    //fCorrelationFunctionsTEST[t][0][0][xyz]->Fill(dBinCenter,dX1); // <X1>

    if(!fFill3pCorrelationFunctions){continue;} // TBI is this really safe

    // 3p correlation terms (note that now binning is vs. Q3 TBI yes, fine for the time being, but this is a landmine clearly...):
    dX1 = singleEventAverageCorrelationsVsQ3[t][0][xyz]->GetBinContent(b+1); // <X1> vs. Q3
    dX1Err = singleEventAverageCorrelationsVsQ3[t][0][xyz]->GetBinError(b+1);
    dX2 = singleEventAverageCorrelationsVsQ3[t][1][xyz]->GetBinContent(b+1); // <X2> vs. Q3
    dX2Err = singleEventAverageCorrelationsVsQ3[t][1][xyz]->GetBinError(b+1);
    dX3 = singleEventAverageCorrelationsVsQ3[t][2][xyz]->GetBinContent(b+1); // <X3> vs. Q3
    dX3Err = singleEventAverageCorrelationsVsQ3[t][2][xyz]->GetBinError(b+1);
    dX1X2 = singleEventAverageCorrelationsVsQ3[t][3][xyz]->GetBinContent(b+1); // <X1X2> vs. Q3
    dX1X2Err = singleEventAverageCorrelationsVsQ3[t][3][xyz]->GetBinError(b+1);
    dX1X3 = singleEventAverageCorrelationsVsQ3[t][4][xyz]->GetBinContent(b+1); // <X1X3> vs. Q3
    dX1X3Err = singleEventAverageCorrelationsVsQ3[t][4][xyz]->GetBinError(b+1);
    dX2X3 = singleEventAverageCorrelationsVsQ3[t][5][xyz]->GetBinContent(b+1); // <X2X3> vs. Q3
    dX2X3Err = singleEventAverageCorrelationsVsQ3[t][5][xyz]->GetBinError(b+1);
    dX1X2X3 = singleEventAverageCorrelationsVsQ3[t][6][xyz]->GetBinContent(b+1); // <X1X2X3> vs. Q3
    dX1X2X3Err = singleEventAverageCorrelationsVsQ3[t][6][xyz]->GetBinError(b+1);
    if(TMath::Abs(dX1)>1.e-14 && TMath::Abs(dX2)>1.e-14 && TMath::Abs(dX3)>1.e-14 &&
       TMath::Abs(dX1X2)>1.e-14 && TMath::Abs(dX1X3)>1.e-14 && TMath::Abs(dX2X3)>1.e-14 &&
       TMath::Abs(dX1X2X3)>1.e-14) // basically, do calculations only if all terms in the 3p cumulant definition are non-zero
    {
     // three 2p cumulants vs. Q3 (why not!?)
     dC12 = dX1X2 - dX1*dX2; // C12 = <X1X2> - <X1><X2> vs. Q3
     dC12Err = 0.; // TBI propagate an error one day
     dC13 = dX1X3 - dX1*dX3; // C13 = <X1X3> - <X1><X3> vs. Q3
     dC13Err = 0.; // TBI propagate an error one day
     dC23 = dX2X3 - dX2*dX3; // C23 = <X2X3> - <X2><X3> vs. Q3
     dC23Err = 0.; // TBI propagate an error one day
     // 3p cumulant vs. Q3:
     dC3 = dX1X2X3 - dX1X2*dX3 - dX1X3*dX2 - dX2X3*dX1 + 2.*dX1*dX2*dX3;
     dC3Err = 0.; // TBI propagate an error one day
     // fill 3p:
     dBinCenter = fCorrelationFunctionsTEST[t][1][0][xyz]->GetBinCenter(b+1); // TBI assuming binning is everywhere the same. TBI it is used also below when filling 3p
     fCorrelationFunctionsTEST[t][1][0][xyz]->Fill(dBinCenter,dX1); // <X1>
     fCorrelationFunctionsTEST[t][1][1][xyz]->Fill(dBinCenter,dX2); // <X2>
     fCorrelationFunctionsTEST[t][1][2][xyz]->Fill(dBinCenter,dX3); // <X3>
     fCorrelationFunctionsTEST[t][1][3][xyz]->Fill(dBinCenter,dX1X2); // <X1X2>
     fCorrelationFunctionsTEST[t][1][4][xyz]->Fill(dBinCenter,dX1X3); // <X1X3>
     fCorrelationFunctionsTEST[t][1][5][xyz]->Fill(dBinCenter,dX2X3); // <X2X3>
     fCorrelationFunctionsTEST[t][1][6][xyz]->Fill(dBinCenter,dX1X2X3); // <X1X2X3>
     fSignalCumulantsTEST[t][1][0][xyz]->Fill(dBinCenter,dC12); // <X1X2> - <X1><X2> vs. Q3
     fSignalCumulantsTEST[t][1][1][xyz]->Fill(dBinCenter,dC13); // <X1X3> - <X1><X3> vs. Q3
     fSignalCumulantsTEST[t][1][2][xyz]->Fill(dBinCenter,dC23); // <X2X3> - <X2><X3> vs. Q3
     fSignalCumulantsTEST[t][1][3][xyz]->Fill(dBinCenter,dC3); // C3 = <X1X2X3> - <X1X2><X3> - <X1X3><X2> - <X2X3><X1> + 2<X1><X2><X3> vs. Q3
    } // if(TMath::Abs(dX1)>1.e-14 && TMath::Abs(dX2)>1.e-14 && TMath::Abs(dX3)>1.e-14 && ...
   } // for(Int_t xyz=0;xyz<3;xyz++)
  } // for(Int_t b=0;b<nBins;b++) // TBI at the moment, I use same binning for Q2 and Q3


  if(2==t) // for "test 2" I just need correlations, and build cumulants at the end of the day. [x = pion yield][y = kaon yield][z = proton yield]
  {
   for(Int_t b=1;b<=nBins;b++) // TBI at the moment, I use same binning for Q2 and Q3. Note that here I start a loop over 1, not like in previous two tests
   {
    Double_t dBinCenter = fCorrelationFunctionsTEST[t][1][0][0]->GetBinCenter(b); // TBI assuming binning is everywhere the same. TBI it is used also below when filling 3p

    Int_t counter[3] = {0}; // [0=pions][1=kaons][2=protons]
    for(Int_t pkp=0;pkp<3;pkp++)
    {
     TObjArray *oa = TString(dycVsQ3[pkp][b]->Data()).Tokenize(":");
     while(oa->At(counter[pkp]))
     {
      counter[pkp]++;
     }
    }

    Int_t dX1 = counter[0]; // X1 (not the average!)
    Int_t dX2 = counter[1]; // X2 (not the average!)
    Int_t dX3 = counter[2]; // X3 (not the average!)
    Int_t dX1X2 = counter[0]*counter[1]; // X1X2
    Int_t dX1X3 = counter[0]*counter[2]; // X1X3
    Int_t dX2X3 = counter[1]*counter[2]; // X2X3
    Int_t dX1X2X3 = counter[0]*counter[1]*counter[2]; // X1X2X3
    fCorrelationFunctionsTEST[2][1][0][0]->Fill(dBinCenter,dX1); // <X1>
    fCorrelationFunctionsTEST[2][1][1][0]->Fill(dBinCenter,dX2); // <X2>
    fCorrelationFunctionsTEST[2][1][2][0]->Fill(dBinCenter,dX3); // <X3>
    fCorrelationFunctionsTEST[2][1][3][0]->Fill(dBinCenter,dX1X2); // <X1X2>
    fCorrelationFunctionsTEST[2][1][4][0]->Fill(dBinCenter,dX1X3); // <X1X3>
    fCorrelationFunctionsTEST[2][1][5][0]->Fill(dBinCenter,dX2X3); // <X2X3>
    fCorrelationFunctionsTEST[2][1][6][0]->Fill(dBinCenter,dX1X2X3); // <X1X2X3>

   } // for(Int_t b=1;b<=nBins;b++)

  } // if(2==t) // for "test 2" I just need correlations, and build cumulants at the end of the day

  if(7==t) // for "test 7" special treatment
  {
 
   for(Int_t b=1;b<=fnQ2bins;b++) // TBI check the boundaries
   {
    Int_t counter[2] = {0,0}; // [0=pions+][1=pions-]

    if(test7Counters[b][0] && !test7Counters[b][0]->EqualTo(""))
    {
     // Count pions+
     TObjArray *oa0 = TString(test7Counters[b][0]->Data()).Tokenize("_");
     while(oa0->At(counter[0]))
     {
      counter[0]++;
     }
     //cout<<Form("Bin %d; N = %d; STRING1: ",b,counter[0])<<test7Counters[b][0]->Data()<<endl;
    } // if(test7Counters[b][0] && !test7Counters[b][0]->EqualTo(""))
  
    if(test7Counters[b][1] && !test7Counters[b][1]->EqualTo(""))
    {
     // Count pions-
     TObjArray *oa1 = TString(test7Counters[b][1]->Data()).Tokenize("_");
     while(oa1->At(counter[1]))
     {
      counter[1]++;
     }
     //cout<<Form("Bin %d; N = %d; STRING2: ",b,counter[1])<<test7Counters[b][1]->Data()<<endl;
    } // if(test7Counters[b][1] && !test7Counters[b][1]->EqualTo("")) 

    Int_t nPionsP = counter[0]; // number of positive pions in bth Q2 bin, in this event, X1 in the cumulant formalism
    Int_t nPionsN = counter[1]; // number of negative pions in bth Q2 bin, in this event, X2 in the cumulant formalism

    if(nPionsP>=1 && nPionsN>=1)
    {
     Double_t dBinCenter = fCorrelationFunctionsTEST[7][0][0][0]->GetBinCenter(b); // TBI assuming same binning everyrhere
     //TProfile *fCorrelationFunctionsTEST[10][2][7][10]; //! [testNo][0=vs Q2, 1=vs Q3][example [0=<x1>][1=<x2>], ...,[6=<x1x2x3>]][differential index, e.g. for test 0 [0=Cx][1=Cy][2=Cz]]
     fCorrelationFunctionsTEST[7][0][0][0]->Fill(dBinCenter,nPionsP); // X1 
     fCorrelationFunctionsTEST[7][0][1][0]->Fill(dBinCenter,nPionsN); // X2
     fCorrelationFunctionsTEST[7][0][2][0]->Fill(dBinCenter,nPionsP*nPionsN); // X1*X2 
    } // if(nPionsP>=1 && nPionsN>=1)
 
   } // for(Int_t b=1;b<=fnQ2bins;b++) // TBI check the boundaries

  } // if(7==t) // for "test 7" special treatment

 } // for(Int_t t=0;t<nTestsMax;t++) // test No


 // f) Delete local TProfile's to hold single-event averages:
 for(Int_t t=0;t<nTestsMax;t++) // test No
 {
  if(!fFillCorrelationFunctionsTEST[t]){continue;}
  for(Int_t ct=0;ct<n2pCumulantTerms;ct++)
  {
   for(Int_t xyz=0;xyz<3;xyz++)
   {
    delete singleEventAverageCorrelationsVsQ2[t][ct][xyz];
   }
  }
 }

 for(Int_t t=0;t<nTestsMax;t++) // test No
 {
  for(Int_t ct=0;ct<n3pCumulantTerms;ct++)
  {
   for(Int_t xyz=0;xyz<3;xyz++)
   {
    if(!fFillCorrelationFunctionsTEST[t]){continue;}
    delete singleEventAverageCorrelationsVsQ3[t][ct][xyz];
   }
  }
 }

 for(Int_t nl=0;nl<nLoopsVsQ2;nl++)
 {
  for(Int_t nme=1;nme<=nMaxEntries;nme++) // yes, from 1, see above
  {
   if(nme > fnQ2bins){break;}
   delete dycVsQ2[nl][nme];
  }
 }
 for(Int_t nl=0;nl<nLoopsVsQ3;nl++)
 {
  for(Int_t nme=1;nme<=nMaxEntries;nme++) // yes, from 1, see above
  {
   if(nme > fnQ3bins){break;}
   delete dycVsQ3[nl][nme];
  }
 }

 for(Int_t q2=0;q2<fnQ2bins;q2++)
 {
  if(test7Counters[q2][0]) delete test7Counters[q2][0];
  if(test7Counters[q2][1]) delete test7Counters[q2][1];
 }

 // c) Three nested loops to calculate C(Q3), just an example; TBI
 //if(fFill3pCorrelationFunctions) this->Calculate3pCorrelationFunctions(aAOD);

 // d) Four nested loops to calculate C(Q4), just an example. TBI
 //if(fFill4pCorrelationFunctions) this->Calculate4pCorrelationFunctions(aAOD);

} // void AliAnalysisTaskMultiparticleFemtoscopy::CalculateCorrelationFunctionsTEST(AliAODEvent *aAOD)

//=======================================================================================================================

Double_t AliAnalysisTaskMultiparticleFemtoscopy::RelativeMomenta(AliAODTrack *agtrack1, AliAODTrack *agtrack2)
{
 // Comment the weather here. TBI

 Double_t k = -44.; // relative momenta k = \frac{1}{2}|\vec{p_1}-\vec{p_2}|

 // \vec{p_1}:
 Double_t p1x = agtrack1->Px();
 Double_t p1y = agtrack1->Py();
 Double_t p1z = agtrack1->Pz();
 // \vec{p_2}:
 Double_t p2x = agtrack2->Px();
 Double_t p2y = agtrack2->Py();
 Double_t p2z = agtrack2->Pz();
 // k:
 k = 0.5*pow(pow(p1x-p2x,2.)+pow(p1y-p2y,2.)+pow(p1z-p2z,2.),0.5);

 if(k<1.e-14) // TBI rething this constraint
 {
  cout<<Form("\nSelf-correlation !!!!")<<endl;
  cout<<Form("p1: %f %f %f",p1x,p1y,p1z)<<endl;
  cout<<Form("p2: %f %f %f\n",p2x,p2y,p2z)<<endl;
  exit(0);
 }

 return k;

} // Double_t AliAnalysisTaskMultiparticleFemtoscopy::RelativeMomenta(AliAODTrack *agtrack1, AliAODTrack *agtrack2)

//=======================================================================================================================

Double_t AliAnalysisTaskMultiparticleFemtoscopy::RelativeMomentaComponent(AliAODTrack *agtrack1, AliAODTrack *agtrack2, const char *component)
{
 // Return the component of relative momenta.

 // Supported:
 // a) "x"|"y"|"z" => Cartesian x,y,z
 // b) TBI add support for other coordinate systems, e.g. o-s-l, spherical, etc.

 TString sMethodName = "Double_t AliAnalysisTaskMultiparticleFemtoscopy::RelativeMomentaComponent(AliAODTrack *agtrack1, AliAODTrack *agtrack2, const char *component)";
 if(!(TString(component).EqualTo("x")||TString(component).EqualTo("y")||TString(component).EqualTo("z"))){Fatal(sMethodName.Data(),"This value for \"component\" is not supported yet: %s",component);}

 Double_t k = -44.; // relative momenta component

 if(TString(component).EqualTo("x"))
 {
  k = agtrack1->Px() - agtrack2->Px();
 }
 else if(TString(component).EqualTo("y"))
 {
  k = agtrack1->Py() - agtrack2->Py();
 }
 else if(TString(component).EqualTo("z"))
 {
  k = agtrack1->Pz() - agtrack2->Pz();
 }

 return k;

} // Double_t AliAnalysisTaskMultiparticleFemtoscopy::RelativeMomentaComponent(AliAODTrack *agtrack1, AliAODTrack *agtrack2, const char *component)


//=======================================================================================================================

Double_t AliAnalysisTaskMultiparticleFemtoscopy::PairVectorComponent(AliAODTrack *agtrack1, AliAODTrack *agtrack2, const char *component)
{
 // Return the component of pair vector \vec{K} = 1/2 (\vec{p1} + \vec{p2})

 // Supported:
 // a) "x"|"y"|"z" => Cartesian x,y,z
 // b) TBI add support for other coordinate systems, e.g. o-s-l, spherical, etc.

 TString sMethodName = "Double_t AliAnalysisTaskMultiparticleFemtoscopy::RelativeMomentaComponent(AliAODTrack *agtrack1, AliAODTrack *agtrack2, const char *component)";
 if(!(TString(component).EqualTo("x")||TString(component).EqualTo("y")||TString(component).EqualTo("z"))){Fatal(sMethodName.Data(),"This value for \"component\" is not supported yet: %s",component);}

 Double_t k = -44.; // pair vector component

 if(TString(component).EqualTo("x"))
 {
  k = 0.5*(agtrack1->Px() + agtrack2->Px());
 }
 else if(TString(component).EqualTo("y"))
 {
  k = 0.5*(agtrack1->Py() + agtrack2->Py());
 }
 else if(TString(component).EqualTo("z"))
 {
  k = 0.5*(agtrack1->Pz() + agtrack2->Pz());
 }

 return k;

} // Double_t AliAnalysisTaskMultiparticleFemtoscopy::PairVectorComponent(AliAODTrack *agtrack1, AliAODTrack *agtrack2, const char *component)

//=======================================================================================================================

Double_t AliAnalysisTaskMultiparticleFemtoscopy::RelativeMomenta(AliAODMCParticle *amcparticle1, AliAODMCParticle *amcparticle2)
{
 // Comment the weather here. TBI

 Double_t k = -44.; // relative momenta k = \frac{1}{2}|\vec{p_1}-\vec{p_2}|

 // \vec{p_1}:
 Double_t p1x = amcparticle1->Px();
 Double_t p1y = amcparticle1->Py();
 Double_t p1z = amcparticle1->Pz();
 // \vec{p_2}:
 Double_t p2x = amcparticle2->Px();
 Double_t p2y = amcparticle2->Py();
 Double_t p2z = amcparticle2->Pz();
 // k:
 k = 0.5*pow(pow(p1x-p2x,2.)+pow(p1y-p2y,2.)+pow(p1z-p2z,2.),0.5);

 if(k<1.e-14) // TBI rething this constraint
 {
  cout<<Form("\nSelf-correlation !!!!")<<endl;
  cout<<Form("p1: %f %f %f",p1x,p1y,p1z)<<endl;
  cout<<Form("p2: %f %f %f\n",p2x,p2y,p2z)<<endl;
  exit(0);
 }

 return k;

} // Double_t AliAnalysisTaskMultiparticleFemtoscopy::RelativeMomenta(AliAODMCParticle *amcparticle1, AliAODMCParticle *amcparticle2)

//=======================================================================================================================

Double_t AliAnalysisTaskMultiparticleFemtoscopy::Q2(AliAODTrack *agtrack1, AliAODTrack *agtrack2)
{
 // Lorentz invariant Q2. Credits: Sir Oliver.

 // p_1:
 Double_t p1x = agtrack1->Px();
 Double_t p1y = agtrack1->Py();
 Double_t p1z = agtrack1->Pz();
 Double_t e1  = agtrack1->E();

 // p_2:
 Double_t p2x = agtrack2->Px();
 Double_t p2y = agtrack2->Py();
 Double_t p2z = agtrack2->Pz();
 Double_t e2  = agtrack2->E();

 // Corresponding energy-momentum four-vectors:
 TLorentzVector track1(p1x,p1y,p1z,e1);
 TLorentzVector track2(p2x,p2y,p2z,e2);

 // Standard gym. to get k*:
 TLorentzVector trackSum = track1 + track2;
 Double_t beta = trackSum.Beta();
 Double_t beta_x = beta*cos(trackSum.Phi())*sin(trackSum.Theta());
 Double_t beta_y = beta*sin(trackSum.Phi())*sin(trackSum.Theta());
 Double_t beta_z = beta*cos(trackSum.Theta());
 TLorentzVector track1_cms = track1;
 TLorentzVector track2_cms = track2;
 track1_cms.Boost(-beta_x,-beta_y,-beta_z);
 track2_cms.Boost(-beta_x,-beta_y,-beta_z);

 TLorentzVector track_relK = track1_cms - track2_cms;
 Double_t Q2 = track_relK.P();

 return Q2;

} // Double_t AliAnalysisTaskMultiparticleFemtoscopy::Q2(AliAODTrack *agtrack1, AliAODTrack *agtrack2)

//=======================================================================================================================

Double_t AliAnalysisTaskMultiparticleFemtoscopy::Q3(AliAODTrack *agtrack1, AliAODTrack *agtrack2, AliAODTrack *agtrack3)
{
 // The Lorentz invariant Q3. TBI validate with some explicit calculation

 Double_t q12 = Q2(agtrack1,agtrack2);
 Double_t q13 = Q2(agtrack1,agtrack3);
 Double_t q23 = Q2(agtrack2,agtrack3);
 Double_t Q3 = pow(pow(q12,2.)+pow(q13,2.)+pow(q23,2.),0.5);

 return Q3;

} // Double_t AliAnalysisTaskMultiparticleFemtoscopy::Q3(AliAODTrack *agtrack1, AliAODTrack *agtrack2, AliAODTrack *agtrack3)

//=======================================================================================================================

Double_t AliAnalysisTaskMultiparticleFemtoscopy::Q4(AliAODTrack *agtrack1, AliAODTrack *agtrack2, AliAODTrack *agtrack3, AliAODTrack *agtrack4)
{
 // The Lorentz invariant Q4. TBI validate with some explicit calculation

 Double_t q12 = Q2(agtrack1,agtrack2);
 Double_t q13 = Q2(agtrack1,agtrack3);
 Double_t q14 = Q2(agtrack1,agtrack4);
 Double_t q23 = Q2(agtrack2,agtrack3);
 Double_t q24 = Q2(agtrack2,agtrack4);
 Double_t q34 = Q2(agtrack3,agtrack4);

 Double_t Q4 = pow(pow(q12,2.)+pow(q13,2.)+pow(q14,2.)+pow(q23,2.)+pow(q24,2.)+pow(q34,2.),0.5);

 return Q4;

} // Double_t AliAnalysisTaskMultiparticleFemtoscopy::Q4(AliAODTrack *agtrack1, AliAODTrack *agtrack2, AliAODTrack *agtrack3, AliAODTrack *agtrack4)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::Calculate2pBackground(TClonesArray *ca1, TClonesArray *ca2)
{
 // Calculate background.

 // TBI this method is not really validated

 // a) Insanity checks;
 // b) Temporary primitive algorithm: correlate particles from 'previous event' to particles from 'current event';
 // c) Two nested loops to calculate B(k);
 // d) Shift [1] -> [0]; TBI I have moved this on another place
 // e) Clean [1]. TBI I have moved this on another place

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::Calculate2pBackground(TClonesArray *ca1, TClonesArray *ca2)";
 if(!ca1){Fatal(sMethodName.Data(),"!ca1");}
 if(!ca2){Fatal(sMethodName.Data(),"!ca2");}

 // b) Temporary primitive algorithm: correlate particles from 'previous event' to particles from 'current event'
 // ...

 // c) Two nested loops to calculate B(k):
 Int_t nTracks1 = ca1->GetEntries();
 Int_t nTracks2 = ca2->GetEntries();

 // Start loop over tracks in the 1st event:
 for(Int_t iTrack1=0;iTrack1<nTracks1;iTrack1++)
 {
  AliAODTrack *atrack1 = dynamic_cast<AliAODTrack*>(ca1->UncheckedAt(iTrack1)); // TBI cross-check UncheckedAt
  // TBI Temporary track insanity checks:
  if(!atrack1){Fatal(sMethodName.Data(),"!atrack1");} // TBI keep this for some time, eventually just continue
  if(atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(atrack1->TestFilterBit(128) && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->TestFiletrBit(128) && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(!PassesCommonTrackCuts(atrack1)){continue;} // TBI re-think
  // Corresponding AOD global track:
  Int_t id1 = atrack1->GetID();
  AliAODTrack *gtrack1 = dynamic_cast<AliAODTrack*>(id1>=0 ? ca1->UncheckedAt(fGlobalTracksAOD[1]->GetValue(id1)) : ca1->UncheckedAt(fGlobalTracksAOD[1]->GetValue(-(id1+1))));
  if(!gtrack1){Fatal(sMethodName.Data(),"!gtrack1");} // TBI keep this for some time, eventually just continue
  // Common track selection criteria for all "normal" global tracks:
  if(!PassesGlobalTrackCuts(gtrack1)){continue;}

  // Start loop over tracks in the 2nd event:
  for(Int_t iTrack2=0;iTrack2<nTracks2;iTrack2++)
  {
   AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(ca2->UncheckedAt(iTrack2));
   // TBI Temporary track insanity checks:
   if(!atrack2){Fatal(sMethodName.Data(),"!atrack2");} // TBI keep this for some time, eventually just continue
   if(atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(atrack2->TestFilterBit(128) && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->TestFiletrBit(128) && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(!PassesCommonTrackCuts(atrack2)){continue;} // TBI re-think
   // Corresponding AOD global track:
   Int_t id2 = atrack2->GetID();
   AliAODTrack *gtrack2 = dynamic_cast<AliAODTrack*>(id2>=0 ? ca2->UncheckedAt(fGlobalTracksAOD[2]->GetValue(id2)) : ca2->UncheckedAt(fGlobalTracksAOD[2]->GetValue(-(id2+1))));
   if(!gtrack2){Fatal(sMethodName.Data(),"!gtrack2");} // TBI keep this for some time, eventually just continue
   // Common track selection criteria for all "normal" global tracks:
   if(!PassesGlobalTrackCuts(gtrack2)){continue;}

   // Okay... So we have two tracks from two different events ready. let's decide whether kinemetics from 'atrack' or 'gtrack' is taken, let's check PID, and calculate the background:
   AliAODTrack *agtrack1 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
   AliAODTrack *agtrack2 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
   if(fFillControlHistogramsWithGlobalTrackInfo)
   {
    agtrack1 = gtrack1;
    agtrack2 = gtrack2;
   }
   else
   {
    agtrack1 = atrack1;
    agtrack2 = atrack2;
   } // if(fFillControlHistogramsWithGlobalTrackInfo)
   if(!agtrack1){Fatal(sMethodName.Data(),"!agtrack1");}
   if(!agtrack2){Fatal(sMethodName.Data(),"!agtrack2");}

   // 1.) Same particle species:

   // a) pion-pion:
   //  a1) pi+pi+ [2][2]:
   if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE))
   {
    f2pBackground[2][2]->Fill(Q2(agtrack1,agtrack2));
   }
   //  a2) pi-pi- [7][7]:
   if(Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE))
   {
    f2pBackground[7][7]->Fill(Q2(agtrack1,agtrack2));
   }
   //  a3) pi+pi- || pi-pi+ [2][7]:
   if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,1,kTRUE)))
   {
    f2pBackground[2][7]->Fill(Q2(agtrack1,agtrack2));
   }

   // b) kaon-kaon:
   //  b1) K+K+ [3][3]:
   if(Kaon(gtrack1,1,kTRUE) && Kaon(gtrack2,1,kTRUE))
   {
    f2pBackground[3][3]->Fill(Q2(agtrack1,agtrack2));
   }
   //  b2) K-K- [8][8]:
   if(Kaon(gtrack1,-1,kTRUE) && Kaon(gtrack2,-1,kTRUE))
   {
    f2pBackground[8][8]->Fill(Q2(agtrack1,agtrack2));
   }
   //  b3) K+K- || K-K+ [3][8]:
   if((Kaon(gtrack1,1,kTRUE) && Kaon(gtrack2,-1,kTRUE)) || (Kaon(gtrack1,-1,kTRUE) && Kaon(gtrack2,1,kTRUE)))
   {
    f2pBackground[3][8]->Fill(Q2(agtrack1,agtrack2));
   }

   // c) proton-proton:
   //  c1) p+p+ [4][4]:
   if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE))
   {
    f2pBackground[4][4]->Fill(Q2(agtrack1,agtrack2));
   }
   //  c2) p-p- [9][9]:
   if(Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    f2pBackground[9][9]->Fill(Q2(agtrack1,agtrack2));
   }
   //  c3) p+p- || p-p+ [4][9]:
   if((Proton(gtrack1,1,kTRUE) && Proton(gtrack2,-1,kTRUE)) || (Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,1,kTRUE)))
   {
    f2pBackground[4][9]->Fill(Q2(agtrack1,agtrack2));
   }

   // 2.) Mixed particle species:
   // a) pion-kaon
   //  a1) pi+K+ [2][3]:
   if(Pion(gtrack1,1,kTRUE) && Kaon(gtrack2,1,kTRUE))
   {
    f2pBackground[2][3]->Fill(Q2(agtrack1,agtrack2));
   }
   //  a2) pi+K- [2][8]
   if(Pion(gtrack1,1,kTRUE) && Kaon(gtrack2,-1,kTRUE))
   {
    f2pBackground[2][8]->Fill(Q2(agtrack1,agtrack2));
   }
   //  a3) K+pi- [3][7]
   if(Kaon(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE))
   {
    f2pBackground[3][7]->Fill(Q2(agtrack1,agtrack2));
   }
   //  a4) pi-K- [7][8]
   if(Pion(gtrack1,-1,kTRUE) && Kaon(gtrack2,-1,kTRUE))
   {
    f2pBackground[7][8]->Fill(Q2(agtrack1,agtrack2));
   }
   // b) pion-proton
   //  b1) pi+p+ [2][4]:
   if(Pion(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE))
   {
    f2pBackground[2][4]->Fill(Q2(agtrack1,agtrack2));
   }
   //  b2) pi+p- [2][9]
   if(Pion(gtrack1,1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    f2pBackground[2][9]->Fill(Q2(agtrack1,agtrack2));
   }
   //  b3) p+pi- [4][7]
   if(Proton(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE))
   {
    f2pBackground[4][7]->Fill(Q2(agtrack1,agtrack2));
   }
   //  b4) pi-p- [7][9]
   if(Pion(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    f2pBackground[7][9]->Fill(Q2(agtrack1,agtrack2));
   }
   // c) kaon-proton
   //  c1) K+p+ [3][4]:
   if(Kaon(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE))
   {
    f2pBackground[3][4]->Fill(Q2(agtrack1,agtrack2));
   }
   //  c2) K+p- [3][9]
   if(Kaon(gtrack1,1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    f2pBackground[3][9]->Fill(Q2(agtrack1,agtrack2));
   }
   //  c3) p+K- [4][8]
   if(Proton(gtrack1,1,kTRUE) && Kaon(gtrack2,-1,kTRUE))
   {
    f2pBackground[4][8]->Fill(Q2(agtrack1,agtrack2));
   }
   //  c4) K-p- [8][9]
   if(Kaon(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    f2pBackground[8][9]->Fill(Q2(agtrack1,agtrack2));
   }

  } // for(Int_t iTrack2=0;iTrack2<nTracks2;iTrack2++)
 } // for(Int_t iTrack1=0;iTrack1<nTracks1;iTrack1++)

} // void AliAnalysisTaskMultiparticleFemtoscopy::Calculate2pBackground(TClonesArray *ca1, TClonesArray *ca2)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::Calculate2pBackground(TClonesArray *ca1, TClonesArray *ca2, TExMap *em1, TExMap *em2)
{
 // Calculate background.

 // a) Insanity checks;
 // b) Two nested loops to calculate B(k);
 // c) Shift [1] -> [0]; TBI I have moved this on another place
 // d) Clean [1]. TBI I have moved this on another place

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::Calculate2pBackground(TClonesArray *ca1, TClonesArray *ca2, TExMap *em1, TExMap *em2)";
 if(!ca1){Fatal(sMethodName.Data(),"!ca1");}
 if(!ca2){Fatal(sMethodName.Data(),"!ca2");}
 if(!em1){Fatal(sMethodName.Data(),"!em1");}
 if(!em2){Fatal(sMethodName.Data(),"!em2");}

 // b) Two nested loops to calculate B(k):
 Int_t nTracks1 = ca1->GetEntries();
 Int_t nTracks2 = ca2->GetEntries();

 // Start loop over tracks in the 1st event:
 for(Int_t iTrack1=0;iTrack1<nTracks1;iTrack1++)
 {
  AliAODTrack *atrack1 = dynamic_cast<AliAODTrack*>(ca1->UncheckedAt(iTrack1)); // TBI cross-check UncheckedAt
  // TBI Temporary track insanity checks:
  if(!atrack1){Fatal(sMethodName.Data(),"!atrack1");} // TBI keep this for some time, eventually just continue
  if(atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(atrack1->TestFilterBit(128) && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->TestFiletrBit(128) && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(!PassesCommonTrackCuts(atrack1)){continue;} // TBI re-think
  // Corresponding AOD global track:
  Int_t id1 = atrack1->GetID();
  AliAODTrack *gtrack1 = dynamic_cast<AliAODTrack*>(id1>=0 ? ca1->UncheckedAt(em1->GetValue(id1)) : ca1->UncheckedAt(em1->GetValue(-(id1+1))));
  if(!gtrack1){Fatal(sMethodName.Data(),"!gtrack1");} // TBI keep this for some time, eventually just continue
  // Common track selection criteria for all "normal" global tracks:
  if(!PassesGlobalTrackCuts(gtrack1)){continue;}

  // Start loop over tracks in the 2nd event:
  for(Int_t iTrack2=0;iTrack2<nTracks2;iTrack2++)
  {
   AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(ca2->UncheckedAt(iTrack2));
   // TBI Temporary track insanity checks:
   if(!atrack2){Fatal(sMethodName.Data(),"!atrack2");} // TBI keep this for some time, eventually just continue
   if(atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(atrack2->TestFilterBit(128) && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->TestFiletrBit(128) && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(!PassesCommonTrackCuts(atrack2)){continue;} // TBI re-think
   // Corresponding AOD global track:
   Int_t id2 = atrack2->GetID();
   AliAODTrack *gtrack2 = dynamic_cast<AliAODTrack*>(id2>=0 ? ca2->UncheckedAt(em2->GetValue(id2)) : ca2->UncheckedAt(em2->GetValue(-(id2+1))));
   if(!gtrack2){Fatal(sMethodName.Data(),"!gtrack2");} // TBI keep this for some time, eventually just continue
   // Common track selection criteria for all "normal" global tracks:
   if(!PassesGlobalTrackCuts(gtrack2)){continue;}

   // Okay... So we have two tracks from two different events ready. let's decide whether kinemetics from 'atrack' or 'gtrack' is taken, let's check PID, and calculate the background:
   AliAODTrack *agtrack1 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
   AliAODTrack *agtrack2 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
   if(fFillControlHistogramsWithGlobalTrackInfo)
   {
    agtrack1 = gtrack1;
    agtrack2 = gtrack2;
   }
   else
   {
    agtrack1 = atrack1;
    agtrack2 = atrack2;
   } // if(fFillControlHistogramsWithGlobalTrackInfo)
   if(!agtrack1){Fatal(sMethodName.Data(),"!agtrack1");}
   if(!agtrack2){Fatal(sMethodName.Data(),"!agtrack2");}

   // 1.) Same particle species:
   // a) pion-pion:
   //  a1) pi+pi+ [2][2]:
   if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE))
   {
    f2pBackground[2][2]->Fill(Q2(agtrack1,agtrack2));
   }

   //  a2) pi-pi- [7][7]:
   if(Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE))
   {
    f2pBackground[7][7]->Fill(Q2(agtrack1,agtrack2));
   }
   //  a3) pi+pi- || pi-pi+ [2][7]:
   if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,1,kTRUE)))
   {
    f2pBackground[2][7]->Fill(Q2(agtrack1,agtrack2));
   }

   // b) kaon-kaon:
   //  b1) K+K+ [3][3]:
   if(Kaon(gtrack1,1,kTRUE) && Kaon(gtrack2,1,kTRUE))
   {
    f2pBackground[3][3]->Fill(Q2(agtrack1,agtrack2));
   }
   //  b2) K-K- [8][8]:
   if(Kaon(gtrack1,-1,kTRUE) && Kaon(gtrack2,-1,kTRUE))
   {
    f2pBackground[8][8]->Fill(Q2(agtrack1,agtrack2));
   }
   //  b3) K+K- || K-K+ [3][8]:
   if((Kaon(gtrack1,1,kTRUE) && Kaon(gtrack2,-1,kTRUE)) || (Kaon(gtrack1,-1,kTRUE) && Kaon(gtrack2,1,kTRUE)))
   {
    f2pBackground[3][8]->Fill(Q2(agtrack1,agtrack2));
   }

   // c) proton-proton:
   //  c1) p+p+ [4][4]:
   if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE))
   {
    f2pBackground[4][4]->Fill(Q2(agtrack1,agtrack2));
   }
   //  c2) p-p- [9][9]:
   if(Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    f2pBackground[9][9]->Fill(Q2(agtrack1,agtrack2));
   }
   //  c3) p+p- || p-p+ [4][9]:
   if((Proton(gtrack1,1,kTRUE) && Proton(gtrack2,-1,kTRUE)) || (Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,1,kTRUE)))
   {
    f2pBackground[4][9]->Fill(Q2(agtrack1,agtrack2));
   }

   // 2.) Mixed particle species:
   // a) pion-kaon
   //  a1) pi+K+ [2][3]:
   if(Pion(gtrack1,1,kTRUE) && Kaon(gtrack2,1,kTRUE))
   {
    f2pBackground[2][3]->Fill(Q2(agtrack1,agtrack2));
   }
   //  a2) pi+K- [2][8]
   if(Pion(gtrack1,1,kTRUE) && Kaon(gtrack2,-1,kTRUE))
   {
    f2pBackground[2][8]->Fill(Q2(agtrack1,agtrack2));
   }
   //  a3) K+pi- [3][7]
   if(Kaon(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE))
   {
    f2pBackground[3][7]->Fill(Q2(agtrack1,agtrack2));
   }
   //  a4) pi-K- [7][8]
   if(Pion(gtrack1,-1,kTRUE) && Kaon(gtrack2,-1,kTRUE))
   {
    f2pBackground[7][8]->Fill(Q2(agtrack1,agtrack2));
   }
   // b) pion-proton
   //  b1) pi+p+ [2][4]:
   if(Pion(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE))
   {
    f2pBackground[2][4]->Fill(Q2(agtrack1,agtrack2));
   }
   //  b2) pi+p- [2][9]
   if(Pion(gtrack1,1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    f2pBackground[2][9]->Fill(Q2(agtrack1,agtrack2));
   }
   //  b3) p+pi- [4][7]
   if(Proton(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE))
   {
    f2pBackground[4][7]->Fill(Q2(agtrack1,agtrack2));
   }
   //  b4) pi-p- [7][9]
   if(Pion(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    f2pBackground[7][9]->Fill(Q2(agtrack1,agtrack2));
   }
   // c) kaon-proton
   //  c1) K+p+ [3][4]:
   if(Kaon(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE))
   {
    f2pBackground[3][4]->Fill(Q2(agtrack1,agtrack2));
   }
   //  c2) K+p- [3][9]
   if(Kaon(gtrack1,1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    f2pBackground[3][9]->Fill(Q2(agtrack1,agtrack2));
   }
   //  c3) p+K- [4][8]
   if(Proton(gtrack1,1,kTRUE) && Kaon(gtrack2,-1,kTRUE))
   {
    f2pBackground[4][8]->Fill(Q2(agtrack1,agtrack2));
   }
   //  c4) K-p- [8][9]
   if(Kaon(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    f2pBackground[8][9]->Fill(Q2(agtrack1,agtrack2));
   }

  } // for(Int_t iTrack2=0;iTrack2<nTracks2;iTrack2++)
 } // for(Int_t iTrack1=0;iTrack1<nTracks1;iTrack1++)

} // void AliAnalysisTaskMultiparticleFemtoscopy::Calculate2pBackground(TClonesArray *ca1, TClonesArray *ca2, TExMap *em1, TExMap *em2)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::Calculate2pBackgroundTEST(TClonesArray *ca1, TClonesArray *ca2, TExMap *em1, TExMap *em2)
{
 // Calculate 2p TEST background.

 // a) Insanity checks;
 // b) Book local TProfile's to hold 1x1 event averages vs. Q2;
 // c) Calculate averages <X1> (1st event), <X2> (2nd event) and <X1X2> (mixed-event), vs. Q2.
 // d) Build all-event correlations and cumulants from single-event averages.

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::Calculate2pBackgroundTEST(TClonesArray *ca1, TClonesArray *ca2, TExMap *em1, TExMap *em2)";
 if(!ca1){Fatal(sMethodName.Data(),"!ca1");}
 if(!ca2){Fatal(sMethodName.Data(),"!ca2");}
 if(!em1){Fatal(sMethodName.Data(),"!em1");}
 if(!em2){Fatal(sMethodName.Data(),"!em2");}

 // b) Book local TProfile's to hold 1x1 event averages vs. Q2:
 const Int_t nTestsMax = 10; // see what is hardwired in .h file
 TString sXYZ[3] = {"x","y","z"};
 const Int_t n2pCumulantTerms = 3;
 TString s2pCumulantTerms[n2pCumulantTerms] = {"#LTX_{1}#GT","#LTX_{2}#GT","#LTX_{1}X_{2}#GT"};
 TProfile *singleEventAverageCorrelationsVsQ2[nTestsMax][n2pCumulantTerms][3] = {{{NULL}}}; // [x,y,z components] TBI rename
 for(Int_t t=0;t<nTestsMax;t++)
 {
  if(!fFillBackgroundTEST[t]){continue;}
  for(Int_t ct=0;ct<n2pCumulantTerms;ct++)
  {
   for(Int_t xyz=0;xyz<3;xyz++)
   {
    singleEventAverageCorrelationsVsQ2[t][ct][xyz] = new TProfile(Form("singleEventAverageCorrelationsVsQ2[%d][%d][%d]",t,ct,xyz),"single-event averages",fnQ2bins,fnQ2min,fnQ2max);
    singleEventAverageCorrelationsVsQ2[t][ct][xyz]->GetXaxis()->SetTitle("Q_{2}");
    singleEventAverageCorrelationsVsQ2[t][ct][xyz]->GetYaxis()->SetTitle(Form("%s_{%s}",s2pCumulantTerms[ct].Data(),sXYZ[xyz].Data()));
   } // for(Int_t xyz=0;xyz<3;xyz++)
  } // for(Int_t ct=0;ct<n2pCumulantTerms;ct++)
 } // for(Int_t t=0;t<nTestsMax;t++)

 // ) Counters for N_\pi+ and N_\pi^-, for each Q2 bin (needed only for "Test 7" at the moment):
 TString *test7Counters[10000][2] = {{NULL}}; // N[44][0] is then N_\pi+ in the 44th Q2 bin in this event, while N[44][1] is then N_\pi- in the 44th Q2 bin in this event TBI hardwired 10000
 for(Int_t q2=0;q2<fnQ2bins;q2++)
 {
  test7Counters[q2][0] = new TString();
  *test7Counters[q2][0] = "";
  
  test7Counters[q2][1] = new TString();
  *test7Counters[q2][1] = "";
 }

 // c) Calculate averages <X1> (1st event), <X2> (2nd event) and <X1X2> (mixed-event), vs. Q2.
 Int_t nTracks1 = ca1->GetEntries();
 Int_t nTracks2 = ca2->GetEntries();
 for(Int_t iTrack1=0;iTrack1<nTracks1;iTrack1++)
 {
  AliAODTrack *atrack1 = dynamic_cast<AliAODTrack*>(ca1->UncheckedAt(iTrack1)); // TBI cross-check UncheckedAt
  // TBI Temporary track insanity checks:
  if(!atrack1){Fatal(sMethodName.Data(),"!atrack1");} // TBI keep this for some time, eventually just continue
  if(atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(atrack1->TestFilterBit(128) && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->TestFiletrBit(128) && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(!PassesCommonTrackCuts(atrack1)){continue;} // TBI re-think
  // Corresponding AOD global track:
  Int_t id1 = atrack1->GetID();
  AliAODTrack *gtrack1 = dynamic_cast<AliAODTrack*>(id1>=0 ? ca1->UncheckedAt(em1->GetValue(id1)) : ca1->UncheckedAt(em1->GetValue(-(id1+1))));
  if(!gtrack1){Fatal(sMethodName.Data(),"!gtrack1");} // TBI keep this for some time, eventually just continue
  // Common track selection criteria for all "normal" global tracks:
  if(!PassesGlobalTrackCuts(gtrack1)){continue;}

  // Start loop over tracks in the 2nd event:
  for(Int_t iTrack2=0;iTrack2<nTracks2;iTrack2++)
  {
   AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(ca2->UncheckedAt(iTrack2));
   // TBI Temporary track insanity checks:
   if(!atrack2){Fatal(sMethodName.Data(),"!atrack2");} // TBI keep this for some time, eventually just continue
   if(atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(atrack2->TestFilterBit(128) && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->TestFiletrBit(128) && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(!PassesCommonTrackCuts(atrack2)){continue;} // TBI re-think
   // Corresponding AOD global track:
   Int_t id2 = atrack2->GetID();
   AliAODTrack *gtrack2 = dynamic_cast<AliAODTrack*>(id2>=0 ? ca2->UncheckedAt(em2->GetValue(id2)) : ca2->UncheckedAt(em2->GetValue(-(id2+1))));
   if(!gtrack2){Fatal(sMethodName.Data(),"!gtrack2");} // TBI keep this for some time, eventually just continue
   // Common track selection criteria for all "normal" global tracks:
   if(!PassesGlobalTrackCuts(gtrack2)){continue;}

   // Okay... So we have a track from the second event ready. Let's decide whether kinemetics from 'atrack' or 'gtrack' is taken, let's check PID, and calculate the background:
   AliAODTrack *agtrack1 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
   AliAODTrack *agtrack2 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
   if(fFillControlHistogramsWithGlobalTrackInfo)
   {
    agtrack1 = gtrack1;
    agtrack2 = gtrack2;
   }
   else
   {
    agtrack1 = atrack1;
    agtrack2 = atrack2;
   } // if(fFillControlHistogramsWithGlobalTrackInfo)
   if(!agtrack1){Fatal(sMethodName.Data(),"!agtrack1");}
   if(!agtrack2){Fatal(sMethodName.Data(),"!agtrack2");}

   // Test 0: "Same charge pions, 2p correlations and cumulants projected onto k, for x, y and z components separately, not Lorentz invariant"
   if(fFillBackgroundTEST[0])
   {
    if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE)))
    {

     Double_t dk = RelativeMomenta(agtrack1,agtrack2); // relative momenta k = \frac{1}{2}|\vec{p_1}-\vec{p_2}|

     fBackgroundYieldTEST[0]->Fill(dk);

     // x-components:
     singleEventAverageCorrelationsVsQ2[0][0][0]->Fill(dk,agtrack1->Px()); // <X1>_x
     singleEventAverageCorrelationsVsQ2[0][1][0]->Fill(dk,agtrack2->Px()); // <X2>_x
     singleEventAverageCorrelationsVsQ2[0][2][0]->Fill(dk,agtrack1->Px()*agtrack2->Px()); // <X1X2>_x

     // y-components:
     singleEventAverageCorrelationsVsQ2[0][0][1]->Fill(dk,agtrack1->Py()); // <X1>_y
     singleEventAverageCorrelationsVsQ2[0][1][1]->Fill(dk,agtrack2->Py()); // <X2>_y
     singleEventAverageCorrelationsVsQ2[0][2][1]->Fill(dk,agtrack1->Py()*agtrack2->Py()); // <X1X2>_y

     // z-components:
     singleEventAverageCorrelationsVsQ2[0][0][2]->Fill(dk,agtrack1->Pz()); // <X1>_z
     singleEventAverageCorrelationsVsQ2[0][1][2]->Fill(dk,agtrack2->Pz()); // <X2>_z
     singleEventAverageCorrelationsVsQ2[0][2][2]->Fill(dk,agtrack1->Pz()*agtrack2->Pz()); // <X1X2>_z

    } // if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE)))
   } // if(fFillBackgroundTEST[0])

   // Test 1: "Same charge pions, 2p correlations and cumulants projected onto Lorentz invariant Q2"
   if(fFillBackgroundTEST[1])
   {
    if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE)))
    {

     Double_t dQ2 = Q2(agtrack1,agtrack2); // Lorentz invariant Q2

     // x-components:
     singleEventAverageCorrelationsVsQ2[1][0][0]->Fill(dQ2,agtrack1->Px()); // <X1>_x
     singleEventAverageCorrelationsVsQ2[1][1][0]->Fill(dQ2,agtrack2->Px()); // <X2>_x
     singleEventAverageCorrelationsVsQ2[1][2][0]->Fill(dQ2,agtrack1->Px()*agtrack2->Px()); // <X1X2>_x

     // y-components:
     singleEventAverageCorrelationsVsQ2[1][0][1]->Fill(dQ2,agtrack1->Py()); // <X1>_y
     singleEventAverageCorrelationsVsQ2[1][1][1]->Fill(dQ2,agtrack2->Py()); // <X2>_y
     singleEventAverageCorrelationsVsQ2[1][2][1]->Fill(dQ2,agtrack1->Py()*agtrack2->Py()); // <X1X2>_y

     // z-components:
     singleEventAverageCorrelationsVsQ2[1][0][2]->Fill(dQ2,agtrack1->Pz()); // <X1>_z
     singleEventAverageCorrelationsVsQ2[1][1][2]->Fill(dQ2,agtrack2->Pz()); // <X2>_z
     singleEventAverageCorrelationsVsQ2[1][2][2]->Fill(dQ2,agtrack1->Pz()*agtrack2->Pz()); // <X1X2>_z

    } // if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE)))
   } // if(fFillBackgroundTEST[1])

   // Test 4: "Same charge pions, 2p correlations and cumulants projected onto Lorentz invariant Q2, where energy prefactors have been taken into account, to make correlations and cumulants manifestly Lorenz invariant, even before projection"
   if(fFillBackgroundTEST[4])
   {
    if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE)))
    {

     Double_t E1 = agtrack1->E(); // energy of 1st track TBI check the mass hypothesis
     Double_t E2 = agtrack2->E(); // energy of 2nd track TBI check the mass hypothesis
     Double_t dQ2 = Q2(agtrack1,agtrack2); // Lorentz invariant Q2

     singleEventAverageCorrelationsVsQ2[4][0][0]->Fill(dQ2,E1); // <X1>
     singleEventAverageCorrelationsVsQ2[4][1][0]->Fill(dQ2,E2); // <X2>
     singleEventAverageCorrelationsVsQ2[4][2][0]->Fill(dQ2,E1*E2); // <X1X2>

    } // if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE)))

   } // if(fFillBackgroundTEST[4])

   // Test 6: "Same charge pions, E_{a}^{(b)} projected onto Lorentz invariant Q2 (see the write-up). This distribution shall be manifestly Lorentz invariant!"
   if(fFillBackgroundTEST[6])
   {
    if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE)))
    {
     // p_1:
     Double_t p1x = agtrack1->Px();
     Double_t p1y = agtrack1->Py();
     Double_t p1z = agtrack1->Pz();
     Double_t e1  = agtrack1->E();
     // p_2:
     Double_t p2x = agtrack2->Px();
     Double_t p2y = agtrack2->Py();
     Double_t p2z = agtrack2->Pz();
     Double_t e2  = agtrack2->E();  

     // Corresponding energy-momentum four-vectors:
     TLorentzVector track1(p1x,p1y,p1z,e1);
     TLorentzVector track2(p2x,p2y,p2z,e2);

     // Calculating L.I. E_{a}^{(b)} (i.e. energy of the 1st particle in the (rest) frame of the 2nd):
     TVector3 vb = TVector3(p2x/e2,p2y/e2,p2z/e2); // velocity of 2nd particle
     track1.Boost(-vb); // boosting the 1st track in the rest frame of 2nd
     //track2.Boost(-vb); // boosting the 2nd track in its rest frame 
     Double_t Eab = track1.E(); // energy of 1st track in the rest frame of 2nd. Yes, this is L.I. quantity

     // Fill Eab vs. dQ2 distribution, which is then manifestly L.I.
     Double_t dQ2 = Q2(agtrack1,agtrack2); // Lorentz invariant Q2
     singleEventAverageCorrelationsVsQ2[6][0][0]->Fill(dQ2,Eab); // <E_{a}^{(b)}> vs. Q2

     fEab_TEST6[1]->Fill(dQ2,Eab);

    } // if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE)))
   } // if(fFillBackgroundTEST[6])

   // Test 7: "Differential yield, see the write-up. Instead of counting # of pairs in Q2, I am counting # of particle "a" and # of particle "b" separately. 
   //          Then, the two independent observables are N_a and N_b, and they fluctuate e-b-e. That's how I get the notion of cumulants, and the projection 
   //          onto Q2 shall render it manifestly Lorentz invariant!"
   // To do:
   // a) Condidering in this test \pi+ \pi- only, later this will be generalized
   // b) Check explicitly the Lorentz invariance 
   // c) This clearly will work only for non-identical particles :'( Which is perhaps not than bad, given the fact that I can start to do my thing already with Lp+p-      
   // d) I cannot do anything for a given event, as it seems... :'( 
   if(fFillBackgroundTEST[7])
   {
    //if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE)) // TBI 20170724
    if( (Pion(gtrack1,1,kTRUE) || Pion(gtrack1,-1,kTRUE)) && (Proton(gtrack2,1,kTRUE) || Proton(gtrack2,-1,kTRUE)) )
    {
     // p_1:
     Double_t p1x = agtrack1->Px();
     Double_t p1y = agtrack1->Py();
     Double_t p1z = agtrack1->Pz();
     Double_t e1  = agtrack1->E();
     // p_2:
     Double_t p2x = agtrack2->Px();
     Double_t p2y = agtrack2->Py();
     Double_t p2z = agtrack2->Pz();
     Double_t e2  = agtrack2->E();  

     // Corresponding energy-momentum four-vectors:
     TLorentzVector track1(p1x,p1y,p1z,e1);
     TLorentzVector track2(p2x,p2y,p2z,e2);

     // Pseudo-code: 
     // 0) for each Q2 bin, formed from p_1 and p_2, count number of unique \pi+ and \pi-
     // 1) Then, treat N_\pi+ and N_\pi- and two independent observables, and after I run over all pairs, fill TProfile to get <N_\pi+ N_\pi-> - <N_\pi+> <N_\pi-> for each Q2 bin
     // 2) See if there is something non-trivial happening only at low Q2
     // 3) Do the same thing for background
     // 4) Generalize for 3-particles

     Double_t dQ2 = Q2(agtrack1,agtrack2); // Lorentz invariant Q2

     Int_t nBinQ2Number = this->BinNoForSpecifiedValue(fBackgroundTEST[7][0][0][0],dQ2);
     //cout<<"BIN_NO: "<<nBinQ2Number<<endl;
     //cout<<"ENERGIES: "<<agtrack1->E()<<" "<<agtrack2->E()<<endl;
     if(!test7Counters[nBinQ2Number][0]->Contains(Form("_%.6f_",agtrack1->E())))
     {
      *test7Counters[nBinQ2Number][0]+=Form("_%.6f_",agtrack1->E());
     }
     //cout<<"STRING1: "<<test7Counters[nBinQ2Number][0]->Data()<<endl;

     if(!test7Counters[nBinQ2Number][1]->Contains(Form("_%.6f_",agtrack2->E())))
     {
      *test7Counters[nBinQ2Number][1]+=Form("_%.6f_",agtrack2->E());
     }
     //cout<<"STRING2: "<<test7Counters[nBinQ2Number][1]->Data()<<endl;
  
    } // if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE))
   } // if(fFillBackgroundTEST[7])


  } // for(Int_t iTrack2=0;iTrack2<nTracks2;iTrack2++)
 } // for(Int_t iTrack1=0;iTrack1<nTracks1;iTrack1++)

 // d) Build all-event correlations and cumulants from single-event averages:
 for(Int_t t=0;t<nTestsMax;t++)
 {
  if(!fFillBackgroundTEST[t]){continue;}
  Int_t nBins = singleEventAverageCorrelationsVsQ2[t][0][0]->GetXaxis()->GetNbins();
  for(Int_t b=0;b<nBins;b++) // TBI at the moment, I use same binning for Q2 and Q3
  {
   for(Int_t xyz=0;xyz<3;xyz++)
   {
    Double_t dX1 = 0., dX1Err = 0.; // <X1>
    Double_t dX2 = 0., dX2Err = 0.; // <X2>
    Double_t dX1X2 = 0., dX1X2Err = 0.; // <X1X2>
    Double_t dC2 = 0., dC2Err = 0.; // C2 = <X1X2> - <X1><X2>
    Double_t dBinCenter = 0.; // TBI assuming binning is everywhere the same
    // 2p correlation terms:
    dX1 = singleEventAverageCorrelationsVsQ2[t][0][xyz]->GetBinContent(b+1); // <X1> vs. Q2
    dX1Err = singleEventAverageCorrelationsVsQ2[t][0][xyz]->GetBinError(b+1);
    dX2 = singleEventAverageCorrelationsVsQ2[t][1][xyz]->GetBinContent(b+1); // <X2> vs. Q2
    dX2Err = singleEventAverageCorrelationsVsQ2[t][1][xyz]->GetBinError(b+1);
    dX1X2 = singleEventAverageCorrelationsVsQ2[t][2][xyz]->GetBinContent(b+1); // <X1X2> vs. Q2
    dX1X2Err = singleEventAverageCorrelationsVsQ2[t][2][xyz]->GetBinError(b+1);
    if(TMath::Abs(dX1)>1.e-14 && TMath::Abs(dX2)>1.e-14 && TMath::Abs(dX1X2)>1.e-14) // basically, do calculations only if all terms in the 2p cumulant definition are non-zero
    {
     // 2p cumulant:
     dC2 = dX1X2 - dX1*dX2;
     dC2Err = 0.; // TBI propagate an error one day
     // fill 2p:
     dBinCenter = fBackgroundTEST[t][0][0][xyz]->GetBinCenter(b+1); // TBI assuming binning is everywhere the same. TBI it is used also below when filling 3p
     fBackgroundTEST[t][0][0][xyz]->Fill(dBinCenter,dX1); // <X1>
     fBackgroundTEST[t][0][1][xyz]->Fill(dBinCenter,dX2); // <X2>
     fBackgroundTEST[t][0][2][xyz]->Fill(dBinCenter,dX1X2); // <X1X2>
     fBackgroundCumulantsTEST[t][0][0][xyz]->Fill(dBinCenter,dC2); // C2 = <X1X2> - <X1><X2>
     //cout<<Form("bc: %.2f, xyz: %d, dX1: %.16f, dX2: %.16f, dX1dX2: %.16f, dC2: %.16f",dBinCenter,xyz,dX1,dX2,dX1X2,dC2)<<endl;
    }

    // TBI this was needed for 'Test 6', commenting out now temporarily:
    //dBinCenter = fBackgroundTEST[t][0][0][xyz]->GetBinCenter(b+1); // TBI assuming binning is everywhere the same. TBI it is used also below when filling 3p
    //fBackgroundTEST[t][0][0][xyz]->Fill(dBinCenter,dX1); // <X1>

   } // for(Int_t xyz=0;xyz<3;xyz++)
  } // for(Int_t b=0;b<nBins;b++) // TBI at the moment, I use same binning for Q2 and Q3

  if(7==t) // for "test 7" special treatment
  {
 
   for(Int_t b=1;b<=fnQ2bins;b++) // TBI check the boundaries
   {
    Int_t counter[2] = {0,0}; // [0=pions+][1=pions-]

    if(test7Counters[b][0] && !test7Counters[b][0]->EqualTo(""))
    {
     // Count pions+
     TObjArray *oa0 = TString(test7Counters[b][0]->Data()).Tokenize("_");
     while(oa0->At(counter[0]))
     {
      counter[0]++;
     }
     //cout<<Form("Bin %d; N = %d; STRING1: ",b,counter[0])<<test7Counters[b][0]->Data()<<endl;
    } // if(test7Counters[b][0] && !test7Counters[b][0]->EqualTo(""))
  
    if(test7Counters[b][1] && !test7Counters[b][1]->EqualTo(""))
    {
     // Count pions-
     TObjArray *oa1 = TString(test7Counters[b][1]->Data()).Tokenize("_");
     while(oa1->At(counter[1]))
     {
      counter[1]++;
     }
     //cout<<Form("Bin %d; N = %d; STRING2: ",b,counter[1])<<test7Counters[b][1]->Data()<<endl;
    } // if(test7Counters[b][1] && !test7Counters[b][1]->EqualTo("")) 

    Int_t nPionsP = counter[0]; // number of positive pions in bth Q2 bin, in this event, X1 in the cumulant formalism
    Int_t nPionsN = counter[1]; // number of negative pions in bth Q2 bin, in this event, X2 in the cumulant formalism
    if(nPionsP>=1 && nPionsN>=1)
    {
     Double_t dBinCenter = fBackgroundTEST[7][0][0][0]->GetBinCenter(b); // TBI assuming same binning everyrhere
     //TProfile *fBackgroundTEST[10][2][7][10]; //! [testNo][0=vs Q2, 1=vs Q3][example [0=<x1>][1=<x2>], ...,[6=<x1x2x3>]][differential index, e.g. for test 0 [0=Cx][1=Cy][2=Cz]]
     fBackgroundTEST[7][0][0][0]->Fill(dBinCenter,nPionsP); // X1 
     fBackgroundTEST[7][0][1][0]->Fill(dBinCenter,nPionsN); // X2
     fBackgroundTEST[7][0][2][0]->Fill(dBinCenter,nPionsP*nPionsN); // X1*X2 
    } // if(nPionsP>=1 && nPionsN>=1)
 
   } // for(Int_t b=1;b<=fnQ2bins;b++) // TBI check the boundaries

  } // if(7==t) // for "test 7" special treatment


 } // for(Int_t t=0;t<nTestsMax;t++)


 // e) Delete local TProfile's to hold single-event averages:
 for(Int_t t=0;t<nTestsMax;t++)
 {
  if(!fFillBackgroundTEST[t]){continue;}
  for(Int_t ct=0;ct<n2pCumulantTerms;ct++)
  {
   for(Int_t xyz=0;xyz<3;xyz++)
   {
    delete singleEventAverageCorrelationsVsQ2[t][ct][xyz];
   }
  }
 }

} // void AliAnalysisTaskMultiparticleFemtoscopy::Calculate2pBackgroundTEST(TClonesArray *ca1, TClonesArray *ca2, TExMap *em1, TExMap *em2)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::Calculate3pBackgroundTEST(TClonesArray *ca1, TClonesArray *ca2, TClonesArray *ca3, TExMap *em1, TExMap *em2, TExMap *em3)
{
 // Calculate 3p TEST background.

 // a) Insanity checks;
 // b) Book local TProfile's to hold 1x1x1 event averages vs. Q3;
 // c) Calculate averages all mixed-event averages vs. Q3;
 // d) Build all-event correlations and cumulants from single-event averages.

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::Calculate3pBackgroundTEST(TClonesArray *ca1, TClonesArray *ca2, TClonesArray *ca3, TExMap *em1, TExMap *em2, TExMap *em3)";
 if(!ca1){Fatal(sMethodName.Data(),"!ca1");}
 if(!ca2){Fatal(sMethodName.Data(),"!ca2");}
 if(!ca3){Fatal(sMethodName.Data(),"!ca3");}
 if(!em1){Fatal(sMethodName.Data(),"!em1");}
 if(!em2){Fatal(sMethodName.Data(),"!em2");}
 if(!em3){Fatal(sMethodName.Data(),"!em3");}

 // b) Book local TProfile's to hold 1x1x1 event averages vs. Q3:
 const Int_t nTestsMax = 10; // see what is hardwired in .h file
 TString sXYZ[3] = {"x","y","z"};
 const Int_t n3pCumulantTerms = 7;
 TString s3pCumulantTerms[n3pCumulantTerms] = {"#LTX_{1}#GT","#LTX_{2}#GT","#LTX_{3}#GT","#LTX_{1}X_{2}#GT","#LTX_{1}X_{3}#GT","#LTX_{2}X_{3}#GT","#LTX_{1}X_{2}X_{3}#GT"};
 TProfile *singleEventAverageCorrelationsVsQ3[nTestsMax][n3pCumulantTerms][3] = {{{NULL}}}; // [x,y,z components]
 for(Int_t t=0;t<nTestsMax;t++)
 {
  if(!fFillBackgroundTEST[t]){continue;}
  for(Int_t ct=0;ct<n3pCumulantTerms;ct++)
  {
   for(Int_t xyz=0;xyz<3;xyz++)
   {
    singleEventAverageCorrelationsVsQ3[t][ct][xyz] = new TProfile(Form("singleEventAverageCorrelationsVsQ3[%d][%d][%d]",t,ct,xyz),"single-event averages",fnQ3bins,fnQ3min,fnQ3max);
    singleEventAverageCorrelationsVsQ3[t][ct][xyz]->GetXaxis()->SetTitle("Q_{3}");
    singleEventAverageCorrelationsVsQ3[t][ct][xyz]->GetYaxis()->SetTitle(Form("%s_{%s}",s3pCumulantTerms[ct].Data(),sXYZ[xyz].Data()));
   } // for(Int_t xyz=0;xyz<3;xyz++)
  } // for(Int_t ct=0;ct<n3pCumulantTerms;ct++)
 }

 // ) Differential yield counters (needed only for "Test 2" at the moment):
 const Int_t nLoopsVsQ2 = 2;
 const Int_t nLoopsVsQ3 = 3;
 const Int_t nMaxEntries = 10000;
 TString *dycVsQ2[nLoopsVsQ2][nMaxEntries] = {{NULL}};
 TString *dycVsQ3[nLoopsVsQ3][nMaxEntries] = {{NULL}};
 for(Int_t nl=0;nl<nLoopsVsQ2;nl++)
 {
  for(Int_t nme=1;nme<=nMaxEntries;nme++) // yes, from 1, since I do not have later bin 0
  {
   if(nme > fnQ2bins){break;}
   dycVsQ2[nl][nme] = new TString();
  }
 }
 for(Int_t nl=0;nl<nLoopsVsQ3;nl++)
 {
  for(Int_t nme=1;nme<=nMaxEntries;nme++) // yes, from 1, since I do not have later bin 0
  {
   if(nme > fnQ3bins){break;}
   dycVsQ3[nl][nme] = new TString();
  }
 }

 // c) Calculate averages <X1> (1st event), <X2> (2nd event) and <X1X2> (mixed-event), vs. Q2.
 Int_t nTracks1 = ca1->GetEntries();
 Int_t nTracks2 = ca2->GetEntries();
 Int_t nTracks3 = ca3->GetEntries();
 for(Int_t iTrack1=0;iTrack1<nTracks1;iTrack1++)
 {
  AliAODTrack *atrack1 = dynamic_cast<AliAODTrack*>(ca1->UncheckedAt(iTrack1)); // TBI cross-check UncheckedAt
  // TBI Temporary track insanity checks:
  if(!atrack1){Fatal(sMethodName.Data(),"!atrack1");} // TBI keep this for some time, eventually just continue
  if(atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(atrack1->TestFilterBit(128) && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->TestFiletrBit(128) && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(!PassesCommonTrackCuts(atrack1)){continue;} // TBI re-think
  // Corresponding AOD global track:
  Int_t id1 = atrack1->GetID();
  AliAODTrack *gtrack1 = dynamic_cast<AliAODTrack*>(id1>=0 ? ca1->UncheckedAt(em1->GetValue(id1)) : ca1->UncheckedAt(em1->GetValue(-(id1+1))));
  if(!gtrack1){Fatal(sMethodName.Data(),"!gtrack1");} // TBI keep this for some time, eventually just continue
  // Common track selection criteria for all "normal" global tracks:
  if(!PassesGlobalTrackCuts(gtrack1)){continue;}

  // Start loop over tracks in the 2nd event:
  for(Int_t iTrack2=0;iTrack2<nTracks2;iTrack2++)
  {
   AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(ca2->UncheckedAt(iTrack2));
   // TBI Temporary track insanity checks:
   if(!atrack2){Fatal(sMethodName.Data(),"!atrack2");} // TBI keep this for some time, eventually just continue
   if(atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(atrack2->TestFilterBit(128) && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->TestFiletrBit(128) && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(!PassesCommonTrackCuts(atrack2)){continue;} // TBI re-think
   // Corresponding AOD global track:
   Int_t id2 = atrack2->GetID();
   AliAODTrack *gtrack2 = dynamic_cast<AliAODTrack*>(id2>=0 ? ca2->UncheckedAt(em2->GetValue(id2)) : ca2->UncheckedAt(em2->GetValue(-(id2+1))));
   if(!gtrack2){Fatal(sMethodName.Data(),"!gtrack2");} // TBI keep this for some time, eventually just continue
   // Common track selection criteria for all "normal" global tracks:
   if(!PassesGlobalTrackCuts(gtrack2)){continue;}

   // Start loop over tracks in the 3rd event:
   for(Int_t iTrack3=0;iTrack3<nTracks3;iTrack3++)
   {
    AliAODTrack *atrack3 = dynamic_cast<AliAODTrack*>(ca3->UncheckedAt(iTrack3));
    // TBI Temporary track insanity checks:
    if(!atrack3){Fatal(sMethodName.Data(),"!atrack3");} // TBI keep this for some time, eventually just continue
    if(atrack3->GetID()>=0 && atrack3->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack3->GetID()>=0 && atrack3->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
    if(atrack3->TestFilterBit(128) && atrack3->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack3->TestFiletrBit(128) && atrack3->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
    if(!PassesCommonTrackCuts(atrack3)){continue;} // TBI re-think
    // Corresponding AOD global track:
    Int_t id3 = atrack3->GetID();
    AliAODTrack *gtrack3 = dynamic_cast<AliAODTrack*>(id3>=0 ? ca3->UncheckedAt(em3->GetValue(id3)) : ca3->UncheckedAt(em3->GetValue(-(id3+1))));
    if(!gtrack3){Fatal(sMethodName.Data(),"!gtrack3");} // TBI keep this for some time, eventually just continue
    // Common track selection criteria for all "normal" global tracks:
    if(!PassesGlobalTrackCuts(gtrack3)){continue;}

    // Okay... So we have three tracks from three different events ready. Let's decide whether kinemetics from 'atrack' or 'gtrack' is taken, let's check PID, and calculate the background:
    AliAODTrack *agtrack1 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
    AliAODTrack *agtrack2 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
    AliAODTrack *agtrack3 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
    if(fFillControlHistogramsWithGlobalTrackInfo)
    {
     agtrack1 = gtrack1;
     agtrack2 = gtrack2;
     agtrack3 = gtrack3;
    }
    else
    {
     agtrack1 = atrack1;
     agtrack2 = atrack2;
     agtrack3 = atrack3;
    } // if(fFillControlHistogramsWithGlobalTrackInfo)
    if(!agtrack1){Fatal(sMethodName.Data(),"!agtrack1");}
    if(!agtrack2){Fatal(sMethodName.Data(),"!agtrack2");}
    if(!agtrack3){Fatal(sMethodName.Data(),"!agtrack3");}

    // Test 0: "Same charge pions, 2p correlations and cumulants projected onto Q2, for x, y and z components separately"
    if(fFillBackgroundTEST[0])
    {
     if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE)))
     {
      // Q3:
      Double_t dQ3 = Q3(agtrack1,agtrack2,agtrack3);

      fBackgroundYieldTEST[1]->Fill(dQ3);

      // x-components:
      singleEventAverageCorrelationsVsQ3[0][0][0]->Fill(dQ3,agtrack1->Px()); // <X1>_x
      singleEventAverageCorrelationsVsQ3[0][1][0]->Fill(dQ3,agtrack2->Px()); // <X2>_x
      singleEventAverageCorrelationsVsQ3[0][2][0]->Fill(dQ3,agtrack3->Px()); // <X3>_x
      singleEventAverageCorrelationsVsQ3[0][3][0]->Fill(dQ3,agtrack1->Px()*agtrack2->Px()); // <X1X2>_x
      singleEventAverageCorrelationsVsQ3[0][4][0]->Fill(dQ3,agtrack1->Px()*agtrack3->Px()); // <X1X3>_x
      singleEventAverageCorrelationsVsQ3[0][5][0]->Fill(dQ3,agtrack2->Px()*agtrack3->Px()); // <X2X3>_x
      singleEventAverageCorrelationsVsQ3[0][6][0]->Fill(dQ3,agtrack1->Px()*agtrack2->Px()*agtrack3->Px()); // <X1X2X3>_x

      // y-components:
      singleEventAverageCorrelationsVsQ3[0][0][1]->Fill(dQ3,agtrack1->Py()); // <X1>_y
      singleEventAverageCorrelationsVsQ3[0][1][1]->Fill(dQ3,agtrack2->Py()); // <X2>_y
      singleEventAverageCorrelationsVsQ3[0][2][1]->Fill(dQ3,agtrack3->Py()); // <X3>_y
      singleEventAverageCorrelationsVsQ3[0][3][1]->Fill(dQ3,agtrack1->Py()*agtrack2->Py()); // <X1X2>_y
      singleEventAverageCorrelationsVsQ3[0][4][1]->Fill(dQ3,agtrack1->Py()*agtrack3->Py()); // <X1X3>_y
      singleEventAverageCorrelationsVsQ3[0][5][1]->Fill(dQ3,agtrack2->Py()*agtrack3->Py()); // <X2X3>_y
      singleEventAverageCorrelationsVsQ3[0][6][1]->Fill(dQ3,agtrack1->Py()*agtrack2->Py()*agtrack3->Py()); // <X1X2X3>_y

      // z-components:
      singleEventAverageCorrelationsVsQ3[0][0][2]->Fill(dQ3,agtrack1->Pz()); // <X1>_z
      singleEventAverageCorrelationsVsQ3[0][1][2]->Fill(dQ3,agtrack2->Pz()); // <X2>_z
      singleEventAverageCorrelationsVsQ3[0][2][2]->Fill(dQ3,agtrack3->Pz()); // <X3>_z
      singleEventAverageCorrelationsVsQ3[0][3][2]->Fill(dQ3,agtrack1->Pz()*agtrack2->Pz()); // <X1X2>_z
      singleEventAverageCorrelationsVsQ3[0][4][2]->Fill(dQ3,agtrack1->Pz()*agtrack3->Pz()); // <X1X3>_z
      singleEventAverageCorrelationsVsQ3[0][5][2]->Fill(dQ3,agtrack2->Pz()*agtrack3->Pz()); // <X2X3>_z
      singleEventAverageCorrelationsVsQ3[0][6][2]->Fill(dQ3,agtrack1->Pz()*agtrack2->Pz()*agtrack3->Pz()); // <X1X2X3>_z

     } // if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE)))
    } // if(fFillBackgroundTEST[0])

    // Test 1: "Same charge pions, 2p correlations and cumulants projected onto Lorentz invariant Q3"
    if(fFillBackgroundTEST[1])
    {
     if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE)))
     {
      // Q3:
      Double_t dQ3 = Q3(agtrack1,agtrack2,agtrack3); // This is Lorentz invariant

      // x-components:
      singleEventAverageCorrelationsVsQ3[1][0][0]->Fill(dQ3,agtrack1->Px()); // <X1>_x
      singleEventAverageCorrelationsVsQ3[1][1][0]->Fill(dQ3,agtrack2->Px()); // <X2>_x
      singleEventAverageCorrelationsVsQ3[1][2][0]->Fill(dQ3,agtrack3->Px()); // <X3>_x
      singleEventAverageCorrelationsVsQ3[1][3][0]->Fill(dQ3,agtrack1->Px()*agtrack2->Px()); // <X1X2>_x
      singleEventAverageCorrelationsVsQ3[1][4][0]->Fill(dQ3,agtrack1->Px()*agtrack3->Px()); // <X1X3>_x
      singleEventAverageCorrelationsVsQ3[1][5][0]->Fill(dQ3,agtrack2->Px()*agtrack3->Px()); // <X2X3>_x
      singleEventAverageCorrelationsVsQ3[1][6][0]->Fill(dQ3,agtrack1->Px()*agtrack2->Px()*agtrack3->Px()); // <X1X2X3>_x

      // y-components:
      singleEventAverageCorrelationsVsQ3[1][0][1]->Fill(dQ3,agtrack1->Py()); // <X1>_y
      singleEventAverageCorrelationsVsQ3[1][1][1]->Fill(dQ3,agtrack2->Py()); // <X2>_y
      singleEventAverageCorrelationsVsQ3[1][2][1]->Fill(dQ3,agtrack3->Py()); // <X3>_y
      singleEventAverageCorrelationsVsQ3[1][3][1]->Fill(dQ3,agtrack1->Py()*agtrack2->Py()); // <X1X2>_y
      singleEventAverageCorrelationsVsQ3[1][4][1]->Fill(dQ3,agtrack1->Py()*agtrack3->Py()); // <X1X3>_y
      singleEventAverageCorrelationsVsQ3[1][5][1]->Fill(dQ3,agtrack2->Py()*agtrack3->Py()); // <X2X3>_y
      singleEventAverageCorrelationsVsQ3[1][6][1]->Fill(dQ3,agtrack1->Py()*agtrack2->Py()*agtrack3->Py()); // <X1X2X3>_y

      // z-components:
      singleEventAverageCorrelationsVsQ3[1][0][2]->Fill(dQ3,agtrack1->Pz()); // <X1>_z
      singleEventAverageCorrelationsVsQ3[1][1][2]->Fill(dQ3,agtrack2->Pz()); // <X2>_z
      singleEventAverageCorrelationsVsQ3[1][2][2]->Fill(dQ3,agtrack3->Pz()); // <X3>_z
      singleEventAverageCorrelationsVsQ3[1][3][2]->Fill(dQ3,agtrack1->Pz()*agtrack2->Pz()); // <X1X2>_z
      singleEventAverageCorrelationsVsQ3[1][4][2]->Fill(dQ3,agtrack1->Pz()*agtrack3->Pz()); // <X1X3>_z
      singleEventAverageCorrelationsVsQ3[1][5][2]->Fill(dQ3,agtrack2->Pz()*agtrack3->Pz()); // <X2X3>_z
      singleEventAverageCorrelationsVsQ3[1][6][2]->Fill(dQ3,agtrack1->Pz()*agtrack2->Pz()*agtrack3->Pz()); // <X1X2X3>_z

     } // if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE)))
    } // if(fFillBackgroundTEST[1])

    // Test 2: "Same charge pions-kaons-protons, 2p correlations and cumulants projected onto Lorentz invariant Q3, using yield, and not momenta like in 'Test 0' and 'Test 1'"
    //         the first particle in the loop is always pion, the second is kaon, the third is proton. [xyz] is always se to 0
    if(fFillBackgroundTEST[2])
    {
     if((Pion(gtrack1,1,kTRUE) && Kaon(gtrack2,1,kTRUE) && Proton(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Kaon(gtrack2,-1,kTRUE) && Proton(gtrack3,-1,kTRUE)))
     //if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE)))
     {
      // Count the yield differentially:
      Double_t dQ3 = Q3(agtrack1,agtrack2,agtrack3); // This is Lorentz invariant
      Int_t binNoQ3 = this->BinNoForSpecifiedValue(singleEventAverageCorrelationsVsQ3[2][0][0],dQ3); // TBI cross-check with some other profile
      if(binNoQ3>0) // otherwise it's either underflow or overflow
      {
       TString pattern1 = Form(":track%d:",iTrack1);
       TString pattern2 = Form(":track%d:",iTrack2);
       TString pattern3 = Form(":track%d:",iTrack3);
       if(!dycVsQ3[0][binNoQ3]->Contains(pattern1.Data()))
       {
        *dycVsQ3[0][binNoQ3] = Form("%s%s",dycVsQ3[0][binNoQ3]->Data(),pattern1.Data());
       }
       if(!dycVsQ3[1][binNoQ3]->Contains(pattern2.Data()))
       {
        *dycVsQ3[1][binNoQ3] = Form("%s%s",dycVsQ3[1][binNoQ3]->Data(),pattern2.Data());
       }
       if(!dycVsQ3[2][binNoQ3]->Contains(pattern3.Data()))
       {
        *dycVsQ3[2][binNoQ3] = Form("%s%s",dycVsQ3[2][binNoQ3]->Data(),pattern3.Data());
       }
      } // if(binNoQ3>0)

     } // if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE)))
    } // if(fFillBackgroundTEST[2])


   // Test 4: "Same charge pions, 2p correlations and cumulants projected onto Lorentz invariant Q2, where energy prefactors have been taken into account, to make correlations and cumulants manifestly Lorenz invariant, even before projection"
    if(fFillBackgroundTEST[4])
    {
     if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE)))
     {
      // Q3:
      Double_t E1 = agtrack1->E(); // energy of 1st track TBI check the mass hypothesis
      Double_t E2 = agtrack2->E(); // energy of 2nd track TBI check the mass hypothesis
      Double_t E3 = agtrack3->E(); // energy of 2nd track TBI check the mass hypothesis
      Double_t dQ3 = Q3(agtrack1,agtrack2,agtrack3); // This is Lorentz invariant

      singleEventAverageCorrelationsVsQ3[4][0][0]->Fill(dQ3,E1); // <X1>
      singleEventAverageCorrelationsVsQ3[4][1][0]->Fill(dQ3,E2); // <X2>
      singleEventAverageCorrelationsVsQ3[4][2][0]->Fill(dQ3,E3); // <X3>
      singleEventAverageCorrelationsVsQ3[4][3][0]->Fill(dQ3,E1*E2); // <X1X2>
      singleEventAverageCorrelationsVsQ3[4][4][0]->Fill(dQ3,E1*E3); // <X1X3>
      singleEventAverageCorrelationsVsQ3[4][5][0]->Fill(dQ3,E2*E3); // <X2X3>
      singleEventAverageCorrelationsVsQ3[4][6][0]->Fill(dQ3,E1*E2*E3); // <X1X2X3>
     } // if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE)))
    } // if(fFillBackgroundTEST[4])

   } // for(Int_t iTrack3=0;iTrack3<nTracks3;iTrack3++)
  } // for(Int_t iTrack2=0;iTrack2<nTracks2;iTrack2++)
 } // for(Int_t iTrack1=0;iTrack1<nTracks1;iTrack1++)

 // d) Build all-event correlations and cumulants from single-event averages:
 for(Int_t t=0;t<nTestsMax;t++)
 {
  if(!fFillBackgroundTEST[t]){continue;}
  Int_t nBins = singleEventAverageCorrelationsVsQ3[t][0][0]->GetXaxis()->GetNbins();
  for(Int_t b=0;b<nBins;b++) // TBI at the moment, I use same binning for Q2 and Q3
  {
   for(Int_t xyz=0;xyz<3;xyz++)
   {
    Double_t dX1 = 0., dX1Err = 0.; // <X1>
    Double_t dX2 = 0., dX2Err = 0.; // <X2>
    Double_t dX3 = 0., dX3Err = 0.; // <X3>
    Double_t dX1X2 = 0., dX1X2Err = 0.; // <X1X2>
    Double_t dX1X3 = 0., dX1X3Err = 0.; // <X1X3>
    Double_t dX2X3 = 0., dX2X3Err = 0.; // <X2X3>
    Double_t dX1X2X3 = 0., dX1X2X3Err = 0.; // <X1X2X3>
    Double_t dC12 = 0., dC12Err = 0.; // C12 = <X1X2> - <X1><X2> vs. Q3
    Double_t dC13 = 0., dC13Err = 0.; // C13 = <X1X3> - <X1><X3> vs. Q3
    Double_t dC23 = 0., dC23Err = 0.; // C23 = <X2X3> - <X2><X3> vs. Q3
    Double_t dC3 = 0., dC3Err = 0.; // C3 = <X1X2X3> - <X1X2><X3> - <X1X3><X2> - <X2X3><X1> + 2<X1><X2><X3>
    Double_t dBinCenter = 0.; // TBI assuming binning is everywhere the same

    // 3p correlation terms (note that now binning is vs. Q3 TBI yes, fine for the time being, but this is a landmine clearly...):
    dX1 = singleEventAverageCorrelationsVsQ3[t][0][xyz]->GetBinContent(b+1); // <X1> vs. Q3
    dX1Err = singleEventAverageCorrelationsVsQ3[t][0][xyz]->GetBinError(b+1);
    dX2 = singleEventAverageCorrelationsVsQ3[t][1][xyz]->GetBinContent(b+1); // <X2> vs. Q3
    dX2Err = singleEventAverageCorrelationsVsQ3[t][1][xyz]->GetBinError(b+1);
    dX3 = singleEventAverageCorrelationsVsQ3[t][2][xyz]->GetBinContent(b+1); // <X3> vs. Q3
    dX3Err = singleEventAverageCorrelationsVsQ3[t][2][xyz]->GetBinError(b+1);
    dX1X2 = singleEventAverageCorrelationsVsQ3[t][3][xyz]->GetBinContent(b+1); // <X1X2> vs. Q3
    dX1X2Err = singleEventAverageCorrelationsVsQ3[t][3][xyz]->GetBinError(b+1);
    dX1X3 = singleEventAverageCorrelationsVsQ3[t][4][xyz]->GetBinContent(b+1); // <X1X3> vs. Q3
    dX1X3Err = singleEventAverageCorrelationsVsQ3[t][4][xyz]->GetBinError(b+1);
    dX2X3 = singleEventAverageCorrelationsVsQ3[t][5][xyz]->GetBinContent(b+1); // <X2X3> vs. Q3
    dX2X3Err = singleEventAverageCorrelationsVsQ3[t][5][xyz]->GetBinError(b+1);
    dX1X2X3 = singleEventAverageCorrelationsVsQ3[t][6][xyz]->GetBinContent(b+1); // <X1X2X3> vs. Q3
    dX1X2X3Err = singleEventAverageCorrelationsVsQ3[t][6][xyz]->GetBinError(b+1);
    if(TMath::Abs(dX1)>1.e-14 && TMath::Abs(dX2)>1.e-14 && TMath::Abs(dX3)>1.e-14 &&
       TMath::Abs(dX1X2)>1.e-14 && TMath::Abs(dX1X3)>1.e-14 && TMath::Abs(dX2X3)>1.e-14 &&
       TMath::Abs(dX1X2X3)>1.e-14) // basically, do calculations only if all terms in the 3p cumulant definition are non-zero
    {
     // 2p cumulants vs. Q3:
     dC12 = dX1X2 - dX1*dX2; // C12 = <X1X2> - <X1><X2> vs. Q3
     dC12Err = 0.; // TBI propagate an error one day
     dC13 = dX1X3 - dX1*dX3; // C13 = <X1X3> - <X1><X3> vs. Q3
     dC13Err = 0.; // TBI propagate an error one day
     dC23 = dX2X3 - dX2*dX3; // C23 = <X2X3> - <X2><X3> vs. Q3
     dC23Err = 0.; // TBI propagate an error one day
     // 3p cumulant vs. Q3:
     dC3 = dX1X2X3 - dX1X2*dX3 - dX1X3*dX2 - dX2X3*dX1 + 2.*dX1*dX2*dX3;
     dC3Err = 0.; // TBI propagate an error one day
     // fill 3p:

     dBinCenter = fBackgroundTEST[t][1][0][xyz]->GetBinCenter(b+1); // TBI assuming binning is everywhere the same. TBI it is used also below when filling 3p
     fBackgroundTEST[t][1][0][xyz]->Fill(dBinCenter,dX1); // <X1> vs. Q3
     fBackgroundTEST[t][1][1][xyz]->Fill(dBinCenter,dX2); // <X2> vs. Q3
     fBackgroundTEST[t][1][2][xyz]->Fill(dBinCenter,dX3); // <X3> vs. Q3
     fBackgroundTEST[t][1][3][xyz]->Fill(dBinCenter,dX1X2); // <X1X2> vs. Q3
     fBackgroundTEST[t][1][4][xyz]->Fill(dBinCenter,dX1X3); // <X1X3> vs. Q3
     fBackgroundTEST[t][1][5][xyz]->Fill(dBinCenter,dX2X3); // <X2X3> vs. Q3
     fBackgroundTEST[t][1][6][xyz]->Fill(dBinCenter,dX1X2X3); // <X1X2X3> vs. Q3
     fBackgroundCumulantsTEST[t][1][0][xyz]->Fill(dBinCenter,dC12); // <X1X2> - <X1><X2> vs. Q3
     fBackgroundCumulantsTEST[t][1][1][xyz]->Fill(dBinCenter,dC13); // <X1X3> - <X1><X3> vs. Q3
     fBackgroundCumulantsTEST[t][1][2][xyz]->Fill(dBinCenter,dC23); // <X2X3> - <X2><X3> vs. Q3
     fBackgroundCumulantsTEST[t][1][3][xyz]->Fill(dBinCenter,dC3); // C3 = <X1X2X3> - <X1X2><X3> - <X1X3><X2> - <X2X3><X1> + 2<X1><X2><X3> vs. Q3
    } // if(TMath::Abs(dX1)>1.e-14 && TMath::Abs(dX2)>1.e-14 && TMath::Abs(dX3)>1.e-14 && ...
   } // for(Int_t xyz=0;xyz<3;xyz++)
  } // for(Int_t b=0;b<nBins;b++) // TBI at the moment, I use same binning for Q2 and Q3


  if(2==t) // for "test 2" I just need correlations, and build cumulants at the end of the day. [x = pion yield][y = kaon yield][z = proton yield]
  {
   for(Int_t b=1;b<=nBins;b++) // TBI at the moment, I use same binning for Q2 and Q3. Note that here I start a loop over 1, not like in previous two tests
   {
    Double_t dBinCenter = fBackgroundTEST[2][1][0][0]->GetBinCenter(b); // TBI assuming binning is everywhere the same. TBI it is used also below when filling 3p

    Int_t counter[3] = {0}; // [0=pions][1=kaons][2=protons]
    for(Int_t pkp=0;pkp<3;pkp++)
    {
     TObjArray *oa = TString(dycVsQ3[pkp][b]->Data()).Tokenize(":");
     while(oa->At(counter[pkp]))
     {
      counter[pkp]++;
     }
    }

    Int_t dX1 = counter[0]; // X1 (not the average!)
    Int_t dX2 = counter[1]; // X2 (not the average!)
    Int_t dX3 = counter[2]; // X3 (not the average!)
    Int_t dX1X2 = counter[0]*counter[1]; // X1X2
    Int_t dX1X3 = counter[0]*counter[2]; // X1X3
    Int_t dX2X3 = counter[1]*counter[2]; // X2X3
    Int_t dX1X2X3 = counter[0]*counter[1]*counter[2]; // X1X2X3
    fBackgroundTEST[2][1][0][0]->Fill(dBinCenter,dX1); // <X1>
    fBackgroundTEST[2][1][1][0]->Fill(dBinCenter,dX2); // <X2>
    fBackgroundTEST[2][1][2][0]->Fill(dBinCenter,dX3); // <X3>
    fBackgroundTEST[2][1][3][0]->Fill(dBinCenter,dX1X2); // <X1X2>
    fBackgroundTEST[2][1][4][0]->Fill(dBinCenter,dX1X3); // <X1X3>
    fBackgroundTEST[2][1][5][0]->Fill(dBinCenter,dX2X3); // <X2X3>
    fBackgroundTEST[2][1][6][0]->Fill(dBinCenter,dX1X2X3); // <X1X2X3>

   } // for(Int_t b=1;b<=nBins;b++)

  } // if(2==t) // for "test 2" I just need correlations, and build cumulants at the end of the day


 } // for(Int_t t=0;t<nTestsMax;t++)

 // e) Delete local TProfile's to hold single-event averages:
 for(Int_t t=0;t<nTestsMax;t++)
 {
  if(!fFillBackgroundTEST[t]){continue;}
  for(Int_t ct=0;ct<n3pCumulantTerms;ct++)
  {
   for(Int_t xyz=0;xyz<3;xyz++)
   {
    delete singleEventAverageCorrelationsVsQ3[t][ct][xyz];
   }
  }
 } // for(Int_t t=0;t<nTestsMax;t++)

 for(Int_t nl=0;nl<nLoopsVsQ2;nl++)
 {
  for(Int_t nme=1;nme<=nMaxEntries;nme++) // yes, from 1, see above
  {
   if(nme > fnQ2bins){break;}
   delete dycVsQ2[nl][nme];
  }
 }
 for(Int_t nl=0;nl<nLoopsVsQ3;nl++)
 {
  for(Int_t nme=1;nme<=nMaxEntries;nme++) // yes, from 1, see above
  {
   if(nme > fnQ3bins){break;}
   delete dycVsQ3[nl][nme];
  }
 }

} // void AliAnalysisTaskMultiparticleFemtoscopy::Calculate3pBackgroundTEST(TClonesArray *ca1, TClonesArray *ca2, TClonesArray *ca3, TExMap *em1, TExMap *em2, TExMap *em3)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::Calculate3pBackground(TClonesArray *ca1, TClonesArray *ca2, TClonesArray *ca3)
{
 // Calculate background.

 // a) Insanity checks;
 // b) Temporary primitive algorithm: correlate particles from 'previous event' to particles from 'current event';
 // c) Three nested loops to calculate B(k);

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::Calculate3pBackground(TClonesArray *ca1, TClonesArray *ca2, TClonesArray *ca3)";
 if(!ca1){Fatal(sMethodName.Data(),"!ca1");}
 if(!ca2){Fatal(sMethodName.Data(),"!ca2");}
 if(!ca3){Fatal(sMethodName.Data(),"!ca3");}

 // b) Temporary primitive algorithm: correlate particles from 'previous event' to particles from 'current event'
 // ...

 // c) Three nested loops to calculate B(k):
 Int_t nTracks1 = ca1->GetEntries();
 Int_t nTracks2 = ca2->GetEntries();
 Int_t nTracks3 = ca3->GetEntries();

 // Start loop over tracks in the 1st event:
 for(Int_t iTrack1=0;iTrack1<nTracks1;iTrack1++)
 {
  AliAODTrack *atrack1 = dynamic_cast<AliAODTrack*>(ca1->UncheckedAt(iTrack1)); // TBI cross-check UncheckedAt
  // TBI Temporary track insanity checks:
  if(!atrack1){Fatal(sMethodName.Data(),"!atrack1");} // TBI keep this for some time, eventually just continue
  if(atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(atrack1->TestFilterBit(128) && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->TestFiletrBit(128) && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(!PassesCommonTrackCuts(atrack1)){continue;} // TBI re-think
  // Corresponding AOD global track:
  Int_t id1 = atrack1->GetID();
  AliAODTrack *gtrack1 = dynamic_cast<AliAODTrack*>(id1>=0 ? ca1->UncheckedAt(fGlobalTracksAOD[1]->GetValue(id1)) : ca1->UncheckedAt(fGlobalTracksAOD[1]->GetValue(-(id1+1))));
  if(!gtrack1){Fatal(sMethodName.Data(),"!gtrack1");} // TBI keep this for some time, eventually just continue
  // Common track selection criteria for all "normal" global tracks:
  if(!PassesGlobalTrackCuts(gtrack1)){continue;}

  // Start loop over tracks in the 2nd event:
  for(Int_t iTrack2=0;iTrack2<nTracks2;iTrack2++)
  {
   AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(ca2->UncheckedAt(iTrack2));
   // TBI Temporary track insanity checks:
   if(!atrack2){Fatal(sMethodName.Data(),"!atrack2");} // TBI keep this for some time, eventually just continue
   if(atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(atrack2->TestFilterBit(128) && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->TestFiletrBit(128) && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(!PassesCommonTrackCuts(atrack2)){continue;} // TBI re-think
   // Corresponding AOD global track:
   Int_t id2 = atrack2->GetID();
   AliAODTrack *gtrack2 = dynamic_cast<AliAODTrack*>(id2>=0 ? ca2->UncheckedAt(fGlobalTracksAOD[2]->GetValue(id2)) : ca2->UncheckedAt(fGlobalTracksAOD[2]->GetValue(-(id2+1))));
   if(!gtrack2){Fatal(sMethodName.Data(),"!gtrack2");} // TBI keep this for some time, eventually just continue
   // Common track selection criteria for all "normal" global tracks:
   if(!PassesGlobalTrackCuts(gtrack2)){continue;}

   // Start loop over tracks in the 3rd event:
   for(Int_t iTrack3=0;iTrack3<nTracks3;iTrack3++)
   {
    AliAODTrack *atrack3 = dynamic_cast<AliAODTrack*>(ca3->UncheckedAt(iTrack3));
    // TBI Temporary track insanity checks:
    if(!atrack3){Fatal(sMethodName.Data(),"!atrack3");} // TBI keep this for some time, eventually just continue
    if(atrack3->GetID()>=0 && atrack3->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack3->GetID()>=0 && atrack3->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
    if(atrack3->TestFilterBit(128) && atrack3->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack3->TestFiletrBit(128) && atrack3->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
    if(!PassesCommonTrackCuts(atrack3)){continue;} // TBI re-think
    // Corresponding AOD global track:
    Int_t id3 = atrack3->GetID();
    AliAODTrack *gtrack3 = dynamic_cast<AliAODTrack*>(id3>=0 ? ca3->UncheckedAt(fGlobalTracksAOD[3]->GetValue(id3)) : ca3->UncheckedAt(fGlobalTracksAOD[3]->GetValue(-(id3+1))));
    if(!gtrack3){Fatal(sMethodName.Data(),"!gtrack3");} // TBI keep this for some time, eventually just continue
    // Common track selection criteria for all "normal" global tracks:
    if(!PassesGlobalTrackCuts(gtrack3)){continue;}

    // Okay... So we have three tracks from three different events ready. let's decide whether kinemetics from 'atrack' or 'gtrack' is taken, let's check PID, and calculate the background:
    AliAODTrack *agtrack1 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
    AliAODTrack *agtrack2 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
    AliAODTrack *agtrack3 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
    if(fFillControlHistogramsWithGlobalTrackInfo)
    {
     agtrack1 = gtrack1;
     agtrack2 = gtrack2;
     agtrack3 = gtrack3;
    }
    else
    {
     agtrack1 = atrack1;
     agtrack2 = atrack2;
     agtrack3 = atrack3;
    } // if(fFillControlHistogramsWithGlobalTrackInfo)
    if(!agtrack1){Fatal(sMethodName.Data(),"!agtrack1");}
    if(!agtrack2){Fatal(sMethodName.Data(),"!agtrack2");}
    if(!agtrack3){Fatal(sMethodName.Data(),"!agtrack3");}

    // Cases of interest (needs to be synchronized with the calculations of correlations functions):
    // a) Same species and same charge;
    // b) Same species but different charge combinations (modulo permutations);
    // c) Two pions + something else (modulo permutations);
    // d) Two nucleons + something else (modulo permutations).

    // a) Same species:
    // pi+pi+pi+
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE))
    {
     f3pBackground[2][2][2]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi-pi-pi-
    if(Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE))
    {
     f3pBackground[7][7][7]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // K+K+K+
    if(Kaon(gtrack1,1,kTRUE) && Kaon(gtrack2,1,kTRUE) && Kaon(gtrack3,1,kTRUE))
    {
     f3pBackground[3][3][3]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // K-K-K-
    if(Kaon(gtrack1,-1,kTRUE) && Kaon(gtrack2,-1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pBackground[8][8][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p+p+p+
    if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE) && Proton(gtrack3,1,kTRUE))
    {
     f3pBackground[4][4][4]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p-p-p-
    if(Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE) && Proton(gtrack3,-1,kTRUE))
    {
     f3pBackground[9][9][9]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }

    // b) Same species but different charge combinations (modulo permutations):
    // pi+pi+pi-
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,-1,kTRUE))
    {
     f3pBackground[2][2][7]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi-pi-
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE))
    {
     f3pBackground[2][7][7]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // K+K+K-
    if(Kaon(gtrack1,1,kTRUE) && Kaon(gtrack2,1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pBackground[3][3][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // K+K-K-
    if(Kaon(gtrack1,1,kTRUE) && Kaon(gtrack2,-1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pBackground[3][8][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p+p+p-
    if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE) && Proton(gtrack3,-1,kTRUE))
    {
     f3pBackground[4][4][9]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p+p-p-
    if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,-1,kTRUE) && Proton(gtrack3,-1,kTRUE))
    {
     f3pBackground[4][9][9]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }

    // c) Two pions + something else (modulo permutations):
    // pi+pi+K+
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Kaon(gtrack3,1,kTRUE))
    {
     f3pBackground[2][2][3]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi+K-
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pBackground[2][2][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi-pi-K+
    if(Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Kaon(gtrack3,1,kTRUE))
    {
     f3pBackground[7][7][3]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi-pi-K-
    if(Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pBackground[7][7][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi-K+
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Kaon(gtrack3,1,kTRUE))
    {
     f3pBackground[2][7][3]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi-K-
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pBackground[2][7][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi+p+
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Proton(gtrack3,1,kTRUE))
    {
     f3pBackground[2][2][4]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi+p-
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Proton(gtrack3,-1,kTRUE))
    {
     f3pBackground[2][2][9]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi-pi-p+
    if(Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Proton(gtrack3,1,kTRUE))
    {
     f3pBackground[7][7][4]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi-pi-p-
    if(Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Proton(gtrack3,-1,kTRUE))
    {
     f3pBackground[7][7][9]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi-p+
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Proton(gtrack3,1,kTRUE))
    {
     f3pBackground[2][7][4]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi-p-
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Proton(gtrack3,-1,kTRUE))
    {
     f3pBackground[2][7][9]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }

    // d) Two nucleons + something else (modulo permutations):
    // p+p+pi+
    if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE))
    {
     f3pBackground[4][4][2]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p+p+pi-
    if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE) && Pion(gtrack3,-1,kTRUE))
    {
     f3pBackground[4][4][7]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p+p+K+
    if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE) && Kaon(gtrack3,1,kTRUE))
    {
     f3pBackground[4][4][3]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p+p+K-
    if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pBackground[4][4][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p-p-pi+
    if(Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE) && Pion(gtrack3,1,kTRUE))
    {
     f3pBackground[9][9][2]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p-p-pi-
    if(Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE))
    {
     f3pBackground[9][9][7]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p-p-K+
    if(Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE) && Kaon(gtrack3,1,kTRUE))
    {
     f3pBackground[9][9][3]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p-p-K-
    if(Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pBackground[9][9][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }

   } // for(Int_t iTrack3=0;iTrack3<nTracks3;iTrack3++)
  } // for(Int_t iTrack2=0;iTrack2<nTracks2;iTrack2++)
 } // for(Int_t iTrack1=0;iTrack1<nTracks1;iTrack1++)

} // void AliAnalysisTaskMultiparticleFemtoscopy::Calculate3pBackground(TClonesArray *ca1, TClonesArray *ca2, TClonesArray *ca3)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::Calculate3pBackground(TClonesArray *ca1, TClonesArray *ca2, TClonesArray *ca3, TExMap *em1, TExMap *em2, TExMap *em3)
{
 // Calculate background.

 // a) Insanity checks;
 // b) Three nested loops to calculate B(k);

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::Calculate3pBackground(TClonesArray *ca1, TClonesArray *ca2, TClonesArray *ca3, TExMap *em1, TExMap *em2, TExMap *em3)";
 if(!ca1){Fatal(sMethodName.Data(),"!ca1");}
 if(!ca2){Fatal(sMethodName.Data(),"!ca2");}
 if(!ca3){Fatal(sMethodName.Data(),"!ca3");}
 if(!em1){Fatal(sMethodName.Data(),"!em1");}
 if(!em2){Fatal(sMethodName.Data(),"!em2");}
 if(!em3){Fatal(sMethodName.Data(),"!em3");}

 // b) Three nested loops to calculate B(k):
 Int_t nTracks1 = ca1->GetEntries();
 Int_t nTracks2 = ca2->GetEntries();
 Int_t nTracks3 = ca3->GetEntries();

 // Start loop over tracks in the 1st event:
 for(Int_t iTrack1=0;iTrack1<nTracks1;iTrack1++)
 {
  AliAODTrack *atrack1 = dynamic_cast<AliAODTrack*>(ca1->UncheckedAt(iTrack1)); // TBI cross-check UncheckedAt
  // TBI Temporary track insanity checks:
  if(!atrack1){Fatal(sMethodName.Data(),"!atrack1");} // TBI keep this for some time, eventually just continue
  if(atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(atrack1->TestFilterBit(128) && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->TestFiletrBit(128) && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(!PassesCommonTrackCuts(atrack1)){continue;} // TBI re-think
  // Corresponding AOD global track:
  Int_t id1 = atrack1->GetID();
  AliAODTrack *gtrack1 = dynamic_cast<AliAODTrack*>(id1>=0 ? ca1->UncheckedAt(em1->GetValue(id1)) : ca1->UncheckedAt(em1->GetValue(-(id1+1))));
  if(!gtrack1){Fatal(sMethodName.Data(),"!gtrack1");} // TBI keep this for some time, eventually just continue
  // Common track selection criteria for all "normal" global tracks:
  if(!PassesGlobalTrackCuts(gtrack1)){continue;}

  // Start loop over tracks in the 2nd event:
  for(Int_t iTrack2=0;iTrack2<nTracks2;iTrack2++)
  {
   AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(ca2->UncheckedAt(iTrack2));
   // TBI Temporary track insanity checks:
   if(!atrack2){Fatal(sMethodName.Data(),"!atrack2");} // TBI keep this for some time, eventually just continue
   if(atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(atrack2->TestFilterBit(128) && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->TestFiletrBit(128) && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(!PassesCommonTrackCuts(atrack2)){continue;} // TBI re-think
   // Corresponding AOD global track:
   Int_t id2 = atrack2->GetID();
   AliAODTrack *gtrack2 = dynamic_cast<AliAODTrack*>(id2>=0 ? ca2->UncheckedAt(em2->GetValue(id2)) : ca2->UncheckedAt(em2->GetValue(-(id2+1))));
   if(!gtrack2){Fatal(sMethodName.Data(),"!gtrack2");} // TBI keep this for some time, eventually just continue
   // Common track selection criteria for all "normal" global tracks:
   if(!PassesGlobalTrackCuts(gtrack2)){continue;}

   // Start loop over tracks in the 3rd event:
   for(Int_t iTrack3=0;iTrack3<nTracks3;iTrack3++)
   {
    AliAODTrack *atrack3 = dynamic_cast<AliAODTrack*>(ca3->UncheckedAt(iTrack3));
    // TBI Temporary track insanity checks:
    if(!atrack3){Fatal(sMethodName.Data(),"!atrack3");} // TBI keep this for some time, eventually just continue
    if(atrack3->GetID()>=0 && atrack3->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack3->GetID()>=0 && atrack3->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
    if(atrack3->TestFilterBit(128) && atrack3->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack3->TestFiletrBit(128) && atrack3->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
    if(!PassesCommonTrackCuts(atrack3)){continue;} // TBI re-think
    // Corresponding AOD global track:
    Int_t id3 = atrack3->GetID();
    AliAODTrack *gtrack3 = dynamic_cast<AliAODTrack*>(id3>=0 ? ca3->UncheckedAt(em3->GetValue(id3)) : ca3->UncheckedAt(em3->GetValue(-(id3+1))));
    if(!gtrack3){Fatal(sMethodName.Data(),"!gtrack3");} // TBI keep this for some time, eventually just continue
    // Common track selection criteria for all "normal" global tracks:
    if(!PassesGlobalTrackCuts(gtrack3)){continue;}

    // Okay... So we have three tracks from three different events ready. let's decide whether kinemetics from 'atrack' or 'gtrack' is taken, let's check PID, and calculate the background:
    AliAODTrack *agtrack1 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
    AliAODTrack *agtrack2 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
    AliAODTrack *agtrack3 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
    if(fFillControlHistogramsWithGlobalTrackInfo)
    {
     agtrack1 = gtrack1;
     agtrack2 = gtrack2;
     agtrack3 = gtrack3;
    }
    else
    {
     agtrack1 = atrack1;
     agtrack2 = atrack2;
     agtrack3 = atrack3;
    } // if(fFillControlHistogramsWithGlobalTrackInfo)
    if(!agtrack1){Fatal(sMethodName.Data(),"!agtrack1");}
    if(!agtrack2){Fatal(sMethodName.Data(),"!agtrack2");}
    if(!agtrack3){Fatal(sMethodName.Data(),"!agtrack3");}

    // Cases of interest (needs to be synchronized with the calculations of correlations functions):
    // a) Same species and same charge;
    // b) Same species but different charge combinations (modulo permutations);
    // c) Two pions + something else (modulo permutations);
    // d) Two nucleons + something else (modulo permutations).

    // a) Same species:
    // pi+pi+pi+
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE))
    {
     f3pBackground[2][2][2]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi-pi-pi-
    if(Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE))
    {
     f3pBackground[7][7][7]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // K+K+K+
    if(Kaon(gtrack1,1,kTRUE) && Kaon(gtrack2,1,kTRUE) && Kaon(gtrack3,1,kTRUE))
    {
     f3pBackground[3][3][3]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // K-K-K-
    if(Kaon(gtrack1,-1,kTRUE) && Kaon(gtrack2,-1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pBackground[8][8][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p+p+p+
    if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE) && Proton(gtrack3,1,kTRUE))
    {
     f3pBackground[4][4][4]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p-p-p-
    if(Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE) && Proton(gtrack3,-1,kTRUE))
    {
     f3pBackground[9][9][9]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }

    // b) Same species but different charge combinations (modulo permutations):
    // pi+pi+pi-
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,-1,kTRUE))
    {
     f3pBackground[2][2][7]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi-pi-
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE))
    {
     f3pBackground[2][7][7]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // K+K+K-
    if(Kaon(gtrack1,1,kTRUE) && Kaon(gtrack2,1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pBackground[3][3][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // K+K-K-
    if(Kaon(gtrack1,1,kTRUE) && Kaon(gtrack2,-1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pBackground[3][8][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p+p+p-
    if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE) && Proton(gtrack3,-1,kTRUE))
    {
     f3pBackground[4][4][9]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p+p-p-
    if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,-1,kTRUE) && Proton(gtrack3,-1,kTRUE))
    {
     f3pBackground[4][9][9]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }

    // c) Two pions + something else (modulo permutations):
    // pi+pi+K+
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Kaon(gtrack3,1,kTRUE))
    {
     f3pBackground[2][2][3]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi+K-
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pBackground[2][2][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi-pi-K+
    if(Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Kaon(gtrack3,1,kTRUE))
    {
     f3pBackground[7][7][3]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi-pi-K-
    if(Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pBackground[7][7][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi-K+
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Kaon(gtrack3,1,kTRUE))
    {
     f3pBackground[2][7][3]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi-K-
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pBackground[2][7][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi+p+
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Proton(gtrack3,1,kTRUE))
    {
     f3pBackground[2][2][4]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi+p-
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Proton(gtrack3,-1,kTRUE))
    {
     f3pBackground[2][2][9]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi-pi-p+
    if(Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Proton(gtrack3,1,kTRUE))
    {
     f3pBackground[7][7][4]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi-pi-p-
    if(Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Proton(gtrack3,-1,kTRUE))
    {
     f3pBackground[7][7][9]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi-p+
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Proton(gtrack3,1,kTRUE))
    {
     f3pBackground[2][7][4]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // pi+pi-p-
    if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Proton(gtrack3,-1,kTRUE))
    {
     f3pBackground[2][7][9]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }

    // d) Two nucleons + something else (modulo permutations):
    // p+p+pi+
    if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE))
    {
     f3pBackground[4][4][2]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p+p+pi-
    if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE) && Pion(gtrack3,-1,kTRUE))
    {
     f3pBackground[4][4][7]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p+p+K+
    if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE) && Kaon(gtrack3,1,kTRUE))
    {
     f3pBackground[4][4][3]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p+p+K-
    if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pBackground[4][4][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p-p-pi+
    if(Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE) && Pion(gtrack3,1,kTRUE))
    {
     f3pBackground[9][9][2]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p-p-pi-
    if(Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE))
    {
     f3pBackground[9][9][7]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p-p-K+
    if(Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE) && Kaon(gtrack3,1,kTRUE))
    {
     f3pBackground[9][9][3]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }
    // p-p-K-
    if(Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE) && Kaon(gtrack3,-1,kTRUE))
    {
     f3pBackground[9][9][8]->Fill(Q3(agtrack1,agtrack2,agtrack3));
    }

   } // for(Int_t iTrack3=0;iTrack3<nTracks3;iTrack3++)
  } // for(Int_t iTrack2=0;iTrack2<nTracks2;iTrack2++)
 } // for(Int_t iTrack1=0;iTrack1<nTracks1;iTrack1++)

} // void AliAnalysisTaskMultiparticleFemtoscopy::Calculate3pBackground(TClonesArray *ca1, TClonesArray *ca2, TClonesArray *ca3)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::Calculate4pBackground(TClonesArray *ca1, TClonesArray *ca2, TClonesArray *ca3, TClonesArray *ca4)
{
 // TBI

 return;

} // void AliAnalysisTaskMultiparticleFemtoscopy::Calculate4pBackground(TClonesArray *ca1, TClonesArray *ca2, TClonesArray *ca3, TClonesArray *ca4)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::Calculate2pBackground(TClonesArray *ca1, TClonesArray *ca2, Bool_t bMC)
{
 // Calculate background for MC events.

 return; // TBI needs to be re-validated. Shifting is nowadays done elsewhere.

 // TBI this method is not really validated
 // TBI unify with previous method

 // a) Insanity checks;
 // b) Temporary primitive algorithm: correlate particles from 'previous event' to particles from 'current event';
 // c) Two nested loops to calculate B(k);
 // d) Shift [1] -> [0];
 // e) Clean [1].

 TString pattern = ".11.13.211.321.2212."; // TBI this can be done programatically TBI only for these particles I need to estimate background TBI move this somewhere else

 // a) Insanity checks:
 if(!bMC) return; // TBI this is not really needed...
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::Calculate2pBackground(TClonesArray *ca1, TClonesArray *ca2, Bool_t bMC)";
 if(!ca1){Fatal(sMethodName.Data(),"!ca1");}
 if(!ca2){Fatal(sMethodName.Data(),"!ca2");}

 // b) Temporary primitive algorithm: correlate particles from 'previous event' to particles from 'current event'
 // ...

 // c) Two nested loops to calculate B(k):
 Int_t nTracks1 = ca1->GetEntries();
 Int_t nTracks2 = ca2->GetEntries();

 // Start loop over tracks in the 1st event:
 for(Int_t iTrack1=0;iTrack1<nTracks1;iTrack1++)
 {
  AliAODMCParticle *amcparticle1 = dynamic_cast<AliAODMCParticle*>(ca1->UncheckedAt(iTrack1)); // TBI cross-check UncheckedAt
  // TBI Temporary track insanity checks:
  if(!amcparticle1){Fatal(sMethodName.Data(),"!amcparticle1");} // TBI keep this for some time, eventually just continue
  if(!PassesCommonTrackCuts(amcparticle1)){continue;} // TBI re-think

  // Start loop over tracks in the 2nd event:
  for(Int_t iTrack2=0;iTrack2<nTracks2;iTrack2++)
  {
   AliAODMCParticle *amcparticle2 = dynamic_cast<AliAODMCParticle*>(ca2->UncheckedAt(iTrack2)); // TBI cross-check UncheckedAt
   // TBI Temporary track insanity checks:
   if(!amcparticle2){Fatal(sMethodName.Data(),"!amcparticle2");} // TBI keep this for some time, eventually just continue
   if(!PassesCommonTrackCuts(amcparticle2)){continue;} // TBI re-think

   // Okay... So we have two tracks from two different events ready. Let's check PID, and calculate the background:
   if(!(pattern.Contains(Form(".%d.",TMath::Abs(amcparticle1->GetPdgCode()))) && pattern.Contains(Form(".%d.",TMath::Abs(amcparticle2->GetPdgCode()))))){continue;}

   // Determine the indices of correlation function to be filled:
   Int_t index1 = -44;
   Int_t index2 = -44;
   if(fCorrelationFunctionsIndices->GetValue(amcparticle1->GetPdgCode())<=fCorrelationFunctionsIndices->GetValue(amcparticle2->GetPdgCode()))
   {
    index1 = fCorrelationFunctionsIndices->GetValue(amcparticle1->GetPdgCode());
    index2 = fCorrelationFunctionsIndices->GetValue(amcparticle2->GetPdgCode());
   }
   else
   {
    index1 = fCorrelationFunctionsIndices->GetValue(amcparticle2->GetPdgCode());
    index2 = fCorrelationFunctionsIndices->GetValue(amcparticle1->GetPdgCode());
   }
   if(-44 == index1){Fatal(sMethodName.Data(),"-44 == index1");}
   if(-44 == index2){Fatal(sMethodName.Data(),"-44 == index2");}

   // Fill or die:
   f2pBackground[index1][index2]->Fill(RelativeMomenta(amcparticle1,amcparticle2)); // for the relative momenta ordering doesn't matter

  } // for(Int_t iTrack2=0;iTrack2<nTracks2;iTrack2++)
 } // for(Int_t iTrack1=0;iTrack1<nTracks1;iTrack1++)

 // d) Shift [1] -> [0]:
 // TBI re-think the lines below, there should be a better way...
 fMixedEvents0[0] = (TClonesArray*)fMixedEvents0[1]->Clone();

 // e) Clean [1]:
 fMixedEvents0[1] = NULL;

} // void AliAnalysisTaskMultiparticleFemtoscopy::Calculate2pBackground(TClonesArray *ca1, TClonesArray *ca2, Bool_t bMC)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::GetOutputHistograms(TList *histList)
{
 // Get pointers for everything in the base list "fHistList".

 // a) Get pointer for base list fHistList;
 // b) Get pointer for profile holding internal flags and set again all flags (TBI this profile is not implemented yet!!!!);
 // *) TBI ...
 // *) Get pointers for correlation functions;
 // *) Get pointers for background;
 // *) Get pointers for buffers.

 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::GetOutputHistograms(TList *histList)";

 // a) Get pointer for base list fHistList and profile holding internal flags;
 fHistList = histList;
 if(!fHistList){Fatal(sMethodName.Data(),"fHistList is not around today...");}

 // b) Get pointer for profile holding internal flags and set again all flags TBI this profile is not implemented yet
 //fInternalFlagsPro = dynamic_cast<TProfile*>(fHistList->FindObject("fInternalFlagsPro")); TBI this was example from MPC
 //if(!fInternalFlagsPro){Fatal(sMethodName.Data(),"fInternalFlagsPro");} TBI this was example from MPC

 // *) TBI
 //this->GetPointersFor...;

 // *) Get pointers for correlation functions:
 this->GetPointersForCorrelationFunctions();

 // *) Get pointers for background:
 this->GetPointersForBackground();

 // *) Get pointers for buffers:
 this->GetPointersForBuffers();

} // void AliAnalysisTaskMultiparticleFemtoscopy::GetOutputHistograms(TList *histList)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::GetPointersForCorrelationFunctions()
{
 // Get pointers for correlation functions.

 // a) Get pointer for fCorrelationFunctionsList and fCorrelationFunctionsSublist[3];
 // b) Get pointer for fCorrelationFunctionsFlagsPro;
 // c) Set again all flags;
 // d) Get pointers for fCorrelationFunctions[10][10];
 // e) Get pointers for f3pCorrelationFunctions[10][10][10];
 // f) Get pointers for f4pCorrelationFunctions[10][10][10][10].

 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::GetPointersForCorrelationFunctions()";

 // a) Get pointer for fCorrelationFunctionsList and fCorrelationFunctionsSublist[3]:
 fCorrelationFunctionsList = dynamic_cast<TList*>(fHistList->FindObject("CorrelationFunctions"));
 if(!fCorrelationFunctionsList){Fatal(sMethodName.Data(),"fCorrelationFunctionsList");}
 for(Int_t cfs=0;cfs<3;cfs++)
 {
  fCorrelationFunctionsSublist[cfs] = dynamic_cast<TList*>(fCorrelationFunctionsList->FindObject(Form("f%dpCorrelationFunctions",cfs+2)));
  if(!fCorrelationFunctionsSublist[cfs]){Fatal(sMethodName.Data(),"fCorrelationFunctionsSublist[%d]",cfs);}
 } // for(Int_t cfs=0;cfs<3;cfs++)

 // b) Get pointer for fCorrelationFunctionsFlagsPro:
 fCorrelationFunctionsFlagsPro = dynamic_cast<TProfile*>(fCorrelationFunctionsList->FindObject("fCorrelationFunctionsFlagsPro"));
 if(!fCorrelationFunctionsFlagsPro){Fatal(sMethodName.Data(),"fCorrelationFunctionsFlagsPro");}
 f2pCorrelationFunctionsFlagsPro = dynamic_cast<TProfile*>(fCorrelationFunctionsSublist[0]->FindObject("f2pCorrelationFunctionsFlagsPro"));
 if(!f2pCorrelationFunctionsFlagsPro){Fatal(sMethodName.Data(),"f2pCorrelationFunctionsFlagsPro");}
 f3pCorrelationFunctionsFlagsPro = dynamic_cast<TProfile*>(fCorrelationFunctionsSublist[1]->FindObject("f3pCorrelationFunctionsFlagsPro"));
 if(!f3pCorrelationFunctionsFlagsPro){Fatal(sMethodName.Data(),"f3pCorrelationFunctionsFlagsPro");}
 f4pCorrelationFunctionsFlagsPro = dynamic_cast<TProfile*>(fCorrelationFunctionsSublist[2]->FindObject("f4pCorrelationFunctionsFlagsPro"));
 if(!f4pCorrelationFunctionsFlagsPro){Fatal(sMethodName.Data(),"f4pCorrelationFunctionsFlagsPro");}

 // c) Set again all flags:
 fFillCorrelationFunctions = (Bool_t) fCorrelationFunctionsFlagsPro->GetBinContent(1);
 fNormalizeCorrelationFunctions = (Bool_t) fCorrelationFunctionsFlagsPro->GetBinContent(2);
 fFill3pCorrelationFunctions = (Bool_t) fCorrelationFunctionsFlagsPro->GetBinContent(3);
 fFill4pCorrelationFunctions = (Bool_t) fCorrelationFunctionsFlagsPro->GetBinContent(4);
 fNormalizationOption = (Int_t) fCorrelationFunctionsFlagsPro->GetBinContent(5);
 fNormalizationInterval[0] = (Float_t) fCorrelationFunctionsFlagsPro->GetBinContent(6);
 fNormalizationInterval[1] = (Float_t) fCorrelationFunctionsFlagsPro->GetBinContent(7);
 fnMergedBins = (Int_t) fCorrelationFunctionsFlagsPro->GetBinContent(8);

 if(!fFillCorrelationFunctions){return;} // TBI is this safe enough

 // d) Get pointers for fCorrelationFunctions[10][10]:
 const Int_t nParticleSpecies = 5; // Supported at the moment: 0=e,1=mu,2=pi,3=K,4=p
 for(Int_t pid1=0;pid1<2*nParticleSpecies;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pid2=pid1;pid2<2*nParticleSpecies;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   fCorrelationFunctions[pid1][pid2] = dynamic_cast<TH1F*>(fCorrelationFunctionsSublist[0]->FindObject(Form("fCorrelationFunctions[%d][%d]",pid1,pid2)));
   if(!fCorrelationFunctions[pid1][pid2]){Fatal(sMethodName.Data(),"fCorrelationFunctions[%d][%d]",pid1,pid2);}
  } // for(Int_t pid2=0;pid2<5;pid2++) // anti-particle [0=e,1=mu,2=pi,3=K,4=p]
 } // for(Int_t pid=0;pid<5;pid++) // particle [0=e,1=mu,2=pi,3=K,4=p]

 // e) Get pointers for f3pCorrelationFunctions[10][10][10]:
 if(fFill3pCorrelationFunctions)
 {
  for(Int_t pid1=0;pid1<2*nParticleSpecies;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   for(Int_t pid2=0;pid2<2*nParticleSpecies;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
   {
    for(Int_t pid3=0;pid3<2*nParticleSpecies;pid3++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
    {
     f3pCorrelationFunctions[pid1][pid2][pid3] = dynamic_cast<TH1F*>(fCorrelationFunctionsSublist[1]->FindObject(Form("f3pCorrelationFunctions[%d][%d][%d]",pid1,pid2,pid3)));
     if(!f3pCorrelationFunctions[pid1][pid2][pid3]){Fatal(sMethodName.Data(),"f3pCorrelationFunctions[%d][%d][%d]",pid1,pid2,pid3);}
    } // for(Int_t pid3=0;pid3<2*nParticleSpecies;pid3++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
   } // for(Int_t pid2=0;pid2<5;pid2++) // anti-particle [0=e,1=mu,2=pi,3=K,4=p]
  } // for(Int_t pid=0;pid<5;pid++) // particle [0=e,1=mu,2=pi,3=K,4=p]
 } // if(fFill3pCorrelationFunctions)

 // f) Get pointers for f4pCorrelationFunctions[10][10][10][10]:
 if(fFill4pCorrelationFunctions)
 {
  for(Int_t pid1=0;pid1<2*nParticleSpecies;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   for(Int_t pid2=0;pid2<2*nParticleSpecies;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
   {
    for(Int_t pid3=0;pid3<2*nParticleSpecies;pid3++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
    {
     for(Int_t pid4=0;pid4<2*nParticleSpecies;pid4++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
     {
      f4pCorrelationFunctions[pid1][pid2][pid3][pid4] = dynamic_cast<TH1F*>(fCorrelationFunctionsSublist[2]->FindObject(Form("f4pCorrelationFunctions[%d][%d][%d][%d]",pid1,pid2,pid3,pid4)));
      if(!f4pCorrelationFunctions[pid1][pid2][pid3][pid4]){Fatal(sMethodName.Data(),"f4pCorrelationFunctions[%d][%d][%d][%d]",pid1,pid2,pid3,pid4);}
     } // for(Int_t pid4=0;pid4<2*nParticleSpecies;pid4++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
    } // for(Int_t pid3=0;pid3<2*nParticleSpecies;pid3++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
   } // for(Int_t pid2=0;pid2<2*nParticleSpecies;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  } // for(Int_t pid1=0;pid1<2*nParticleSpecies;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 } // if(fFill4pCorrelationFunctions)

} // void AliAnalysisTaskMultiparticleFemtoscopy::GetPointersForCorrelationFunctions()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::GetPointersForBackground()
{
 // Get pointers for background.

 // a) Get pointer for fBackgroundList and fBackgroundSublist[3];
 // b) Get pointer for fBackgroundFlagsPro;
 // c) Set again all flags;
 // d) Get pointers for f2pBackground[10][10];
 // e) Get pointers for f3pBackground[10][10][10].
 // f) TBI 4p

 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::GetPointersForBackground()";

 // a) Get pointer for fBackgroundList and fBackgroundSubist[3];
 fBackgroundList = dynamic_cast<TList*>(fHistList->FindObject("Background"));
 if(!fBackgroundList){Fatal(sMethodName.Data(),"fBackgroundList");}
 for(Int_t bs=0;bs<3;bs++)
 {
  fBackgroundSublist[bs] = dynamic_cast<TList*>(fBackgroundList->FindObject(Form("f%dpBackground",bs+2)));
  if(!fBackgroundSublist[bs]){Fatal(sMethodName.Data(),"fBackgroundSublist[%d]",bs);}
 } // for(Int_t bs=0;bs<3;bs++)

 // b) Get pointer for fBackgroundFlagsPro:
 fBackgroundFlagsPro = dynamic_cast<TProfile*>(fBackgroundList->FindObject("fBackgroundFlagsPro"));
 if(!fBackgroundFlagsPro){Fatal(sMethodName.Data(),"fBackgroundFlagsPro");}
 f2pBackgroundFlagsPro = dynamic_cast<TProfile*>(fBackgroundSublist[0]->FindObject("f2pBackgroundFlagsPro"));
 if(!f2pBackgroundFlagsPro){Fatal(sMethodName.Data(),"f2pBackgroundFlagsPro");}
 f3pBackgroundFlagsPro = dynamic_cast<TProfile*>(fBackgroundSublist[1]->FindObject("f3pBackgroundFlagsPro"));
 if(!f3pBackgroundFlagsPro){Fatal(sMethodName.Data(),"f3pBackgroundFlagsPro");}
 f4pBackgroundFlagsPro = dynamic_cast<TProfile*>(fBackgroundSublist[2]->FindObject("f4pBackgroundFlagsPro"));
 if(!f4pBackgroundFlagsPro){Fatal(sMethodName.Data(),"f4pBackgroundFlagsPro");}

 // c) Set again all flags:
 fBackgroundOption = (Int_t) fBackgroundFlagsPro->GetBinContent(1);
 fEstimate2pBackground = (Bool_t) fBackgroundFlagsPro->GetBinContent(2);
 fEstimate3pBackground = (Bool_t) fBackgroundFlagsPro->GetBinContent(3);
 fEstimate4pBackground = (Bool_t) fBackgroundFlagsPro->GetBinContent(4);

 if(! (fEstimate2pBackground || fEstimate3pBackground || fEstimate4pBackground )){return;} // TBI is this safe enough

 // d) Get pointers for f2pBackground[10][10]:
 const Int_t nParticleSpecies = 5; // Supported at the moment: 0=e,1=mu,2=pi,3=K,4=p
 for(Int_t pid1=0;pid1<2*nParticleSpecies;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pid2=pid1;pid2<2*nParticleSpecies;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   f2pBackground[pid1][pid2] = dynamic_cast<TH1F*>(fBackgroundSublist[0]->FindObject(Form("f2pBackground[%d][%d]",pid1,pid2)));
   if(!f2pBackground[pid1][pid2]){Fatal(sMethodName.Data(),"f2pBackground[%d][%d]",pid1,pid2);}
  } // for(Int_t pid2=0;pid2<5;pid2++) // anti-particle [0=e,1=mu,2=pi,3=K,4=p]
 } // for(Int_t pid=0;pid<5;pid++) // particle [0=e,1=mu,2=pi,3=K,4=p]

 // e) Get pointers for f3pBackground[10][10][10]:
 for(Int_t pid1=0;pid1<2*nParticleSpecies;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pid2=0;pid2<2*nParticleSpecies;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   for(Int_t pid3=0;pid3<2*nParticleSpecies;pid3++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
   {
    f3pBackground[pid1][pid2][pid3] = dynamic_cast<TH1F*>(fBackgroundSublist[1]->FindObject(Form("f3pBackground[%d][%d][%d]",pid1,pid2,pid3)));
    if(!f3pBackground[pid1][pid2][pid3]){Fatal(sMethodName.Data(),"f3pBackground[%d][%d][%d]",pid1,pid2,pid3);}
   }
  } // for(Int_t pid2=0;pid2<5;pid2++) // anti-particle [0=e,1=mu,2=pi,3=K,4=p]
 } // for(Int_t pid=0;pid<5;pid++) // particle [0=e,1=mu,2=pi,3=K,4=p]

 // f) TBI get pointers for 4p background

} // void AliAnalysisTaskMultiparticleFemtoscopy::GetPointersForBackground()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::GetPointersForBuffers()
{
 // Get pointers for buffers.

 // a) Get pointer for fBuffersList;
 // b) Get pointer for fBuffersFlagsPro;
 // c) Set again all flags;
 // d) ...

 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::GetPointersForBuffers()";

 // a) Get pointer for fBuffersList:
 fBuffersList = dynamic_cast<TList*>(fHistList->FindObject("Buffers"));
 if(!fBuffersList){Fatal(sMethodName.Data(),"fBuffersList");}

 // b) Get pointer for fBuffersFlagsPro:
 fBuffersFlagsPro = dynamic_cast<TProfile*>(fBuffersList->FindObject("fBuffersFlagsPro"));
 if(!fBuffersFlagsPro){Fatal(sMethodName.Data(),"fBuffersFlagsPro");}

 // c) Set again all flags:
 fFillBuffers = (Bool_t) fBuffersFlagsPro->GetBinContent(1);
 fMaxBuffer = fBuffersFlagsPro->GetBinContent(2);

} // void AliAnalysisTaskMultiparticleFemtoscopy::GetPointersForBuffers()

//=======================================================================================================================

Int_t AliAnalysisTaskMultiparticleFemtoscopy::InsanityChecksForTracks(AliAODTrack *atrack)
{
 // Insanity checks for each track ('atrack') in AOD.

 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::InsanityChecksForTracks(AliVEvent *ave)";

 // List of insanity checks (enumeration corresponds to return values):
 // 1) Check whether 'normal global' track is also 'global constrained';
 // 2) Check whether 'TPC only' track is also 'global constrained'.

 // 1) Check whether 'normal global' track is also 'global constrained':
 if(atrack->GetID()>=0 && atrack->IsGlobalConstrained()) return 1;

 // 2) Check whether 'TPC only' track is also 'global constrained':
 if(atrack->TestFilterBit(128) && atrack->IsGlobalConstrained()) return 2;

 // TBI work out case below:
   // if(0==atrack->GetFilterMap())
   // TBI a vast majority of such tracks pass
   //   if(0==atrack->GetFilterMap() && atrack->GetType() == AliAODTrack::kFromDecayVtx)
   // but there are exceptions which pass
   //   if(0==atrack->GetFilterMap() && atrack->GetType() == AliAODTrack::kPrimary)
   // So it's a mess...
   // 1) Clearly, for none of such tracks we can test a filter bit
   // 2) Apparently all of them have positive ID, meaning they are global (?)

 return 0;

} // Bool_t AliAnalysisTaskMultiparticleFemtoscopy::InsanityChecksForTracks(AliAODTrack *atrack

//=======================================================================================================================

Int_t AliAnalysisTaskMultiparticleFemtoscopy::InsanityChecksForGlobalTracks(AliAODTrack *atrack)
{
 // Insanity checks only for global tracks ('gtrack') in AOD.

 // List of insanity checks (enumeration corresponds to return values):
 // 1) Wrong ID signature, for 'normal global' tracks ID >=0; TBI this check is most likely an overkill, drop it at some point
 // 2) ...

 // 1) Check the signature of ID;
 if(atrack->GetID()<0) return 1;

 return 0;

} // Bool_t AliAnalysisTaskMultiparticleFemtoscopy::InsanityChecksForGlobalTracks(AliAODTrack *gtrack)

//=======================================================================================================================

Int_t AliAnalysisTaskMultiparticleFemtoscopy::BinNoForSpecifiedValue(TH1F *hist, Double_t value)
{
 // Determine in which bin the specified value would lend. If the call fails, returns -1.

 if(!hist){return -1;}

 Int_t nBins = hist->GetXaxis()->GetNbins();
 Double_t min = hist->GetBinLowEdge(1);
 Double_t max = hist->GetBinLowEdge(nBins+1);
 Double_t bw = (max-min)/nBins; // bin width
 Int_t binNo = 1;
 for(Int_t b=1;b<=nBins;b++)
 {
  if(min+b*bw>value){binNo=b; break;}
 }
 if(value<min||value>=max){binNo=-1;} // underflow or overlflow

 return binNo;

} // Int_t AliAnalysisTaskMultiparticleFemtoscopy::BinNoForSpecifiedValue(TH1D *hist, Double_t value)

//=======================================================================================================================

Int_t AliAnalysisTaskMultiparticleFemtoscopy::BinNoForSpecifiedValue(TProfile *pro, Double_t value)
{
 // Determine in which bin the specified value would lend. If the call fails, returns -1.

 if(!pro){return -1;}

 Int_t nBins = pro->GetXaxis()->GetNbins();
 Double_t min = pro->GetBinLowEdge(1);
 Double_t max = pro->GetBinLowEdge(nBins+1);
 Double_t bw = (max-min)/nBins; // bin width
 Int_t binNo = 1;
 for(Int_t b=1;b<=nBins;b++)
 {
  if(min+b*bw>value){binNo=b; break;}
 }
 if(value<min||value>=max){binNo=-1;} // underflow or overlflow

 return binNo;

} // Int_t AliAnalysisTaskMultiparticleFemtoscopy::BinNoForSpecifiedValue(TH1D *hist, Double_t value)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::HybridApproach1stTerm(AliAODEvent *aAOD)
{
 // Calculate 1st term in the 'hybrid approach'.

 // a) Insanity checks;
 // b) Nested loops to calculate the 3p yield.

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::HybridApproach1stTerm(AliAODEvent *aAOD)";
 if(!aAOD){Fatal(sMethodName.Data(),"!aAOD");}
 if(0 == fGlobalTracksAOD[0]->GetSize()){Fatal(sMethodName.Data(),"0 == fGlobalTracksAOD[0]->GetSize()");} // this case shall be already treated in UserExec by default, so I can re-use it here

 // b) Nested loops to calculate the 3p yield:
 Int_t nTracks = aAOD->GetNumberOfTracks();
 for(Int_t iTrack1=0;iTrack1<nTracks;iTrack1++)
 {
  AliAODTrack *atrack1 = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack1));
  // TBI Temporary track insanity checks:
  if(!atrack1){Fatal(sMethodName.Data(),"!atrack1");} // TBI keep this for some time, eventually just continue
  if(atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(atrack1->TestFilterBit(128) && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->TestFiletrBit(128) && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(!PassesCommonTrackCuts(atrack1)){continue;} // TBI re-think
  // Corresponding AOD global track:
  Int_t id1 = atrack1->GetID();
  AliAODTrack *gtrack1 = dynamic_cast<AliAODTrack*>(id1>=0 ? aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(id1)) : aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(-(id1+1))));
  if(!gtrack1){Fatal(sMethodName.Data(),"!gtrack1");} // TBI keep this for some time, eventually just continue
  Int_t gid1 = (id1>=0 ? id1 : -(id1+1)); // ID of corresponding global track
  // Common track selection criteria for all "normal" global tracks:
  if(!PassesGlobalTrackCuts(gtrack1)){continue;}

  // Loop over the 2nd particle:
  for(Int_t iTrack2=0;iTrack2<nTracks;iTrack2++) // TBI: I can gain in performance if I do not start from 0
  {
   if(iTrack2<=iTrack1){continue;} // Eliminate self-evident self-correlations, and permutations as well
   //if(iTrack2==iTrack1){continue;} // Eliminate self-evident self-correlations
   AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack2));
   // TBI Temporary track insanity checks:
   if(!atrack2){Fatal(sMethodName.Data(),"!atrack2");} // TBI keep this for some time, eventually just continue
   if(atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(atrack2->TestFilterBit(128) && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->TestFiletrBit(128) && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(!PassesCommonTrackCuts(atrack2)){continue;} // TBI re-think
   // Corresponding AOD global track:
   Int_t id2 = atrack2->GetID();
   AliAODTrack *gtrack2 = dynamic_cast<AliAODTrack*>(id2>=0 ? aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(id2)) : aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(-(id2+1))));
   if(!gtrack2){Fatal(sMethodName.Data(),"!gtrack2");} // TBI keep this for some time, eventually just continue
   Int_t gid2 = (id2>=0 ? id2 : -(id2+1)); // ID of corresponding global track
   if(gid1==gid2){continue;} // Eliminate self-correlations:
   // Common track selection criteria for all "normal" global tracks:
   if(!PassesGlobalTrackCuts(gtrack2)){continue;}

   // Loop over the 3rd particle:
   for(Int_t iTrack3=0;iTrack3<nTracks;iTrack3++)
   {
    if(!fFill3pCorrelationFunctions){break;}
    if(iTrack3<=iTrack2 || iTrack3<=iTrack1){continue;} // Eliminate self-evident self-correlations, and permutations as well
    //if(iTrack3==iTrack2 || iTrack3==iTrack1){continue;} // Eliminate self-evident self-correlations
    AliAODTrack *atrack3 = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack3));
    if(!atrack3){Fatal(sMethodName.Data(),"!atrack3");} // TBI keep this for some time, eventually just continue
    if(atrack3->GetID()>=0 && atrack3->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack3->GetID()>=0 && atrack3->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
    if(atrack3->TestFilterBit(128) && atrack3->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack3->TestFiletrBit(128) && atrack3->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
    if(!PassesCommonTrackCuts(atrack3)){continue;} // TBI re-think
    // Corresponding AOD global track:
    Int_t id3 = atrack3->GetID();
    AliAODTrack *gtrack3 = dynamic_cast<AliAODTrack*>(id3>=0 ? aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(id3)) : aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(-(id3+1))));
    if(!gtrack3){Fatal(sMethodName.Data(),"!gtrack3");} // TBI keep this for some time, eventually just continue
    Int_t gid3 = (id3>=0 ? id3 : -(id3+1)); // ID of corresponding global track
    if(gid3==gid2 || gid3==gid1){continue;} // Eliminate not-so-evident self-correlations
    // Common track selection criteria for all "normal" global tracks:
    if(!PassesGlobalTrackCuts(gtrack3)){continue;}

    // Okay... So we have three track. Let's decide whether kinemetics from 'atrack' or 'gtrack' is taken, let's check PID, and calculate the background:
    AliAODTrack *agtrack1 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
    AliAODTrack *agtrack2 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
    AliAODTrack *agtrack3 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
    if(fFillControlHistogramsWithGlobalTrackInfo)
    {
     agtrack1 = gtrack1;
     agtrack2 = gtrack2;
     agtrack3 = gtrack3;
    }
    else
    {
     agtrack1 = atrack1;
     agtrack2 = atrack2;
     agtrack3 = atrack3;
    } // if(fFillControlHistogramsWithGlobalTrackInfo)
    if(!agtrack1){Fatal(sMethodName.Data(),"!agtrack1");}
    if(!agtrack2){Fatal(sMethodName.Data(),"!agtrack2");}
    if(!agtrack3){Fatal(sMethodName.Data(),"!agtrack3");}

    if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE)))
    {
     // Q3:
     Double_t dQ3 = Q3(agtrack1,agtrack2,agtrack3); // TBI this is Lorentz invariant
     fDistributionsVsQ3[0]->Fill(dQ3);
    } // if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE)))

   } // for(Int_t iTrack3=0;iTrack3<nTracks;iTrack3++)
  } // for(Int_t iTrack2=0;iTrack2<nTracks;iTrack2++)
 } // for(Int_t iTrack1=0;iTrack1<nTracks;iTrack1++)

} // void AliAnalysisTaskMultiparticleFemtoscopy::HybridApproach1stTerm(AliAODEvent *aAOD)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::HybridApproach2ndTerm(AliAODEvent *aAOD, TClonesArray *ca3, TExMap *em3)
{
 // Calculate 2nd term in the 'hybrid approach'.

 // a) Insanity checks;
 // b) Nested loops to calculate the 3p yield for cross-term;

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::HybridApproach2ndTerm(AliAODEvent *aAOD, TClonesArray *ca1, TExMap *em1)";
 if(!aAOD){Fatal(sMethodName.Data(),"!aAOD");}
 if(0 == fGlobalTracksAOD[0]->GetSize()){Fatal(sMethodName.Data(),"0 == fGlobalTracksAOD[0]->GetSize()");} // this case shall be already treated in UserExec by default, so I can re-use it here
 if(!ca3){Fatal(sMethodName.Data(),"!ca1");}
 if(!em3){Fatal(sMethodName.Data(),"!em1");}

 // b) Nested loops to calculate the 3p yield for cross-term:
 Int_t nTracks = aAOD->GetNumberOfTracks(); // applies to first two loops (signal)
 Int_t nTracks3 = ca3->GetEntries(); // loop over background event
 for(Int_t iTrack1=0;iTrack1<nTracks;iTrack1++)
 {
  AliAODTrack *atrack1 = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack1));
  // TBI Temporary track insanity checks:
  if(!atrack1){Fatal(sMethodName.Data(),"!atrack1");} // TBI keep this for some time, eventually just continue
  if(atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(atrack1->TestFilterBit(128) && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->TestFiletrBit(128) && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(!PassesCommonTrackCuts(atrack1)){continue;} // TBI re-think
  // Corresponding AOD global track:
  Int_t id1 = atrack1->GetID();
  AliAODTrack *gtrack1 = dynamic_cast<AliAODTrack*>(id1>=0 ? aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(id1)) : aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(-(id1+1))));
  if(!gtrack1){Fatal(sMethodName.Data(),"!gtrack1");} // TBI keep this for some time, eventually just continue
  Int_t gid1 = (id1>=0 ? id1 : -(id1+1)); // ID of corresponding global track
  // Common track selection criteria for all "normal" global tracks:
  if(!PassesGlobalTrackCuts(gtrack1)){continue;}

  // Loop over the 2nd particle:
  for(Int_t iTrack2=0;iTrack2<nTracks;iTrack2++) // TBI: I can gain in performance if I do not start from 0
  {
   if(iTrack2<=iTrack1){continue;} // Eliminate self-evident self-correlations, and permutations as well
   //if(iTrack2==iTrack1){continue;} // Eliminate self-evident self-correlations
   AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack2));
   // TBI Temporary track insanity checks:
   if(!atrack2){Fatal(sMethodName.Data(),"!atrack2");} // TBI keep this for some time, eventually just continue
   if(atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(atrack2->TestFilterBit(128) && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->TestFiletrBit(128) && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(!PassesCommonTrackCuts(atrack2)){continue;} // TBI re-think
   // Corresponding AOD global track:
   Int_t id2 = atrack2->GetID();
   AliAODTrack *gtrack2 = dynamic_cast<AliAODTrack*>(id2>=0 ? aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(id2)) : aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(-(id2+1))));
   if(!gtrack2){Fatal(sMethodName.Data(),"!gtrack2");} // TBI keep this for some time, eventually just continue
   Int_t gid2 = (id2>=0 ? id2 : -(id2+1)); // ID of corresponding global track
   if(gid1==gid2){continue;} // Eliminate self-correlations:
   // Common track selection criteria for all "normal" global tracks:
   if(!PassesGlobalTrackCuts(gtrack2)){continue;}

   // Start loop over tracks in the 3rd event (background):
   for(Int_t iTrack3=0;iTrack3<nTracks3;iTrack3++)
   {
    // Remark: Since this is a separate event, I do not need to care about permutations, etc.
    AliAODTrack *atrack3 = dynamic_cast<AliAODTrack*>(ca3->UncheckedAt(iTrack3));
    // TBI Temporary track insanity checks:
    if(!atrack3){Fatal(sMethodName.Data(),"!atrack3");} // TBI keep this for some time, eventually just continue
    if(atrack3->GetID()>=0 && atrack3->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack3->GetID()>=0 && atrack3->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
    if(atrack3->TestFilterBit(128) && atrack3->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack3->TestFiletrBit(128) && atrack3->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
    if(!PassesCommonTrackCuts(atrack3)){continue;} // TBI re-think
    // Corresponding AOD global track:
    Int_t id3 = atrack3->GetID();
    AliAODTrack *gtrack3 = dynamic_cast<AliAODTrack*>(id3>=0 ? ca3->UncheckedAt(em3->GetValue(id3)) : ca3->UncheckedAt(em3->GetValue(-(id3+1))));
    if(!gtrack3){Fatal(sMethodName.Data(),"!gtrack3");} // TBI keep this for some time, eventually just continue
    // Common track selection criteria for all "normal" global tracks:
    if(!PassesGlobalTrackCuts(gtrack3)){continue;}

    // Okay... So we have three track, first two from same event, the 3rd from background. Let's decide whether kinemetics from 'atrack' or 'gtrack' is taken, let's check PID, and calculate the background:
    AliAODTrack *agtrack1 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
    AliAODTrack *agtrack2 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
    AliAODTrack *agtrack3 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
    if(fFillControlHistogramsWithGlobalTrackInfo)
    {
     agtrack1 = gtrack1;
     agtrack2 = gtrack2;
     agtrack3 = gtrack3;
    }
    else
    {
     agtrack1 = atrack1;
     agtrack2 = atrack2;
     agtrack3 = atrack3;
    } // if(fFillControlHistogramsWithGlobalTrackInfo)
    if(!agtrack1){Fatal(sMethodName.Data(),"!agtrack1");}
    if(!agtrack2){Fatal(sMethodName.Data(),"!agtrack2");}
    if(!agtrack3){Fatal(sMethodName.Data(),"!agtrack3");}

    if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE)))
    {
     // Q3:
     Double_t dQ3 = Q3(agtrack1,agtrack2,agtrack3); // TBI this is Lorentz invariant
     fDistributionsVsQ3[1]->Fill(dQ3);
    } // if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE)))

   } // for(Int_t iTrack3=0;iTrack3<nTracks3;iTrack3++)
  } // for(Int_t iTrack2=0;iTrack2<nTracks;iTrack2++)
 } // for(Int_t iTrack1=0;iTrack1<nTracks;iTrack1++)

} // void AliAnalysisTaskMultiparticleFemtoscopy::HybridApproach2ndTerm(AliAODEvent *aAOD, TClonesArray *ca1, TExMap *em1)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::HybridApproach5thTerm(TClonesArray *ca1, TClonesArray *ca2, TClonesArray *ca3, TExMap *em1, TExMap *em2, TExMap *em3)
{
 // a) Insanity checks;
 // b) Three nested loops to calculate B(k);

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::HybridApproach5thTerm(TClonesArray *ca1, TClonesArray *ca2, TClonesArray *ca3, TExMap *em1, TExMap *em2, TExMap *em3)";
 if(!ca1){Fatal(sMethodName.Data(),"!ca1");}
 if(!ca2){Fatal(sMethodName.Data(),"!ca2");}
 if(!ca3){Fatal(sMethodName.Data(),"!ca3");}
 if(!em1){Fatal(sMethodName.Data(),"!em1");}
 if(!em2){Fatal(sMethodName.Data(),"!em2");}
 if(!em3){Fatal(sMethodName.Data(),"!em3");}

 // b) Three nested loops to calculate B(k):
 Int_t nTracks1 = ca1->GetEntries();
 Int_t nTracks2 = ca2->GetEntries();
 Int_t nTracks3 = ca3->GetEntries();

 // Start loop over tracks in the 1st event:
 for(Int_t iTrack1=0;iTrack1<nTracks1;iTrack1++)
 {
  AliAODTrack *atrack1 = dynamic_cast<AliAODTrack*>(ca1->UncheckedAt(iTrack1)); // TBI cross-check UncheckedAt
  // TBI Temporary track insanity checks:
  if(!atrack1){Fatal(sMethodName.Data(),"!atrack1");} // TBI keep this for some time, eventually just continue
  if(atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->GetID()>=0 && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(atrack1->TestFilterBit(128) && atrack1->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack1->TestFiletrBit(128) && atrack1->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
  if(!PassesCommonTrackCuts(atrack1)){continue;} // TBI re-think
  // Corresponding AOD global track:
  Int_t id1 = atrack1->GetID();
  AliAODTrack *gtrack1 = dynamic_cast<AliAODTrack*>(id1>=0 ? ca1->UncheckedAt(em1->GetValue(id1)) : ca1->UncheckedAt(em1->GetValue(-(id1+1))));
  if(!gtrack1){Fatal(sMethodName.Data(),"!gtrack1");} // TBI keep this for some time, eventually just continue
  // Common track selection criteria for all "normal" global tracks:
  if(!PassesGlobalTrackCuts(gtrack1)){continue;}

  // Start loop over tracks in the 2nd event:
  for(Int_t iTrack2=0;iTrack2<nTracks2;iTrack2++)
  {
   AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(ca2->UncheckedAt(iTrack2));
   // TBI Temporary track insanity checks:
   if(!atrack2){Fatal(sMethodName.Data(),"!atrack2");} // TBI keep this for some time, eventually just continue
   if(atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->GetID()>=0 && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(atrack2->TestFilterBit(128) && atrack2->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack2->TestFiletrBit(128) && atrack2->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(!PassesCommonTrackCuts(atrack2)){continue;} // TBI re-think
   // Corresponding AOD global track:
   Int_t id2 = atrack2->GetID();
   AliAODTrack *gtrack2 = dynamic_cast<AliAODTrack*>(id2>=0 ? ca2->UncheckedAt(em2->GetValue(id2)) : ca2->UncheckedAt(em2->GetValue(-(id2+1))));
   if(!gtrack2){Fatal(sMethodName.Data(),"!gtrack2");} // TBI keep this for some time, eventually just continue
   // Common track selection criteria for all "normal" global tracks:
   if(!PassesGlobalTrackCuts(gtrack2)){continue;}

   // Start loop over tracks in the 3rd event:
   for(Int_t iTrack3=0;iTrack3<nTracks3;iTrack3++)
   {
    AliAODTrack *atrack3 = dynamic_cast<AliAODTrack*>(ca3->UncheckedAt(iTrack3));
    // TBI Temporary track insanity checks:
    if(!atrack3){Fatal(sMethodName.Data(),"!atrack3");} // TBI keep this for some time, eventually just continue
    if(atrack3->GetID()>=0 && atrack3->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack3->GetID()>=0 && atrack3->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
    if(atrack3->TestFilterBit(128) && atrack3->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack3->TestFiletrBit(128) && atrack3->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
    if(!PassesCommonTrackCuts(atrack3)){continue;} // TBI re-think
    // Corresponding AOD global track:
    Int_t id3 = atrack3->GetID();
    AliAODTrack *gtrack3 = dynamic_cast<AliAODTrack*>(id3>=0 ? ca3->UncheckedAt(em3->GetValue(id3)) : ca3->UncheckedAt(em3->GetValue(-(id3+1))));
    if(!gtrack3){Fatal(sMethodName.Data(),"!gtrack3");} // TBI keep this for some time, eventually just continue
    // Common track selection criteria for all "normal" global tracks:
    if(!PassesGlobalTrackCuts(gtrack3)){continue;}

    // Okay... So we have three tracks from three different events ready. let's decide whether kinemetics from 'atrack' or 'gtrack' is taken, let's check PID, and calculate the background:
    AliAODTrack *agtrack1 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
    AliAODTrack *agtrack2 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
    AliAODTrack *agtrack3 = NULL; // either 'atrack' or 'gtrack', depending on the flag fFillControlHistogramsWithGlobalTrackInfo. By default, agtrack = 'atrack'
    if(fFillControlHistogramsWithGlobalTrackInfo)
    {
     agtrack1 = gtrack1;
     agtrack2 = gtrack2;
     agtrack3 = gtrack3;
    }
    else
    {
     agtrack1 = atrack1;
     agtrack2 = atrack2;
     agtrack3 = atrack3;
    } // if(fFillControlHistogramsWithGlobalTrackInfo)
    if(!agtrack1){Fatal(sMethodName.Data(),"!agtrack1");}
    if(!agtrack2){Fatal(sMethodName.Data(),"!agtrack2");}
    if(!agtrack3){Fatal(sMethodName.Data(),"!agtrack3");}

    if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE)))
    {
     // Q3:
     Double_t dQ3 = Q3(agtrack1,agtrack2,agtrack3); // TBI this is Lorentz invariant
     fDistributionsVsQ3[4]->Fill(dQ3);
    } // if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE) && Pion(gtrack3,1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE) && Pion(gtrack3,-1,kTRUE)))

   } // for(Int_t iTrack3=0;iTrack3<nTracks3;iTrack3++)
  } // for(Int_t iTrack2=0;iTrack2<nTracks2;iTrack2++)
 } // for(Int_t iTrack1=0;iTrack1<nTracks1;iTrack1++)

} // void AliAnalysisTaskMultiparticleFemtoscopy::HybridApproach5thTerm(TClonesArray *ca1, TClonesArray *ca2, TClonesArray *ca3, TExMap *em1, TExMap *em2, TExMap *em3)

//=======================================================================================================================



