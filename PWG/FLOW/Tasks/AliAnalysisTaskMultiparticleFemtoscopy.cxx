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

#include "TCanvas.h" // TBI
#include "TFile.h" // TBI

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskMultiparticleFemtoscopy)

//================================================================================================================

AliAnalysisTaskMultiparticleFemtoscopy::AliAnalysisTaskMultiparticleFemtoscopy(const char *name, Bool_t useParticleWeights): 
 AliAnalysisTaskSE(name), 
 fHistList(NULL),
 fAnalysisType(NULL),
 fPIDResponse(NULL),
 fMaxNoGlobalTracksAOD(3),
 fProcessBothKineAndReco(kFALSE),
 fRejectFakeTracks(kTRUE),
 fMC(NULL),
 // 1.) Control histograms:
 fControlHistogramsList(NULL),
 fControlHistogramsFlagsPro(NULL),
 fFillControlHistograms(kTRUE), // re-think
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
 fControlHistogramsIdentifiedParticlesList(NULL),
 fControlHistogramsIdentifiedParticlesFlagsPro(NULL),
 fFillControlHistogramsIdentifiedParticles(kFALSE),
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
 fFillCorrelationFunctions(kFALSE),
 fNormalizeCorrelationFunctions(kFALSE),
 fCorrelationFunctionsIndices(NULL),
 // 4.) Background:
 fBackgroundList(NULL),
 fBackgroundFlagsPro(NULL),
 fEstimateBackground(kFALSE),
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
  fHistList->SetName("GMlist"); // TBI
  fHistList->SetOwner(kTRUE);

  // Initialize all arrays:
  this->InitializeArrays();
  this->InitializeArraysForControlHistograms();
  this->InitializeArraysForEBEObjects();
  this->InitializeArraysForCorrelationFunctions();
  this->InitializeArraysForBackground();

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
 fMaxNoGlobalTracksAOD(3),
 fProcessBothKineAndReco(kFALSE),
 fRejectFakeTracks(kTRUE),
 fMC(NULL),
 // 1.) Control histograms:
 fControlHistogramsList(NULL),
 fControlHistogramsFlagsPro(NULL),
 fFillControlHistograms(kTRUE), // TBI re-think
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
 fControlHistogramsIdentifiedParticlesList(NULL),
 fControlHistogramsIdentifiedParticlesFlagsPro(NULL),
 fFillControlHistogramsIdentifiedParticles(kFALSE),
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
 fFillCorrelationFunctions(kFALSE),
 fNormalizeCorrelationFunctions(kFALSE),
 fCorrelationFunctionsIndices(NULL),
 // 4.) Background:
 fBackgroundList(NULL),
 fBackgroundFlagsPro(NULL),
 fEstimateBackground(kFALSE),
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
 // c) Book and nest all lists;
 // d) Book all objects;
 // e) Set all flags;
 // f) Trick to avoid name clashes, part 2.

 // a) Insanity checks:
 this->InsanityChecksUserCreateOutputObjects();

 // b) Trick to avoid name clashes, part 1:
 Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
 TH1::AddDirectory(kFALSE);

 // c) Book and nest all lists:
 this->BookAndNestAllLists();

 // d) Book all objects:
 this->BookEverything();
 this->BookEverythingForControlHistograms();
 this->BookEverythingForEBEObjects();
 this->BookEverythingForCorrelationFunctions();
 this->BookEverythingForBackground();

 // e) Set all flags:
 // ...

 // f) Trick to avoid name clashes, part 2:
 TH1::AddDirectory(oldHistAddStatus);

 PostData(1,fHistList);

} // void AliAnalysisTaskMultiparticleFemtoscopy::UserCreateOutputObjects() 

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::UserExec(Option_t *) 
{
 // Main loop (called for each event).

 // a) Determine Ali{MC,ESD,AOD}Event;
 // b) Insanity checks;
 // c) Start analysis over MC, ESD or AODs;
 // d) Reset event-by-event objects;
 // e) PostData;
 // f) Online monitoring.

 // a) Determine Ali{MC,ESD,AOD}Event:
 AliMCEvent *aMC = MCEvent();                                  // from TaskSE
 AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(InputEvent()); // from TaskSE
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(InputEvent()); // from TaskSE

 // b) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::UserExec(Option_t *)";
 if(aMC && aAOD && fProcessBothKineAndReco){*fAnalysisType="MC_AOD";}
 else if(aMC){*fAnalysisType="MC";}
 else if(aESD){*fAnalysisType="ESD";}
 else if(aAOD){*fAnalysisType="AOD";}
 else{Fatal(sMethodName.Data(),"Neither MC, nor ESD, nor AOD, nor MC_AOD. Okay...");}
 InsanityChecksUserExec();

 // c) Start analysis over MC, ESD or AOD:
 if(fProcessBothKineAndReco) // TBI: this is a bit shaky, isn't it...
 {
  fMC = aMC;
  AOD(aAOD);
 }
 else if(aMC){MC(aMC);}
 else if(aESD){ESD(aESD);}
 else if(aAOD){AOD(aAOD);}
 else{Fatal(sMethodName.Data(),"Neither MC, nor ESD, nor AOD. Okay...");} // TBI I have this check already

 // d) Reset event-by-event objects:
 this->ResetEBEObjects();

 // e) PostData:
 PostData(1,fHistList);

 // f) Online monitoring:
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
 fHistList = (TList*)GetOutputData(1);
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

void AliAnalysisTaskMultiparticleFemtoscopy::NormalizeCorrelationFunctions()
{
 // Normalize correlation functions with the background.

 for(Int_t pid1=0;pid1<10;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pid2=0;pid2<10;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   if(fCorrelationFunctions[pid1][pid2] && fCorrelationFunctions[pid1][pid2]->GetEntries() > 0 &&
      fBackground[pid1][pid2] && fBackground[pid1][pid2]->GetEntries() > 0)
   {
    fCorrelationFunctions[pid1][pid2]->Divide(fBackground[pid1][pid2]);
   } // if(..)
  } // for(Int_t pid2=0;pid2<10;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 } // for(Int_t pid1=0;pid1<10;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]

} // void AliAnalysisTaskMultiparticleFemtoscopy::NormalizeCorrelationFunctions()

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsEvent(AliVEvent *ave)
{
 // Fill control histograms for global event observables.

 // TBI: add support for fProcessBothKineAndReco

 // To do:
 // 1) Add support for MC and ESD.

 // a) Determine Ali{MC,ESD,AOD}Event;
 // b) ...

 // a) Determine Ali{MC,ESD,AOD}Event:
 AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 if(aMC)
 {
  // MC event:
  fGetNumberOfTracksHist->Fill(aMC->GetNumberOfTracks());
 } // if(aMC)
 else if(aESD)
 {
  // TBI
 } // else if(aESD)
 else if(aAOD)
 {
  // AOD event:
  fGetNumberOfTracksHist->Fill(aAOD->GetNumberOfTracks()); // TBI not all tracks are unique
  fGetNumberOfGlobalTracksHist->Fill(fGlobalTracksAOD[0]->GetSize()); // this is then my multiplicity...
  fGetNumberOfV0sHist->Fill(aAOD->GetNumberOfV0s()); // TBI some V0s share the daughter
  fGetNumberOfCascadesHist->Fill(aAOD->GetNumberOfCascades()); // TBI validate
  fGetMagneticFieldHist->Fill(aAOD->GetMagneticField());
  fGetEventTypeHist->Fill(aAOD->GetEventType());
  if(aAOD->GetCentrality())
  {
   fGetCentralityHist->Fill(aAOD->GetCentrality()->GetCentralityPercentile("V0M")); // TBI validate this
  }
  // AOD primary vertex:
  AliAODVertex *avtx = (AliAODVertex*)aAOD->GetPrimaryVertex();
  fVertexXYZ[0]->Fill(avtx->GetX());
  fVertexXYZ[1]->Fill(avtx->GetY());
  fVertexXYZ[2]->Fill(avtx->GetZ());
  fGetNContributorsHist->Fill(avtx->GetNContributors());
  fGetChi2perNDFHist->Fill(avtx->GetChi2perNDF());
  fGetNDaughtersHist->Fill(avtx->GetNDaughters());
 } // else if(aAOD)
 
} // void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsEvent(AliVEvent *ave)

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsParticle(AliVEvent *ave)
{
 // Fill control histograms for particles.
 // Remark: The idea is to have one loop over the particles and to fill everything which needs to be filled.

 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsParticle(AliVEvent *ave)";

 // a) Determine Ali{MC,ESD,AOD}Event:
 AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 if(aMC)
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

 else if(aAOD)
 {
  // AOD tracks:
  for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++)
  {
   // a) Determine "atrack";
   // b) Insanity checks for "atrack";
   // c) Determine the corresponding "gtrack" (a.k.a. "normal global" track);
   // d) Fill the control histograms for non-identified particles;
   // e) Fill the control histograms for identified particles;

   // a) Determine "atrack":
   AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack));

   // b) Insanity checks for "atrack":
   if(!atrack){Fatal(sMethodName.Data(),"!atrack");}
   if(atrack->GetID()>=0 && atrack->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack->GetID()>=0 && atrack->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   if(atrack->TestFilterBit(128) && atrack->IsGlobalConstrained()){Fatal(sMethodName.Data(),"atrack->TestFiletrBit(128) && atrack->IsGlobalConstrained()");} // TBI keep this for some time, eventually just continue
   // if(0==atrack->GetFilterMap())
   // TBI a vast majority of such tracks pass
   //   if(0==atrack->GetFilterMap() && atrack->GetType() == AliAODTrack::kFromDecayVtx)
   // but there are exceptions which pass
   //   if(0==atrack->GetFilterMap() && atrack->GetType() == AliAODTrack::kPrimary)
   // So it's a mess...
   // 1) Clearly, for none of such tracks we can test a filter bit
   // 2) Apparently all of them have positive ID, meaning they are global (?)
    // Corresponding AOD global track:

   // c) Determine the corresponding "gtrack" (a.k.a. "normal global" track):
   if(0 == fGlobalTracksAOD[0]->GetSize()){GlobalTracksAOD(aAOD,0);}
   if(0 == fGlobalTracksAOD[0]->GetSize()){return;}
   Int_t id = atrack->GetID();
   AliAODTrack *gtrack = dynamic_cast<AliAODTrack*>(id>=0 ? aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(id)) : aAOD->GetTrack(fGlobalTracksAOD[0]->GetValue(-(id+1))));
   if(!gtrack){Fatal(sMethodName.Data(),"!gtrack");}

   // d) Fill the control histograms for non-identified particles:
   if(fFillControlHistogramsNonIdentifiedParticles){this->FillControlHistogramsNonIdentifiedParticles(gtrack);}

   // e) Fill the control histograms for identified particles:
   if(fFillControlHistogramsIdentifiedParticles){this->FillControlHistogramsIdentifiedParticles(atrack,gtrack);}

  } // for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++)

 } // else if(aAOD)

} // void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsParticle(AliVEvent *ave)

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsNonIdentifiedParticles(AliAODTrack *gtrack)
{
 // Fill control histograms for non-identified particles.
 // Uses PassesCommonGlobalTrackCuts(...) TBI

 // a) Insanity checks;
 // b) Check cut selection criteria;
 // c) Fill control histograms.

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsNonIdentifiedParticles(AliAODTrack *gtrack)";
 if(!gtrack){Fatal(sMethodName.Data(),"!gtrack");} // TBI keep this for some time, eventually just continue

 // b) Check cut selection criteria:
 if(!PassesCommonGlobalTrackCuts(gtrack)){return;} // TBI in the method PassesCommonGlobalTrackCuts track is hardwired to AliAODTrack. Try to generalize

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

 // a) Insanity checks;
 // b) Check cut selection criteria;
 // c) Fill control histograms.

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsNonIdentifiedParticles(AliAODMCParticle *amcparticle)";
 if(!amcparticle){Fatal(sMethodName.Data(),"!amcparticle");}

 // b) Check cut selection criteria:
 //if(!PassesCommonGlobalTrackCuts(gtrack)){return;} // TBI in the method PassesCommonGlobalTrackCuts track is hardwired to AliAODTrack. Try to generalize

 // c) Fill control histograms:
 fChargeHist->Fill(amcparticle->Charge()+0.5); // see how this histogram was booked for convention used
 fPtHist->Fill(amcparticle->Pt());
 fEtaHist->Fill(amcparticle->Eta());
 fPhiHist->Fill(amcparticle->Phi());
 fMassHist->Fill(amcparticle->M());
 fGetPdgCode->Fill(amcparticle->GetPdgCode());

} // void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsNonIdentifiedParticles(AliAODMCParticle *amcparticle)

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsIdentifiedParticles(AliAODTrack *atrack, AliAODTrack *gtrack)
{
 // Fill control histograms for identified particles.
 // Uses PassesCommonTrackCuts(...) TBI yet still it stored the paremeters of corresponding global track... TBI

 // To do:
 // 1) Add support for MC and ESD, now it works only for AliAODTrack.

 // a) Insanity checks;
 // b) Check cut selection criteria;
 // c) Fill control histograms.

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::FillControlHistogramsIdentifiedParticles(AliAODTrack *atrack, AliAODTrack *gtrack)";
 if(!atrack){Fatal(sMethodName.Data(),"!atrack");} // TBI keep this for some time, eventually just continue
 if(!gtrack){Fatal(sMethodName.Data(),"!gtrack");} // TBI keep this for some time, eventually just continue

 // b) Check cut selection criteria:
 if(!PassesCommonTrackCuts(atrack)){return;} // TBI in the method PassesCommonGlobalTrackCuts track is hardwired to AliAODTrack. Try to generalize
 // TBI do I need some check for atrack?

 // c) Fill control histograms:
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
    fPtPIDHist[2][charge][ps]->Fill(gtrack->Pt());
    fMassPIDHist[2][charge][ps]->Fill(gtrack->M());
    fEtaPIDHist[2][charge][ps]->Fill(gtrack->Eta());
    fPhiPIDHist[2][charge][ps]->Fill(gtrack->Phi());
   }
   // c3) Kaons:
   if(Kaon(gtrack,nCharge,bPrimary))
   {
    fPtPIDHist[3][charge][ps]->Fill(gtrack->Pt());
    fMassPIDHist[3][charge][ps]->Fill(gtrack->M());
    fEtaPIDHist[3][charge][ps]->Fill(gtrack->Eta());
    fPhiPIDHist[3][charge][ps]->Fill(gtrack->Phi());
   }
   // c4) Protons:
   if(Proton(gtrack,nCharge,bPrimary))
   {
    fPtPIDHist[4][charge][ps]->Fill(gtrack->Pt());
    fMassPIDHist[4][charge][ps]->Fill(gtrack->M());
    fEtaPIDHist[4][charge][ps]->Fill(gtrack->Eta());
    fPhiPIDHist[4][charge][ps]->Fill(gtrack->Phi());
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
  if(!PassesMixedEventCuts(aMC)){return;} // TBI this is empty at the moment for MC
  Int_t nTracks = aMC->GetNumberOfTracks();
  if(!fMixedEvents[0])
  {
   TClonesArray ca0("AliAODMCParticle"); // temporary holder, its entries will be copied to fMixedEvents[0]
   Int_t ca0Counter = 0;
   for(Int_t iTrack=0;iTrack<nTracks;iTrack++)
   {
    AliAODMCParticle *amcparticle = dynamic_cast<AliAODMCParticle*>(aMC->GetTrack(iTrack));
    if(!amcparticle){Fatal(sMethodName.Data(),"!amcparticle");}
    if(!PassesCommonTrackCuts(amcparticle)){continue;} // TBI re-think, see implemntation of this method
    if(!(pattern.Contains(Form(".%d.",TMath::Abs(amcparticle->GetPdgCode()))))){continue;}
    // Fill fMixedEvents[0] with tracks which survived all checks:
    ca0[ca0Counter++] = amcparticle;
   } // for(Int_t iTrack=0;iTrack<nTracks;iTrack++)
   fMixedEvents[0] = (TClonesArray*)ca0.Clone();
   if(!fMixedEvents[0]){Fatal(sMethodName.Data(),"!fMixedEvents[0]");}
  }
  else if(!fMixedEvents[1])
  {
   TClonesArray ca1("AliAODMCParticle"); // temporary holder, its entries will be copied to fMixedEvents[1]
   Int_t ca1Counter = 0;
   for(Int_t iTrack=0;iTrack<nTracks;iTrack++)
   {
    AliAODMCParticle *amcparticle = dynamic_cast<AliAODMCParticle*>(aMC->GetTrack(iTrack));
    if(!amcparticle){Fatal(sMethodName.Data(),"!amcparticle");}
    if(!PassesCommonTrackCuts(amcparticle)){continue;} // TBI re-think, see implemntation of this method
    if(!(pattern.Contains(Form(".%d.",TMath::Abs(amcparticle->GetPdgCode()))))){continue;}
    // Fill fMixedEvents[0] with tracks which survived all checks:
    ca1[ca1Counter++] = amcparticle;
   } // for(Int_t iTrack=0;iTrack<nTracks;iTrack++)
   fMixedEvents[1] = (TClonesArray*)ca1.Clone();
   if(!fMixedEvents[1]){Fatal(sMethodName.Data(),"!fMixedEvents[1]");}
  }
  // Shall we do something?
  if(fMixedEvents[0] && fMixedEvents[1])
  {
   CalculateBackground(fMixedEvents[0],fMixedEvents[1],kTRUE);
  }
 } // if(aMC)
 else if(aESD)
 {
  // TBI
 } // else if(aESD)
 else if(aAOD)
 {
  if(!PassesMixedEventCuts(aAOD)){return;}
  if(!fMixedEvents[0])
  {
   if(fGlobalTracksAOD[1]) fGlobalTracksAOD[1]->Delete();
   if(0 != fGlobalTracksAOD[1]->GetSize()){Fatal(sMethodName.Data(),"0 != fGlobalTracksAOD[1]->GetSize()");}
   GlobalTracksAOD(aAOD,1);
   if(0 == fGlobalTracksAOD[1]->GetSize()){return;}
   fMixedEvents[0] = (TClonesArray*)aAOD->GetTracks()->Clone();
   if(!fMixedEvents[0]){Fatal(sMethodName.Data(),"!fMixedEvents[0]");}
  }
  else if(!fMixedEvents[1])
  {
   if(fGlobalTracksAOD[2]) fGlobalTracksAOD[2]->Delete();
   if(0 != fGlobalTracksAOD[2]->GetSize()){Fatal(sMethodName.Data(),"0 != fGlobalTracksAOD[2]->GetSize()");}
   GlobalTracksAOD(aAOD,2);
   if(0 == fGlobalTracksAOD[2]->GetSize()){return;}
   fMixedEvents[1] = (TClonesArray*)aAOD->GetTracks()->Clone();
   if(!fMixedEvents[1]){Fatal(sMethodName.Data(),"!fMixedEvents[1]");}
  }
  // Shall we do something?
  if(fMixedEvents[0] && fMixedEvents[1])
  {
   CalculateBackground(fMixedEvents[0],fMixedEvents[1]);
  }
 } // else if(aAOD)

} // void AliAnalysisTaskMultiparticleFemtoscopy::EstimateBackground(AliVEvent *ave)

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

 // TBI

 //

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
  if(fProcessBothKineAndReco && fEstimateBackground){Fatal(sMethodName.Data(),"fProcessBothKineAndReco && fEstimateBackground");}
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
 if(fEstimateBackground){this->EstimateBackground(aMC);}

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
 // c) Filter our "normal global" tracks;
 // d) Fill control histogram for global event observables;
 // e) Fill control histogram for particles;
 // f) Calculate correlation functions;
 // g) Calculate background;
 // h) V0s.

 // a) Debugging:
 if(fDoSomeDebugging){this->DoSomeDebugging(aAOD);}

 // b) Common event selection criteria:
 if(!this->PassesCommonEventCuts(aAOD)){return;}

 // c) Filter out "normal global" tracks:
 this->GlobalTracksAOD(aAOD,0); // [0] stands for default analysis
 //if(0 == fGlobalTracksAOD[0]->GetSize()) return;

 // d) Fill control histogram for global event observables:
 if(fFillControlHistogramsEvent){this->FillControlHistogramsEvent(aAOD);}

 // e) Fill control histograms for particles:
 if(fFillControlHistogramsNonIdentifiedParticles || fFillControlHistogramsIdentifiedParticles){this->FillControlHistogramsParticle(aAOD);}

 // f) Calculate correlation functions:
 if(fFillCorrelationFunctions){this->CalculateCorrelationFunctions(aAOD);}

 // g) Calculate background:
 if(fEstimateBackground){this->EstimateBackground(aAOD);}


 return;


 // h) V0s:
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

 for(Int_t pid1=0;pid1<10;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pid2=0;pid2<10;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   fCorrelationFunctions[pid1][pid2] = NULL;
  }
 }

} // void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForCorrelationFunctions()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForBackground()
{
 // Initialize all arrays for background.

 for(Int_t pid1=0;pid1<10;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pid2=0;pid2<10;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   fBackground[pid1][pid2] = NULL;
  }
 }

 for(Int_t me=0;me<2;me++) // [0] is buffer for 1st event; [1] for 2nd
 {
  fMixedEvents[me] = NULL;
 }

} // void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForBackground()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::BookAndNestAllLists()
{
 // Book and nest all lists nested in the base list fHistList.

 // a) Book and nest lists for control histograms;
 // b) Book and nest lists for eveny-by-event histograms;
 // c) Book and nest lists for correlation functions;
 // d) Book and nest lists for background.

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
 fCorrelationFunctionsList->SetName("Correlation_Functions");
 fCorrelationFunctionsList->SetOwner(kTRUE);
 fHistList->Add(fCorrelationFunctionsList);

 // d) Book and nest lists for background:
 fBackgroundList = new TList();
 fBackgroundList->SetName("Background");
 fBackgroundList->SetOwner(kTRUE);
 fHistList->Add(fBackgroundList);

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

 // Inclusive sigma cuts:
 fInclusiveSigmaCuts[2] = 2.; // i.e. in function Pion(...) the inclusive cut for pions is 2 sigma
 fInclusiveSigmaCuts[3] = 2.; // i.e. in function Kaon(...) the inclusive cut for kaons is 2 sigma
 fInclusiveSigmaCuts[4] = 2.; // i.e. in function Proton(...) the inclusive cut for protons is 2 sigma

 // Exclusive sigma cuts:
 // Pion(...)
 fExclusiveSigmaCuts[2][3] = 4.; // i.e. in function Pion(...) the exclusive cut for kaons is 4 sigma
 fExclusiveSigmaCuts[2][4] = 4.; // i.e. in function Pion(...) the exclusive cut for protons is 4 sigma
 // Kaon(...)
 fExclusiveSigmaCuts[3][2] = 4.; // i.e. in function Kaon(...) the exclusive cut for pions is 4 sigma
 fExclusiveSigmaCuts[3][4] = 4.; // i.e. in function Kaon(...) the exclusive cut for protons is 4 sigma
 // Proton(...)
 fExclusiveSigmaCuts[4][2] = 4.; // i.e. in function Proton(...) the exclusive cut for pions is 4 sigma
 fExclusiveSigmaCuts[4][3] = 4.; // i.e. in function Proton(...) the exclusive cut for kaons is 4 sigma

} // void AliAnalysisTaskMultiparticleFemtoscopy::BookEverything()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForControlHistograms()
{
 // Book all the stuff for control histograms.

 // a) Book the profile holding all the flags for control histograms;
 // b) Common vaiables;
 // c) Book all control histograms...
 //  c0) Event;
 //  c1) Non-identified particles;
 //  c2) Identified particles;
 //  c3) V0s.

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

 if(!fFillControlHistograms){return;} // TBI is this safe?

 // b) Common vaiables:
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
  fGetMagneticFieldHist = new TH1D("fGetMagneticFieldHist","aAOD->GetMagneticField()",20,-10.,10.);
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

 //  c2) Identified particles:
 // Book the profile holding all the flags for TBI:
 fControlHistogramsIdentifiedParticlesFlagsPro = new TProfile("fControlHistogramsIdentifiedParticlesFlagsPro","Flags and settings for TBI",1,0,1);
 fControlHistogramsIdentifiedParticlesFlagsPro->SetTickLength(-0.01,"Y");
 fControlHistogramsIdentifiedParticlesFlagsPro->SetMarkerStyle(25);
 fControlHistogramsIdentifiedParticlesFlagsPro->SetLabelSize(0.04);
 fControlHistogramsIdentifiedParticlesFlagsPro->SetLabelOffset(0.02,"Y");
 fControlHistogramsIdentifiedParticlesFlagsPro->SetStats(kFALSE);
 fControlHistogramsIdentifiedParticlesFlagsPro->SetFillColor(kGray);
 fControlHistogramsIdentifiedParticlesFlagsPro->SetLineColor(kBlack);
 fControlHistogramsIdentifiedParticlesFlagsPro->GetXaxis()->SetBinLabel(1,"fFillControlHistogramsIdentifiedParticles"); fControlHistogramsIdentifiedParticlesFlagsPro->Fill(0.5,fFillControlHistogramsIdentifiedParticles);
 fControlHistogramsIdentifiedParticlesList->Add(fControlHistogramsIdentifiedParticlesFlagsPro);
 if(fFillControlHistogramsIdentifiedParticles)
 {
  for(Int_t pid=0;pid<5;pid++) // [0=e,1=mu,2=pi,3=K,4=p]
  {
   for(Int_t pa=0;pa<2;pa++) // particle(+q)/antiparticle(-q)
   {
    for(Int_t ps=0;ps<2;ps++) // kPrimary/kFromDecayVtx
    {
     fMassPIDHist[pid][pa][ps] = new TH1F(Form("fMassPIDHist[%d][%d][%d]",pid,pa,ps),Form("fMassPIDHist[%d][%d][%d] (%s)",pid,pa,ps,sParticleLabel[pid].Data()),10000,0.,10.);
     fMassPIDHist[pid][pa][ps]->GetXaxis()->SetTitle("m [GeV/c^{2}]");
     for(Int_t nm=0;nm<5;nm++) // nominal masses
     {
      fMassPIDHist[pid][pa][ps]->GetXaxis()->SetBinLabel(fMassPIDHist[pid][pa][ps]->FindBin(dNominalMass[nm]),Form("m_{%s}",sParticleLabel[nm].Data()));
     }
     fMassPIDHist[pid][pa][ps]->SetFillColor(kBlue-10);
     fControlHistogramsIdentifiedParticlesList->Add(fMassPIDHist[pid][pa][ps]);
     fPtPIDHist[pid][pa][ps] = new TH1F(Form("fPtPIDHist[%d][%d][%d]",pid,pa,ps),Form("fPtPIDHist[%d][%d][%d]",pid,pa,ps),1000,0.,10.);
     fPtPIDHist[pid][pa][ps]->GetXaxis()->SetTitle("p_{T} [TBI units]");
     fPtPIDHist[pid][pa][ps]->SetFillColor(kBlue-10);
     fControlHistogramsIdentifiedParticlesList->Add(fPtPIDHist[pid][pa][ps]);
     fEtaPIDHist[pid][pa][ps] = new TH1F(Form("fEtaPIDHist[%d][%d][%d]",pid,pa,ps),Form("fEtaPIDHist[%d][%d][%d]",pid,pa,ps),200000,-2.,2.);
     fEtaPIDHist[pid][pa][ps]->SetFillColor(kBlue-10);
     fControlHistogramsIdentifiedParticlesList->Add(fEtaPIDHist[pid][pa][ps]);
     fPhiPIDHist[pid][pa][ps] = new TH1F(Form("fPhiPIDHist[%d][%d][%d]",pid,pa,ps),Form("fPhiPIDHist[%d][%d][%d]",pid,pa,ps),360,0.,TMath::TwoPi());
     fPhiPIDHist[pid][pa][ps]->GetXaxis()->SetTitle("#phi");
     fPhiPIDHist[pid][pa][ps]->SetFillColor(kBlue-10);
     fControlHistogramsIdentifiedParticlesList->Add(fPhiPIDHist[pid][pa][ps]);
    } // for(Int_t ps=0;ps<2;ps++) // kPrimary/kFromDecayVtx
   } // for(Int_t pa=0;pa<2;pa++) // particle(+q)/antiparticle(-q)
  } // for(Int_t pid=0;pid<4;pid++) // [0=e,1=mu,2=pi,3=K,4=p]
 } // if(fFillControlHistogramsIdentifiedParticles)

 //  c3) V0s:
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
 // b) Book all correlation functions;
 // c) Book TExMap *fCorrelationFunctionsIndices.

 // a) Book the profile holding all the flags for correlation functions objects:
 fCorrelationFunctionsFlagsPro = new TProfile("fCorrelationFunctionsFlagsPro","Flags and settings for correlation functions histograms",2,0,2);
 fCorrelationFunctionsFlagsPro->SetTickLength(-0.01,"Y");
 fCorrelationFunctionsFlagsPro->SetMarkerStyle(25);
 fCorrelationFunctionsFlagsPro->SetLabelSize(0.04);
 fCorrelationFunctionsFlagsPro->SetLabelOffset(0.02,"Y");
 fCorrelationFunctionsFlagsPro->SetStats(kFALSE);
 fCorrelationFunctionsFlagsPro->SetFillColor(kGray);
 fCorrelationFunctionsFlagsPro->SetLineColor(kBlack);
 fCorrelationFunctionsFlagsPro->GetXaxis()->SetBinLabel(1,"fFillCorrelationFunctions"); fCorrelationFunctionsFlagsPro->Fill(0.5,fFillCorrelationFunctions);
 fCorrelationFunctionsFlagsPro->GetXaxis()->SetBinLabel(2,"fNormalizeCorrelationFunctions"); fCorrelationFunctionsFlagsPro->Fill(1.5,fNormalizeCorrelationFunctions);
 fCorrelationFunctionsList->Add(fCorrelationFunctionsFlagsPro);

 if(!fFillCorrelationFunctions){return;} // TBI is this safe?

 const Int_t nParticleSpecies = 5; // Supported at the moment: 0=e,1=mu,2=pi,3=K,4=p
 TString sParticles[2*nParticleSpecies] = {"e^{+}","#mu^{+}","#pi^{+}","K^{+}","p^{+}","e^{-}","#mu^{-}","#pi^{-}","K^{-}","p^{-}"};

 // b) Book all correlation functions:
 // Remark 0: First particle in the pair is always the one with positive charge.
 // Remark 1: Diagonal elements are particle/antiparticle pairs.
 for(Int_t pid1=0;pid1<2*nParticleSpecies;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pid2=pid1;pid2<2*nParticleSpecies;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   // Correlation functions:
   fCorrelationFunctions[pid1][pid2] = new TH1F(Form("fCorrelationFunctions[%d][%d]",pid1,pid2),Form("fCorrelationFunctions[%d][%d] = (%s,%s)",pid1,pid2,sParticles[pid1].Data(),sParticles[pid2].Data()),10000,0.,10.); // TBI rethink the boundaries and nbins
   fCorrelationFunctions[pid1][pid2]->SetStats(kFALSE);
   fCorrelationFunctions[pid1][pid2]->SetFillColor(kBlue-10);
   fCorrelationFunctions[pid1][pid2]->SetXTitle("k");
   fCorrelationFunctions[pid1][pid2]->SetYTitle("C(k)");
   fCorrelationFunctionsList->Add(fCorrelationFunctions[pid1][pid2]);
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

} // void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForCorrelationFunctions()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForBackground()
{
 // Book all the stuff for background objects.

 // a) Book the profile holding all the flags for background objects;
 // b) Book all histograms for background;
 // c) Book buffer objects.

 // a) Book the profile holding all the flags for correlation functions objects:
 fBackgroundFlagsPro = new TProfile("fBackgroundFlagsPro","Flags and settings for background histograms",1,0,1);
 fBackgroundFlagsPro->SetTickLength(-0.01,"Y");
 fBackgroundFlagsPro->SetMarkerStyle(25);
 fBackgroundFlagsPro->SetLabelSize(0.04);
 fBackgroundFlagsPro->SetLabelOffset(0.02,"Y");
 fBackgroundFlagsPro->SetStats(kFALSE);
 fBackgroundFlagsPro->SetFillColor(kGray);
 fBackgroundFlagsPro->SetLineColor(kBlack);
 fBackgroundFlagsPro->GetXaxis()->SetBinLabel(1,"fEstimateBackground"); fBackgroundFlagsPro->Fill(0.5,fEstimateBackground);
 fBackgroundList->Add(fBackgroundFlagsPro);

 if(!fEstimateBackground){return;} // TBI is this safe?

 const Int_t nParticleSpecies = 5; // Supported at the moment: 0=e,1=mu,2=pi,3=K,4=p
 TString sParticles[2*nParticleSpecies] = {"e^{+}","#mu^{+}","#pi^{+}","K^{+}","p^{+}","e^{-}","#mu^{-}","#pi^{-}","K^{-}","p^{-}"};

 // b) Book all histograms for background:
 // Remark 0: First particle in the pair is always the one with positive charge.
 // Remark 1: Diagonal elements are particle/antiparticle pairs.
 for(Int_t pid1=0;pid1<2*nParticleSpecies;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pid2=pid1;pid2<2*nParticleSpecies;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   // Background:
   fBackground[pid1][pid2] = new TH1F(Form("fBackground[%d][%d]",pid1,pid2),Form("fBackground[%d][%d] = (%s,%s)",pid1,pid2,sParticles[pid1].Data(),sParticles[pid2].Data()),10000,0.,10.); // TBI rethink the boundaries and nbins
   fBackground[pid1][pid2]->SetStats(kFALSE);
   fBackground[pid1][pid2]->SetFillColor(kBlue-10);
   fBackground[pid1][pid2]->SetXTitle("k");
   fBackground[pid1][pid2]->SetYTitle("B(k)");
   fBackgroundList->Add(fBackground[pid1][pid2]);
  } // for(Int_t pid2=0;pid2<5;pid2++) // anti-particle [0=e,1=mu,2=pi,3=K,4=p]
 } // for(Int_t pid=0;pid<5;pid++) // particle [0=e,1=mu,2=pi,3=K,4=p]

 // c) Book buffer objects:
 // TBI not sure why I need this here again. Re-think...
 for(Int_t me=0;me<2;me++) // [0] is buffer for 1st event; [1] for 2nd
 {
  fMixedEvents[me] = NULL;
 }

} // void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForBackground()

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

 for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++)
 {
  AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack));
  if(aodTrack)
  {
   Int_t id = aodTrack->GetID();
   //if(id>=0 && aodTrack->GetFilterMap()>0 && !aodTrack->IsGlobalConstrained()) // TBI rethink this
   if(id>=0 && !aodTrack->IsGlobalConstrained()) // TBI rethink this, it seems that id>=0 is just enough, the second constraint is most likely just an overkill
   {
    fGlobalTracksAOD[index]->Add(id,iTrack);
   } // if(id>=0 && !aodTrack->IsGlobalConstrained())
  } // if(aodTrack)
 } // for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++)

} // void AliAnalysisTaskMultiparticleFemtoscopy::GlobalTracksAOD(AliAODEvent *aAOD)

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
Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesCommonGlobalTrackCuts(AliAODTrack *gtrack)
{
 // Check if the track passes common global track cuts (irrespectively of PID).

 TString sMethodName = "Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesCommonGlobalTrackCuts(AliAODTrack *gtrack)";
 if(!gtrack){Fatal(sMethodName.Data(),"!gtrack");}

 // To do: add data members and corresponding setters:
 // fPtMin, fPtMax
 // fEtaMin, fEtaMax
 // fPhiMin, fPhiMax
 // fTPCNclsMin, fTPCNclsMax
 // fTPCsignalNMin, fTPCsignalNMax

 if(gtrack->Pt()<0.2) return kFALSE;
 if(gtrack->Pt()>=5.0) return kFALSE;
 if(gtrack->Eta()<-0.8) return kFALSE;
 if(gtrack->Eta()>=0.8) return kFALSE;
 //if(gtrack->Phi()<-0.6) return kFALSE;
 //if(gtrack->Phi()>=0.6) return kFALSE;
 //if(gtrack->GetTPCNcls()<70) return kFALSE;
 //if(gtrack->GetTPCsignalN()<70) return kFALSE;

 if(fRejectFakeTracks && gtrack->GetLabel()<0) return kFALSE; // TBI is it more precise <=0 ?

 return kTRUE;

} // Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesCommonGlobalTrackCuts(AliAODTrack *gtrack)

//=======================================================================================================================

Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesCommonTrackCuts(AliAODTrack *atrack)
{
 // Check if the track passes common analysis specific track (e.g. TPC-only) cuts.
 // Thereforem we can cut independetly on global track parameters, and on TPC-only cut parameters.

 TString sMethodName = "Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesCommonTrackCuts(AliAODTrack *atrack)";
 if(!atrack){Fatal(sMethodName.Data(),"!atrack");}

 if(!atrack->TestFilterBit(128)) return kFALSE; // TPC-only TBI setter

 if(fRejectFakeTracks && atrack->GetLabel()<0) return kFALSE; // TBI is it more precise <=0 ?

 return kTRUE; 

} // Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesCommonTrackCuts(AliAODTrack *atrack)

//=======================================================================================================================

Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesCommonTrackCuts(AliAODMCParticle *amcparticle)
{
 // TBI this method applies only to MC, make it uniform with AOD.

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
  if(!aAOD->GetPrimaryVertex()) return kFALSE;
  if(TMath::Abs(aAOD->GetMagneticField())<0.001) return kFALSE;
  AliAODVertex *avtx = (AliAODVertex*)aAOD->GetPrimaryVertex();
  //if(TMath::Abs(avtx->GetX())>10.0) return kFALSE;
  //if(TMath::Abs(avtx->GetY())>10.0) return kFALSE;
  if(TMath::Abs(avtx->GetZ())>10.0) return kFALSE; // TBI setter
  if(avtx->GetNContributors()<=2) return kFALSE; // TBI setter
 } // else if(aAOD)

 return kTRUE;

} // Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesCommonEventCuts(AliAODEvent *aAOD)

//=======================================================================================================================

Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesMixedEventCuts(AliVEvent *ave)
{
 // Check if the event passes mixed event cuts.

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
  if(!aAOD->GetPrimaryVertex()) return kFALSE;
  if(TMath::Abs(aAOD->GetMagneticField())<0.001) return kFALSE;
  AliAODVertex *avtx = (AliAODVertex*)aAOD->GetPrimaryVertex();
  //if(TMath::Abs(avtx->GetX())>10.0) return kFALSE;
  //if(TMath::Abs(avtx->GetY())>10.0) return kFALSE;
  if(TMath::Abs(avtx->GetZ())>10.0) return kFALSE; // TBI setter
  if(avtx->GetNContributors()<=2) return kFALSE; // TBI setter
 }

 return kTRUE;

} // Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesMixedEventCuts(AliVEvent *ave)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::CalculateCorrelationFunctions(AliAODEvent *aAOD)
{
 // Calculate correlation functions.

 // a) Insanity checks;
 // b) Two nested loops to calculate C(k), just an example.

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::CalculateCorrelationFunctions(AliAODEvent *aAOD)";
 if(!aAOD){Fatal(sMethodName.Data(),"!aAOD");}
 if(0 == fGlobalTracksAOD[0]->GetSize()){GlobalTracksAOD(aAOD,0);}
 if(0 == fGlobalTracksAOD[0]->GetSize()){return;}

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
  if(!PassesCommonGlobalTrackCuts(gtrack1)){continue;}

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
   if(!PassesCommonGlobalTrackCuts(gtrack2)){continue;}

   // Okay, so we have two tracks, let's check PID, and fill the correlation functions:

   // 1.) Same particle species:

   // a) pion-pion:
   //  a1) pi+pi+ [2][2]:
   if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE))
   {
    fCorrelationFunctions[2][2]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  a2) pi-pi- [7][7]:
   if(Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[7][7]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  a3) pi+pi- || pi-pi+ [2][7]:
   if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,1,kTRUE)))
   {
    fCorrelationFunctions[2][7]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }

   // b) kaon-kaon:
   //  b1) K+K+ [3][3]:
   if(Kaon(gtrack1,1,kTRUE) && Kaon(gtrack2,1,kTRUE))
   {
    fCorrelationFunctions[3][3]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  b2) K-K- [8][8]:
   if(Kaon(gtrack1,-1,kTRUE) && Kaon(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[8][8]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  b3) K+K- || K-K+ [3][8]:
   if((Kaon(gtrack1,1,kTRUE) && Kaon(gtrack2,-1,kTRUE)) || (Kaon(gtrack1,-1,kTRUE) && Kaon(gtrack2,1,kTRUE)))
   {
    fCorrelationFunctions[3][8]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }

   // c) proton-proton:
   //  c1) p+p+ [4][4]:
   if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE))
   {
    fCorrelationFunctions[4][4]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  c2) p-p- [9][9]:
   if(Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[9][9]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  c3) p+p- || p-p+ [4][9]:
   if((Proton(gtrack1,1,kTRUE) && Proton(gtrack2,-1,kTRUE)) || (Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,1,kTRUE)))
   {
    fCorrelationFunctions[4][9]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }

   // 2.) Mixed particle species:
   // a) pion-kaon
   //  a1) pi+K+ [2][3]:
   if(Pion(gtrack1,1,kTRUE) && Kaon(gtrack2,1,kTRUE))
   {
    fCorrelationFunctions[2][3]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  a2) pi+K- [2][8]
   if(Pion(gtrack1,1,kTRUE) && Kaon(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[2][8]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  a3) K+pi- [3][7]
   if(Kaon(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[3][7]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  a4) pi-K- [7][8]
   if(Pion(gtrack1,-1,kTRUE) && Kaon(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[7][8]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   // b) pion-proton
   //  b1) pi+p+ [2][4]:
   if(Pion(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE))
   {
    fCorrelationFunctions[2][4]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  b2) pi+p- [2][9]
   if(Pion(gtrack1,1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[2][9]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  b3) p+pi- [4][7]
   if(Proton(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[4][7]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  b4) pi-p- [7][9]
   if(Pion(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[7][9]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   // c) kaon-proton
   //  c1) K+p+ [3][4]:
   if(Kaon(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE))
   {
    fCorrelationFunctions[3][4]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  c2) K+p- [3][9]
   if(Kaon(gtrack1,1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[3][9]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  c3) p+K- [4][8]
   if(Proton(gtrack1,1,kTRUE) && Kaon(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[4][8]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  c4) K-p- [8][9]
   if(Kaon(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    fCorrelationFunctions[8][9]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
  } // for(Int_t iTrack2=0;iTrack2<nTracks;iTrack2++)
 } // for(Int_t iTrack1=0;iTrack1<nTracks;iTrack1++)

} // void AliAnalysisTaskMultiparticleFemtoscopy::CalculateCorrelationFunctions()

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

void AliAnalysisTaskMultiparticleFemtoscopy::CalculateBackground(TClonesArray *ca1, TClonesArray *ca2)
{
 // Calculate background.

 // TBI this method is not really validated

 // a) Insanity checks;
 // b) Temporary primitive algorithm: correlate particles from 'previous event' to particles from 'current event';
 // c) Two nested loops to calculate B(k);
 // d) Shift [1] -> [0];
 // e) Clean [1].

 // a) Insanity checks:
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::CalculateBackground(TClonesArray *ca1, TClonesArray *ca2)";
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
  if(!PassesCommonGlobalTrackCuts(gtrack1)){continue;}

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
   if(!PassesCommonGlobalTrackCuts(gtrack2)){continue;}

   // Okay... So we have two tracks from two different events ready. Let's check PID, and calculate the background:

   // 1.) Same particle species:

   // a) pion-pion:
   //  a1) pi+pi+ [2][2]:
   if(Pion(gtrack1,1,kTRUE) && Pion(gtrack2,1,kTRUE))
   {
    fBackground[2][2]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  a2) pi-pi- [7][7]:
   if(Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,-1,kTRUE))
   {
    fBackground[7][7]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  a3) pi+pi- || pi-pi+ [2][7]:
   if((Pion(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE)) || (Pion(gtrack1,-1,kTRUE) && Pion(gtrack2,1,kTRUE)))
   {
    fBackground[2][7]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }

   // b) kaon-kaon:
   //  b1) K+K+ [3][3]:
   if(Kaon(gtrack1,1,kTRUE) && Kaon(gtrack2,1,kTRUE))
   {
    fBackground[3][3]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  b2) K-K- [8][8]:
   if(Kaon(gtrack1,-1,kTRUE) && Kaon(gtrack2,-1,kTRUE))
   {
    fBackground[8][8]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  b3) K+K- || K-K+ [3][8]:
   if((Kaon(gtrack1,1,kTRUE) && Kaon(gtrack2,-1,kTRUE)) || (Kaon(gtrack1,-1,kTRUE) && Kaon(gtrack2,1,kTRUE)))
   {
    fBackground[3][8]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }

   // c) proton-proton:
   //  c1) p+p+ [4][4]:
   if(Proton(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE))
   {
    fBackground[4][4]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  c2) p-p- [9][9]:
   if(Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    fBackground[9][9]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  c3) p+p- || p-p+ [4][9]:
   if((Proton(gtrack1,1,kTRUE) && Proton(gtrack2,-1,kTRUE)) || (Proton(gtrack1,-1,kTRUE) && Proton(gtrack2,1,kTRUE)))
   {
    fBackground[4][9]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }

   // 2.) Mixed particle species:
   // a) pion-kaon
   //  a1) pi+K+ [2][3]:
   if(Pion(gtrack1,1,kTRUE) && Kaon(gtrack2,1,kTRUE))
   {
    fBackground[2][3]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  a2) pi+K- [2][8]
   if(Pion(gtrack1,1,kTRUE) && Kaon(gtrack2,-1,kTRUE))
   {
    fBackground[2][8]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  a3) K+pi- [3][7]
   if(Kaon(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE))
   {
    fBackground[3][7]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  a4) pi-K- [7][8]
   if(Pion(gtrack1,-1,kTRUE) && Kaon(gtrack2,-1,kTRUE))
   {
    fBackground[7][8]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   // b) pion-proton
   //  b1) pi+p+ [2][4]:
   if(Pion(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE))
   {
    fBackground[2][4]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  b2) pi+p- [2][9]
   if(Pion(gtrack1,1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    fBackground[2][9]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  b3) p+pi- [4][7]
   if(Proton(gtrack1,1,kTRUE) && Pion(gtrack2,-1,kTRUE))
   {
    fBackground[4][7]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  b4) pi-p- [7][9]
   if(Pion(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    fBackground[7][9]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   // c) kaon-proton
   //  c1) K+p+ [3][4]:
   if(Kaon(gtrack1,1,kTRUE) && Proton(gtrack2,1,kTRUE))
   {
    fBackground[3][4]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  c2) K+p- [3][9]
   if(Kaon(gtrack1,1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    fBackground[3][9]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  c3) p+K- [4][8]
   if(Proton(gtrack1,1,kTRUE) && Kaon(gtrack2,-1,kTRUE))
   {
    fBackground[4][8]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }
   //  c4) K-p- [8][9]
   if(Kaon(gtrack1,-1,kTRUE) && Proton(gtrack2,-1,kTRUE))
   {
    fBackground[8][9]->Fill(RelativeMomenta(gtrack1,gtrack2));
   }

  } // for(Int_t iTrack2=0;iTrack2<nTracks2;iTrack2++)
 } // for(Int_t iTrack1=0;iTrack1<nTracks1;iTrack1++)

 // d) Shift [1] -> [0]:
 // TBI re-think the lines below, there should be a better way...
 fMixedEvents[0] = (TClonesArray*)fMixedEvents[1]->Clone();
 fGlobalTracksAOD[1] = (TExMap*)fGlobalTracksAOD[2]->Clone();

 // e) Clean [1]:
 fMixedEvents[1] = NULL;
 fGlobalTracksAOD[2]->Delete(); // TBI or = NULL ?

} // void AliAnalysisTaskMultiparticleFemtoscopy::CalculateBackground(TClonesArray *ca1, TClonesArray *ca2)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::CalculateBackground(TClonesArray *ca1, TClonesArray *ca2, Bool_t bMC)
{
 // Calculate background for MC events.

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
 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::CalculateBackground(TClonesArray *ca1, TClonesArray *ca2, Bool_t bMC)";
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
   fBackground[index1][index2]->Fill(RelativeMomenta(amcparticle1,amcparticle2)); // for the relative momenta ordering doesn't matter

  } // for(Int_t iTrack2=0;iTrack2<nTracks2;iTrack2++)
 } // for(Int_t iTrack1=0;iTrack1<nTracks1;iTrack1++)

 // d) Shift [1] -> [0]:
 // TBI re-think the lines below, there should be a better way...
 fMixedEvents[0] = (TClonesArray*)fMixedEvents[1]->Clone();

 // e) Clean [1]:
 fMixedEvents[1] = NULL;

} // void AliAnalysisTaskMultiparticleFemtoscopy::CalculateBackground(TClonesArray *ca1, TClonesArray *ca2, Bool_t bMC)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::GetOutputHistograms(TList *histList)
{
 // Get pointers for everything in the base list "fHistList".

 // a) Get pointer for base list fHistList;
 // b) Get pointer for profile holding internal flags and set again all flags TBI this profile is not implemented yet;
 // *) TBI ...
 // *) Get pointers for correlation functions;
 // *) Get pointers for background.

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

} // void AliAnalysisTaskMultiparticleFemtoscopy::GetOutputHistograms(TList *histList)

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::GetPointersForCorrelationFunctions()
{
 // Get pointers for correlation functions.

 // a) Get pointer for fCorrelationFunctionsList;
 // b) Get pointer for fCorrelationFunctionsFlagsPro;
 // c) Set again all flags;
 // d) Get pointers for fCorrelationFunctions[10][10].

 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::GetPointersForCorrelationFunctions()";

 // a) Get pointer for fCorrelationFunctionsList:
 fCorrelationFunctionsList = dynamic_cast<TList*>(fHistList->FindObject("Correlation_Functions"));
 if(!fCorrelationFunctionsList){Fatal(sMethodName.Data(),"fCorrelationFunctionsList");}

 // b) Get pointer for fCorrelationFunctionsFlagsPro:
 fCorrelationFunctionsFlagsPro = dynamic_cast<TProfile*>(fCorrelationFunctionsList->FindObject("fCorrelationFunctionsFlagsPro"));
 if(!fCorrelationFunctionsFlagsPro){Fatal(sMethodName.Data(),"fCorrelationFunctionsFlagsPro");}

 // c) Set again all flags:
 fFillCorrelationFunctions = (Bool_t) fCorrelationFunctionsFlagsPro->GetBinContent(1);
 fNormalizeCorrelationFunctions = (Bool_t) fCorrelationFunctionsFlagsPro->GetBinContent(2);

 if(!fFillCorrelationFunctions){return;} // TBI is this safe enough

 // d) Get pointers for fCorrelationFunctions[10][10]:
 const Int_t nParticleSpecies = 5; // Supported at the moment: 0=e,1=mu,2=pi,3=K,4=p
 for(Int_t pid1=0;pid1<2*nParticleSpecies;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pid2=pid1;pid2<2*nParticleSpecies;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   fCorrelationFunctions[pid1][pid2] = dynamic_cast<TH1F*>(fCorrelationFunctionsList->FindObject(Form("fCorrelationFunctions[%d][%d]",pid1,pid2)));
   if(!fCorrelationFunctions[pid1][pid2]){Fatal(sMethodName.Data(),"fCorrelationFunctions[%d][%d]",pid1,pid2);}
  } // for(Int_t pid2=0;pid2<5;pid2++) // anti-particle [0=e,1=mu,2=pi,3=K,4=p]
 } // for(Int_t pid=0;pid<5;pid++) // particle [0=e,1=mu,2=pi,3=K,4=p]

} // void AliAnalysisTaskMultiparticleFemtoscopy::GetPointersForCorrelationFunctions()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::GetPointersForBackground()
{
 // Get pointers for background.

 // a) Get pointer for fBackgroundList;
 // b) Get pointer for fBackgroundFlagsPro;
 // c) Set again all flags;
 // d) Get pointers for fBackground[10][10].

 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::GetPointersForBackground()";

 // a) Get pointer for fBackgroundList:
 fBackgroundList = dynamic_cast<TList*>(fHistList->FindObject("Background"));
 if(!fBackgroundList){Fatal(sMethodName.Data(),"fBackgroundList");}

 // b) Get pointer for fBackgroundFlagsPro:
 fBackgroundFlagsPro = dynamic_cast<TProfile*>(fBackgroundList->FindObject("fBackgroundFlagsPro"));
 if(!fBackgroundFlagsPro){Fatal(sMethodName.Data(),"fBackgroundFlagsPro");}

 // c) Set again all flags:
 fEstimateBackground = (Bool_t) fBackgroundFlagsPro->GetBinContent(1);

 if(!fEstimateBackground){return;} // TBI is this safe enough

 // d) Get pointers for fBackground[10][10]:
 const Int_t nParticleSpecies = 5; // Supported at the moment: 0=e,1=mu,2=pi,3=K,4=p
 for(Int_t pid1=0;pid1<2*nParticleSpecies;pid1++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pid2=pid1;pid2<2*nParticleSpecies;pid2++) // [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p]
  {
   fBackground[pid1][pid2] = dynamic_cast<TH1F*>(fBackgroundList->FindObject(Form("fBackground[%d][%d]",pid1,pid2)));
   if(!fBackground[pid1][pid2]){Fatal(sMethodName.Data(),"fBackground[%d][%d]",pid1,pid2);}
  } // for(Int_t pid2=0;pid2<5;pid2++) // anti-particle [0=e,1=mu,2=pi,3=K,4=p]
 } // for(Int_t pid=0;pid<5;pid++) // particle [0=e,1=mu,2=pi,3=K,4=p]

} // void AliAnalysisTaskMultiparticleFemtoscopy::GetPointersForBackground()

//=======================================================================================================================


