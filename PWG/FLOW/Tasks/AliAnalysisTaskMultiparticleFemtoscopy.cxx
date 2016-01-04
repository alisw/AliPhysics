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

#include "TCanvas.h" // TBI
#include "TFile.h" // TBI

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskMultiparticleFemtoscopy)

//================================================================================================================

AliAnalysisTaskMultiparticleFemtoscopy::AliAnalysisTaskMultiparticleFemtoscopy(const char *name, Bool_t useParticleWeights): 
 AliAnalysisTaskSE(name), 
 fHistList(NULL),
 fPIDResponse(NULL),
 fGlobalTracksAOD(NULL),
 fUseInternalFlags(kFALSE),
 // 1.) Control histograms:
 fControlHistogramsList(NULL),
 fControlHistogramsFlagsPro(NULL),
 fFillControlHistograms(kFALSE),
 fControlHistogramsEventList(NULL),
 fControlHistogramsEventFlagsPro(NULL),
 fFillControlHistogramsEvent(kFALSE),
 fGetNumberOfTracksHist(NULL),
 fGetNumberOfV0sHist(NULL),
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
 // *.) Debugging:
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
  this->InitializeArraysForControlHistograms();
  this->InitializeArraysForEBEObjects();

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
  
  fPIDResponse = new AliPIDResponse();
  fGlobalTracksAOD = new TExMap();

  DefineOutput(1, TList::Class());  

} // AliAnalysisTaskMultiparticleFemtoscopy::AliAnalysisTaskMultiparticleFemtoscopy(const char *name, Bool_t useParticleWeights): 

//================================================================================================================

AliAnalysisTaskMultiparticleFemtoscopy::AliAnalysisTaskMultiparticleFemtoscopy(): 
 AliAnalysisTaskSE(),
 fHistList(NULL),
 fPIDResponse(NULL),
 fGlobalTracksAOD(NULL),
 fUseInternalFlags(kFALSE),
 // 1.) Control histograms:
 fControlHistogramsList(NULL),
 fControlHistogramsFlagsPro(NULL),
 fFillControlHistograms(kFALSE),
 fControlHistogramsEventList(NULL),
 fControlHistogramsEventFlagsPro(NULL),
 fFillControlHistogramsEvent(kFALSE),
 fGetNumberOfTracksHist(NULL),
 fGetNumberOfV0sHist(NULL),
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
 // *.) Debugging:
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
 if(fGlobalTracksAOD) delete fGlobalTracksAOD;

} // AliAnalysisTaskMultiparticleFemtoscopy::~AliAnalysisTaskMultiparticleFemtoscopy()

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::UserCreateOutputObjects() 
{
 // Called at every worker node to initialize.

 // a) Trick to avoid name clashes, part 1; 
 // b) Book and nest all lists;
 // c) Book all objects;
 // d) Set all flags;
 // e) Trick to avoid name clashes, part 2. 
  
 // a) Trick to avoid name clashes, part 1: 
 Bool_t oldHistAddStatus = TH1::AddDirectoryStatus(); 
 TH1::AddDirectory(kFALSE);

 // b) Book and nest all lists:
 this->BookAndNestAllLists();

 // c) Book all objects:
 this->BookEverythingForControlHistograms();
 this->BookEverythingForEBEObjects();

 // d) Set all flags:
 // ...

 // e) Trick to avoid name clashes, part 2:
 TH1::AddDirectory(oldHistAddStatus);

 PostData(1,fHistList);

} // void AliAnalysisTaskMultiparticleFemtoscopy::UserCreateOutputObjects() 

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::UserExec(Option_t *) 
{
 // Main loop (called for each event).

 // *) Insanity checks;
 // *) Reset event-by-event objects;

 //fEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(0));

 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::UserExec(Option_t *)";
 
 // *) Insanity checks:
 InsanityChecksUserExec();

 // ... TBI
 AliMCEvent *aMC = MCEvent();                                  // from TaskSE
 AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(InputEvent()); // from TaskSE
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(InputEvent()); // from TaskSE

 if(aMC)
 {
  cout<<"aMC"<<endl;  
 }
 else if(aESD)
 {
  cout<<"aESD"<<endl;
 }
 else if(aAOD)
 {
  // Debugging:
  if(fWaitForSpecifiedEvent)
  {
   cout<<Form("aAOD->GetRunNumber() = %d",aAOD->GetRunNumber())<<endl;
   cout<<Form("aAOD->GetBunchCrossNumber() = %d",aAOD->GetBunchCrossNumber())<<endl;
   cout<<Form("aAOD->GetOrbitNumber() = %d",aAOD->GetOrbitNumber())<<endl;
   cout<<Form("aAOD->GetPeriodNumber() = %d",aAOD->GetPeriodNumber())<<endl;
   if(!SpecifiedEvent(aAOD->GetRunNumber(),aAOD->GetBunchCrossNumber(),aAOD->GetOrbitNumber(),aAOD->GetPeriodNumber())){return;}
  } // if(fWaitForSpecifiedEvent)

  // Filter out normal global tracks:
  GlobalTracksAOD(aAOD);
  if(0 == fGlobalTracksAOD->GetSize()) return;

  // Common event selection criteria:
  if(!PassesCommonEventCuts(aAOD)){return;}

  // AOD event:
  fGetNumberOfTracksHist->Fill(aAOD->GetNumberOfTracks()); // TBI not all tracks are unique
  fGetNumberOfV0sHist->Fill(aAOD->GetNumberOfV0s()); // TBI some V0s share the daughter

  // AOD primary vertex:
  AliAODVertex *avtx = (AliAODVertex*)aAOD->GetPrimaryVertex();
  fVertexXYZ[0]->Fill(avtx->GetX());
  fVertexXYZ[1]->Fill(avtx->GetY());
  fVertexXYZ[2]->Fill(avtx->GetZ());
  fGetNContributorsHist->Fill(avtx->GetNContributors());
  fGetChi2perNDFHist->Fill(avtx->GetChi2perNDF());
  fGetNDaughtersHist->Fill(avtx->GetNDaughters());
 
  // AOD tracks:
  for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++)
  {
   AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack));
   if(!atrack){Fatal(sMethodName.Data(),"!atrack");} // TBI keep this for some time, eventually just continue
   if(0==atrack->GetFilterMap()){continue;} // TBI add comment

   // Corresponding AOD global track:
   Int_t id = atrack->GetID();
   AliAODTrack *gtrack = dynamic_cast<AliAODTrack*>(id>=0 ? aAOD->GetTrack(fGlobalTracksAOD->GetValue(id)) : aAOD->GetTrack(fGlobalTracksAOD->GetValue(-(id+1)))); 
   if(!gtrack){Fatal(sMethodName.Data(),"!gtrack");} // TBI keep this for some time, eventually just continue

   // Common track selection critera:
   if(!PassesCommonTrackCuts(gtrack)){continue;} // TBI I am applying them to global tracks... rethink

   // Fill the histograms:
   fChargeHist->Fill(gtrack->Charge()+0.5); // see how this histogram was booked for convention used
   fGetTPCNclsHist->Fill(gtrack->GetTPCNcls());
   fGetTPCsignalNHist->Fill(gtrack->GetTPCsignalN());
   fGetITSNclsHist->Fill(gtrack->GetITSNcls());
   fdEdxVsPtHist->Fill(gtrack->GetTPCmomentum(),gtrack->GetTPCsignal());
   fPtHist->Fill(gtrack->Pt());
   fEtaHist->Fill(gtrack->Eta());
   fPhiHist->Fill(gtrack->Phi());
   fMassHist->Fill(gtrack->M());
   
   // PID:
   // Protons:
   if(Proton(gtrack,1,kTRUE))
   {
    fPtPIDHist[4][0][0]->Fill(gtrack->Pt());
    fMassPIDHist[4][0][0]->Fill(gtrack->M());
    fEtaPIDHist[4][0][0]->Fill(gtrack->Eta());
    fPhiPIDHist[4][0][0]->Fill(gtrack->Phi());
   } 
   else if(Proton(gtrack,1,kFALSE))
   {
    fPtPIDHist[4][0][1]->Fill(gtrack->Pt());
    fMassPIDHist[4][0][1]->Fill(gtrack->M());
    fEtaPIDHist[4][0][1]->Fill(gtrack->Eta());
    fPhiPIDHist[4][0][1]->Fill(gtrack->Phi());
   } 
   else if(Proton(gtrack,-1,kTRUE))
   {
    fPtPIDHist[4][1][0]->Fill(gtrack->Pt());
    fMassPIDHist[4][1][0]->Fill(gtrack->M());
    fEtaPIDHist[4][1][0]->Fill(gtrack->Eta());
    fPhiPIDHist[4][1][0]->Fill(gtrack->Phi());
   } 
   else if(Proton(gtrack,-1,kFALSE))
   {
    fPtPIDHist[4][1][1]->Fill(gtrack->Pt());
    fMassPIDHist[4][1][1]->Fill(gtrack->M());
    fEtaPIDHist[4][1][1]->Fill(gtrack->Eta());
    fPhiPIDHist[4][1][1]->Fill(gtrack->Phi());
   }          

   // Pions:
   if(Pion(gtrack,1,kTRUE))
   {
    fPtPIDHist[2][0][0]->Fill(gtrack->Pt());
    fMassPIDHist[2][0][0]->Fill(gtrack->M());
    fEtaPIDHist[2][0][0]->Fill(gtrack->Eta());
    fPhiPIDHist[2][0][0]->Fill(gtrack->Phi());
   } 
   else if(Pion(gtrack,1,kFALSE))
   {
    fPtPIDHist[2][0][1]->Fill(gtrack->Pt());
    fMassPIDHist[2][0][1]->Fill(gtrack->M());
    fEtaPIDHist[2][0][1]->Fill(gtrack->Eta());
    fPhiPIDHist[2][0][1]->Fill(gtrack->Phi());
   } 
   else if(Pion(gtrack,-1,kTRUE))
   {
    fPtPIDHist[2][1][0]->Fill(gtrack->Pt());
    fMassPIDHist[2][1][0]->Fill(gtrack->M());
    fEtaPIDHist[2][1][0]->Fill(gtrack->Eta());
    fPhiPIDHist[2][1][0]->Fill(gtrack->Phi());
   } 
   else if(Pion(gtrack,-1,kFALSE))
   {
    fPtPIDHist[2][1][1]->Fill(gtrack->Pt());
    fMassPIDHist[2][1][1]->Fill(gtrack->M());
    fEtaPIDHist[2][1][1]->Fill(gtrack->Eta());
    fPhiPIDHist[2][1][1]->Fill(gtrack->Phi());
   }          

   // Kaons:
   if(Kaon(gtrack,1,kTRUE))
   {
    fPtPIDHist[3][0][0]->Fill(gtrack->Pt());
    fMassPIDHist[3][0][0]->Fill(gtrack->M());
    fEtaPIDHist[3][0][0]->Fill(gtrack->Eta());
    fPhiPIDHist[3][0][0]->Fill(gtrack->Phi());
   } 
   else if(Kaon(gtrack,1,kFALSE))
   {
    fPtPIDHist[3][0][1]->Fill(gtrack->Pt());
    fMassPIDHist[3][0][1]->Fill(gtrack->M());
    fEtaPIDHist[3][0][1]->Fill(gtrack->Eta());
    fPhiPIDHist[3][0][1]->Fill(gtrack->Phi());
   } 
   else if(Kaon(gtrack,-1,kTRUE))
   {
    fPtPIDHist[3][1][0]->Fill(gtrack->Pt());
    fMassPIDHist[3][1][0]->Fill(gtrack->M());
    fEtaPIDHist[3][1][0]->Fill(gtrack->Eta());
    fPhiPIDHist[3][1][0]->Fill(gtrack->Phi());
   } 
   else if(Kaon(gtrack,-1,kFALSE))
   {
    fPtPIDHist[3][1][1]->Fill(gtrack->Pt());
    fMassPIDHist[3][1][1]->Fill(gtrack->M());
    fEtaPIDHist[3][1][1]->Fill(gtrack->Eta());
    fPhiPIDHist[3][1][1]->Fill(gtrack->Phi());
   }          

  } // for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++)
  

  // AOD V0s:
  TClonesArray *caV0s = aAOD->GetV0s(); 
  Int_t index = 0;
  AliAODv0 *aAODv0 = NULL;
  Int_t nProngs = -44;
  while(caV0s->At(index))
  {
   aAODv0 = (AliAODv0*) caV0s->At(index++);
   if(!aAODv0) break;


   AliAODv0 *temp = (AliAODv0*)fPIDV0sCA[0]->ConstructedAt(index-1);
   temp = (AliAODv0*)aAODv0->Clone();
   
   cout<<Form("index = %d, fPIDV0sCA[0]->GetEntries() = %d",index,fPIDV0sCA[0]->GetEntries())<<endl;

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

   cout<<Form("V0PosID: %d , V0NegID: %d",aAODv0->GetPosID(),aAODv0->GetNegID())<<endl;
   Int_t trackPos = aAODv0->GetPosID()>=0 ? fGlobalTracksAOD->GetValue(aAODv0->GetPosID()) : fGlobalTracksAOD->GetValue(-(aAODv0->GetPosID()+1));
   Int_t trackNeg = aAODv0->GetNegID()>=0 ? fGlobalTracksAOD->GetValue(aAODv0->GetNegID()) : fGlobalTracksAOD->GetValue(-(aAODv0->GetNegID()+1));
   cout<<Form("global : %d, global : %d", trackPos, trackNeg)<<endl; 
   

   if(-1 != fUniqueIDHistEBE->FindFirstBinAbove(1.44,1)) // TBI
   {
    cout<<Form("fUniqueIDHistEBE->FindFirstBinAbove(1.44) %d:",(Int_t)fUniqueIDHistEBE->FindFirstBinAbove(1.44,1))<<endl; // 
   }
  } // while(caV0s->At(index))


  Int_t index2 = 0;
  while(caV0s->At(index2))
  {
   cout<<((AliAODv0*) caV0s->At(index2))->Pt()<<endl;
   cout<<((AliAODv0*) fPIDV0sCA[0]->At(index2))->Pt()<<endl;
   index2++;
   cout<<endl;
  }


   //   AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,globaltrack);
   //   AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,globaltrack);

 } 
 else
 {
  exit(0);
 }



 // *) Reset event-by-event objects:
 this->ResetEBEObjects();

 PostData(1,fHistList);

 if( 0 == ((Int_t)fGetNumberOfTracksHist->GetEntries())%10000 )
 {
  cout<<Form("nEvts: %d",(Int_t)fGetNumberOfTracksHist->GetEntries())<<endl;
  TFile *f = new TFile("AnalysisResults.root","RECREATE"); // TBI remove eventually
  fHistList->Write(fHistList->GetName(),TObject::kSingleKey); // TBI remove eventually
  f->Close(); // TBI remove eventually
 }

} // void AliAnalysisTaskMultiparticleFemtoscopy::UserExec(Option_t *) 

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::Terminate(Option_t *) 
{
 // Accessing the merged output list. 

 fHistList = (TList*)GetOutputData(1);

 if(!fHistList){exit(1);}

 //TDirectoryFile *df = new TDirectoryFile("outputMPFanalysis","");
 //df->Add(fHistList);

 TFile *f = new TFile("AnalysisResults.root","RECREATE");
 fHistList->Write(fHistList->GetName(),TObject::kSingleKey);

 delete f;

} // end of void AliAnalysisTaskMultiparticleFemtoscopy::Terminate(Option_t *)

//================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::InsanityChecksUserExec()
{
 // Insanity...

 TString sMethodName = "void AliAnalysisTaskMultiparticleFemtoscopy::InsanityChecksUserExec()";

 if(!fGlobalTracksAOD){Fatal(sMethodName.Data(),"!fGlobalTracksAOD");}
 if(0 != fGlobalTracksAOD->GetSize()){Fatal(sMethodName.Data(),"0 != fGlobalTracksAOD->GetSize()");}

 for(Int_t pid=0;pid<5;pid++) // [0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pa=0;pa<2;pa++) // particle/antiparticle
  {
   for(Int_t ps=0;ps<2;ps++) // kPrimary/kFromDecayVtx
   {
    if(!fPIDCA[pid][pa][ps]) {Fatal(sMethodName.Data(),"!fPIDCA[pid][pa][ps]");}
    if(0 != fPIDCA[pid][pa][ps]->GetEntriesFast()){Fatal(sMethodName.Data(),"0 != fPIDCA[pid][pa][ps]->GetEntriesFast()"); }
   }
  }
 }

 for(Int_t pid=0;pid<1;pid++) // [0=Lambda,1=...] 
 {
  if(!fPIDV0sCA[pid]){Fatal(sMethodName.Data(),"!fPIDV0sCA[pid]");}
  if(0 != fPIDV0sCA[pid]->GetEntriesFast()){Fatal(sMethodName.Data(),"0 != fPIDV0sCA[pid]->GetEntriesFast()");}  
 }

} // void AliAnalysisTaskMultiparticleFemtoscopy::InsanityChecksUserExec()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::ResetEBEObjects()
{
 // Reset all event-by-event objects.

 if(fUniqueIDHistEBE) fUniqueIDHistEBE->Reset();
 if(fGlobalTracksAOD) fGlobalTracksAOD->Delete();

 // TBI add comment
 for(Int_t pid=0;pid<5;pid++) // [0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pa=0;pa<2;pa++) // particle/antiparticle
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

} // void AliAnalysisTaskMultiparticleFemtoscopy::ResetEBEObjects()

//=======================================================================================================================

Bool_t AliAnalysisTaskMultiparticleFemtoscopy::Pion(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)
{
 // Is this primary/secondary (anti-)pion?

 // a) Insanity checks;
 // b) Trivial checks; // TBI
 // c) Track quality cuts;
 // d) PID. 

 Double_t inc = 1.0; Double_t exl = 3.0; // TBI to be removed eventually

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
  if(!(dSigmaTPCPion < inc && dSigmaTPCProton > exl && dSigmaTPCKaon > exl)) return kFALSE; // TBI hardwired 3. and 4.
 } // if(gtrack->GetTPCmomentum() <= 0.75) // TBI hardwired 0.75
 else if(gtrack->GetTPCmomentum() >= 0.75 )
 {
  // TBI use combined TPC and TOf
  return kFALSE; // TBI remove eventually
 }

 return kTRUE;

} // Bool_t AliAnalysisTaskMultiparticleFemtoscopy::Pion(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)


//=======================================================================================================================

Bool_t AliAnalysisTaskMultiparticleFemtoscopy::Kaon(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)
{
 // Is this primary/secondary (anti-)kaon?

 // a) Insanity checks;
 // b) Trivial checks;
 // c) Track quality cuts;
 // d) PID. 


 Double_t inc = 1.0; Double_t exl = 3.0; // TBI to be removed eventually


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
  if(!(dSigmaTPCKaon < inc && dSigmaTPCProton > exl && dSigmaTPCPion > exl)) return kFALSE; // TBI hardwired 3. and 4.
 } // if(gtrack->GetTPCmomentum() <= 0.75) // TBI hardwired 0.75
 else if(gtrack->GetTPCmomentum() >= 0.75 )
 {
  // TBI use combined TPC and TOf
  return kFALSE; // TBI remove eventually
 }

 return kTRUE;

} // Bool_t AliAnalysisTaskMultiparticleFemtoscopy::Kaon(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)

//=======================================================================================================================

Bool_t AliAnalysisTaskMultiparticleFemtoscopy::Proton(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)
{
 // Is this primary/secondary (anti-)proton?

 // a) Insanity checks;
 // b) Trivial checks;
 // c) Track quality cuts;
 // d) PID. 

 
 Double_t inc = 1.0; Double_t exl = 3.0; // TBI to be removed eventually


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
  if(!(dSigmaTPCProton < inc && dSigmaTPCPion > exl && dSigmaTPCKaon > exl)) return kFALSE; // TBI hardwired 3. and 4.
 } // if(gtrack->GetTPCmomentum() <= 0.75) // TBI hardwired 0.75
 else if(gtrack->GetTPCmomentum() >= 0.75 )
 {
  // TBI use combined TPC and TOf
  return kFALSE; // TBI remove eventually
 }

 return kTRUE;

} // Bool_t AliAnalysisTaskMultiparticleFemtoscopy::Proton(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)

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
  for(Int_t pa=0;pa<2;pa++) // particle/antiparticle
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

} // void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForControlHistograms()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::InitializeArraysForEBEObjects()
{
 // Initialize all arrays for e-b-e objects.

 for(Int_t pid=0;pid<5;pid++) // [0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pa=0;pa<2;pa++) // particle/antiparticle
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

void AliAnalysisTaskMultiparticleFemtoscopy::BookAndNestAllLists()
{
 // Book and nest all lists nested in the base list fHistList.

 // a) Book and nest lists for control histograms;
 // b) Book and nest lists for eveny-by-event histograms;
 // ...
 
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

} // void AliAnalysisTaskMultiparticleFemtoscopy::BookAndNestAllLists()

//=======================================================================================================================

void AliAnalysisTaskMultiparticleFemtoscopy::BookEverythingForControlHistograms()
{
 // Book all the stuff for control histograms.

 // a) Book the profile holding all the flags for control histograms;
 // b) Book all control histograms...
 //  b0) Event;
 //  b1) Non-identified particles;
 //  b2) Identified particles;
 //  b3) V0s. 

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

 //  b0) Event:
 // Book the profile holding all the flags for TBI:
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

 fGetNumberOfTracksHist = new TH1I("fGetNumberOfTracksHist","aAOD->GetNumberOfTracks()",10000,0,10000);
 //fGetNumberOfTracksHist->SetStats(kFALSE);
 fGetNumberOfTracksHist->SetFillColor(kBlue-10);
 fControlHistogramsEventList->Add(fGetNumberOfTracksHist);
 fGetNumberOfV0sHist = new TH1I("fGetNumberOfV0sHist","aAOD->GetNumberOfV0s()",10000,0,10000);
 fGetNumberOfV0sHist->SetStats(kFALSE);
 fGetNumberOfV0sHist->SetFillColor(kBlue-10);
 fControlHistogramsEventList->Add(fGetNumberOfV0sHist);
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

 //  b1) Non-identified particles:
 // Book the profile holding all the flags for TBI:
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
 fControlHistogramsNonIdentifiedParticlesList->Add(fGetTPCNclsHist);
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
 // 
 Double_t dNominalMass[5] = {TDatabasePDG::Instance()->GetParticle(11)->Mass(),TDatabasePDG::Instance()->GetParticle(13)->Mass(),TDatabasePDG::Instance()->GetParticle(211)->Mass(),TDatabasePDG::Instance()->GetParticle(321)->Mass(),TDatabasePDG::Instance()->GetParticle(2212)->Mass()};
 TString sParticleLabel[5] = {"e","#mu","#pi","K","p"};
 fMassHist = new TH1F("fMassHist","atrack->M()",10000,0.,10.);
 fMassHist->SetStats(kFALSE);
 fMassHist->SetFillColor(kBlue-10);
 fMassHist->SetMinimum(0.);
 for(Int_t nm=0;nm<5;nm++) // nominal masses 
 {
  fMassHist->GetXaxis()->SetBinLabel(fMassHist->FindBin(dNominalMass[nm]),Form("m_{%s}",sParticleLabel[nm].Data()));
 }
 fControlHistogramsNonIdentifiedParticlesList->Add(fMassHist);
 // ...

 //  b2) Identified particles:
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
 for(Int_t pid=0;pid<5;pid++) // [0=e,1=mu,2=pi,3=K,4=p]
 {
  for(Int_t pa=0;pa<2;pa++) // particle/antiparticle
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
    fEtaPIDHist[pid][pa][ps]->GetXaxis()->SetTitle("#eta");
    fEtaPIDHist[pid][pa][ps]->SetFillColor(kBlue-10);
    fControlHistogramsIdentifiedParticlesList->Add(fEtaPIDHist[pid][pa][ps]);
    fPhiPIDHist[pid][pa][ps] = new TH1F(Form("fPhiPIDHist[%d][%d][%d]",pid,pa,ps),Form("fPhiPIDHist[%d][%d][%d]",pid,pa,ps),360,0.,TMath::TwoPi());                      
    fPhiPIDHist[pid][pa][ps]->GetXaxis()->SetTitle("#phi");
    fPhiPIDHist[pid][pa][ps]->SetFillColor(kBlue-10);
    fControlHistogramsIdentifiedParticlesList->Add(fPhiPIDHist[pid][pa][ps]);
   } // for(Int_t ps=0;ps<2;ps++) // kPrimary/kFromDecayVtx
  } // for(Int_t pa=0;pa<2;pa++) // particle/antiparticle
 } // for(Int_t pid=0;pid<4;pid++) // [0=e,1=mu,2=pi,3=K,4=p]

 //  b3) V0s:
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
  for(Int_t pa=0;pa<2;pa++) // particle/antiparticle
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

void AliAnalysisTaskMultiparticleFemtoscopy::GlobalTracksAOD(AliAODEvent *aAOD)
{
 // Filter out global tracks in AOD and store them in fGlobalTracksAOD.
 
 // Remark 0: All global tracks have positive ID, the duplicated TPC-only tracks have -(ID+1);
 // Remark 1: The issue here is that there are apparently two sets of global tracks: a) "normal" and b) constrained to primary vertex. 
 //           The latter ones are 

 for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++) 
 {
  AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack));
  if(aodTrack)
  { 
   Int_t id = aodTrack->GetID();
   if(id>=0 && aodTrack->GetFilterMap()>0 && !aodTrack->IsGlobalConstrained()) // TBI rethink this 
   {
    fGlobalTracksAOD->Add(id,iTrack);
   } // if(id>=0 && aodTrack->GetFilterMap()>0)
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
Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesCommonTrackCuts(AliAODTrack *gtrack)
{
 // Check if the track passes common track cuts (irrespectively of PID).

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
 if(gtrack->GetTPCNcls()<70) return kFALSE; 
 if(gtrack->GetTPCsignalN()<70) return kFALSE; 

 return kTRUE; 

} // Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesCommonTrackCuts(AliAODTrack *gtrack)

//=======================================================================================================================

Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesCommonEventCuts(AliAODEvent *aAOD)
{
 // Check if the event passes common event cuts.

 // ...

 return kTRUE;

} // Bool_t AliAnalysisTaskMultiparticleFemtoscopy::PassesCommonEventCuts(AliAODEvent *aAOD)


//=======================================================================================================================
