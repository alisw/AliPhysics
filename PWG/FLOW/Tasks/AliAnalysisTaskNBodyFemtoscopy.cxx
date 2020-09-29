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

/************************************** 
* Femtoscopy task for N bodies * 
*Laura Serksnyte
*laura.serksnyte@cern.ch
**************************************/ 
  
#include "Riostream.h"
#include "AliAnalysisTaskNBodyFemtoscopy.h"
#include "AliLog.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "TFile.h"
#include "AliMultSelection.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskNBodyFemtoscopy)

//================================================================================================================

AliAnalysisTaskNBodyFemtoscopy::AliAnalysisTaskNBodyFemtoscopy(const char *name, Bool_t useParticleWeights): 
 AliAnalysisTaskSE(name), 
 fPIDResponse(NULL),
 fHistList(NULL),
 fRandomIndices(NULL),
 // Control histograms:
 fControlHistogramsList(NULL),
 fPtHist(NULL),
 fCentralityHist(NULL),
 fNumberOfTracksHist(NULL),
 fCentralityHistVzCut(NULL),
 fPtHistEtaCut(NULL),
 fPtHistEtaCutPTCut(NULL),
 fPtHistEtaCutPTCutPhiCut(NULL),
 fNumberOfTracksHistAfterAllCuts(NULL),
 fTestPIDTrueFalsePositive(NULL), 
 fTestTPCOnlyVsGlobal(NULL),
 fNbinsPt(1000),
 fMinBinPt(0.),
 fMaxBinPt(10.),
 fNbinsCentrality(100),
 fMinCentrality(0.),
 fMaxCentrality(100.),
 fNbinsMultiplicity(100.),
 fMinBinMultiplicity(0.),
 fMaxBinMultiplicity(20000.),
 fNbinsTPCOnlyVsGlobal(200),
 fMinBinTPCOnlyVsGlobal(-5.),
 fMaxBinTPCOnlyVsGlobal(5),
 fRejectEventsNoPrimaryVertex(kTRUE),
 fCutOnVertexZ(kTRUE),
 fApplyCommonTrackCuts(kTRUE),
 fUseDefaultInclusiveSigmaCuts(kTRUE),
 fUseDefaultExclusiveSigmaCuts(kTRUE),
 fRejectFakeTracks(kTRUE),
 fProcessBothKineAndReco(kFALSE),
 fCentralityEstimator("V0M"),
 fFilterBit(128),
 fMC(NULL),
 // Final results:
 fFinalResultsList(NULL)
 {
  // Constructor.
 
  AliDebug(2,"AliAnalysisTaskNBodyFemtoscopy::AliAnalysisTaskNBodyFemtoscopy(const char *name, Bool_t useParticleWeights)");

  // Base list:
  fHistList = new TList();
  fHistList->SetName("FemtoNBody");
  fHistList->SetOwner(kTRUE);

  // Initialize all arrays:
  this->InitializeArrays();

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
   // not needed for the time being
  }

  // c) Determine seed for gRandom:
  delete gRandom;
  gRandom = new TRandom3(0); 
} // AliAnalysisTaskNBodyFemtoscopy::AliAnalysisTaskNBodyFemtoscopy(const char *name, Bool_t useParticleWeights): 

//================================================================================================================

AliAnalysisTaskNBodyFemtoscopy::AliAnalysisTaskNBodyFemtoscopy(): 
 AliAnalysisTaskSE(), 
 fPIDResponse(NULL),
 fHistList(NULL),
 fRandomIndices(NULL),
 // Control histograms:
 fControlHistogramsList(NULL),
 fPtHist(NULL),
 fCentralityHist(NULL),
 fNumberOfTracksHist(NULL),
 fCentralityHistVzCut(NULL),
 fPtHistEtaCut(NULL),
 fPtHistEtaCutPTCut(NULL),
 fPtHistEtaCutPTCutPhiCut(NULL),
 fNumberOfTracksHistAfterAllCuts(NULL),
 fTestPIDTrueFalsePositive(NULL), 
 fTestTPCOnlyVsGlobal(NULL),
 fNbinsPt(1000),
 fMinBinPt(0.),
 fMaxBinPt(10.),
 fNbinsCentrality(100),
 fMinCentrality(0.),
 fMaxCentrality(100.),
 fNbinsMultiplicity(100.),
 fMinBinMultiplicity(0.),
 fMaxBinMultiplicity(20000.),
 fNbinsTPCOnlyVsGlobal(200),
 fMinBinTPCOnlyVsGlobal(-5.),
 fMaxBinTPCOnlyVsGlobal(5),
 fRejectEventsNoPrimaryVertex(kTRUE),
 fCutOnVertexZ(kTRUE),
 fApplyCommonTrackCuts(kTRUE),
 fUseDefaultInclusiveSigmaCuts(kTRUE),
 fUseDefaultExclusiveSigmaCuts(kTRUE),
 fRejectFakeTracks(kTRUE),
 fProcessBothKineAndReco(kFALSE),
 fCentralityEstimator("V0M"),
 fFilterBit(128),
 fMC(NULL),
 // Final results:
 fFinalResultsList(NULL)
{
  // Dummy constructor.
 
  AliDebug(2,"AliAnalysisTaskNBodyFemtoscopy::AliAnalysisTaskNBodyFemtoscopy()");
  this->InitializeArrays();

} // AliAnalysisTaskNBodyFemtoscopy::AliAnalysisTaskNBodyFemtoscopy():

//================================================================================================================

AliAnalysisTaskNBodyFemtoscopy::~AliAnalysisTaskNBodyFemtoscopy()
{
 // Destructor.

 if(fHistList) delete fHistList;
 if(fPIDResponse) delete fPIDResponse;
  
} // AliAnalysisTaskNBodyFemtoscopy::~AliAnalysisTaskNBodyFemtoscopy()

//================================================================================================================

void AliAnalysisTaskNBodyFemtoscopy::UserCreateOutputObjects() 
{
 // Called at every worker node to initialize.

 // a) Trick to avoid name clashes, part 1;
 // b) Book and nest all lists;
 // c) Book all objects;
 // *) Trick to avoid name clashes, part 2.
  
 // a) Trick to avoid name clashes, part 1:
 Bool_t oldHistAddStatus = TH1::AddDirectoryStatus(); 
 TH1::AddDirectory(kFALSE);

 // b) Book and nest all lists:
 this->BookAndNestAllLists();

 // c) Book all objects:
 this->BookControlHistograms();
 this->BookFinalResultsHistograms();
 this->BookEverything();
 // *) Trick to avoid name clashes, part 2:
 TH1::AddDirectory(oldHistAddStatus);

 PostData(1,fHistList);

} // void AliAnalysisTaskNBodyFemtoscopy::UserCreateOutputObjects() 

//================================================================================================================

void AliAnalysisTaskNBodyFemtoscopy::UserExec(Option_t *) 
{
 // Main loop (called for each event).
 // a) Get pointer to AOD event, chech multiplicity, apply event cuts;
 // b) Start analysis over AODs;
 // c) Reset event-by-event objects;
 // d) PostData.

 TString sMethodName = "void AliAnalysisTaskNBodyFemtoscopy::UserExec(Option_t *)";
 // a)  get MC event to check PID results
 fMC = MCEvent();
 if(!fMC){Fatal(sMethodName.Data(),"!fMC");}
 //Get pointer to AOD event:
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(InputEvent()); // from TaskSE
 if(!aAOD){Fatal(sMethodName.Data(),"!aAOD");}

 // create the randomizes indexes
 this->CreateRandomIndices(aAOD);


 // Check Multiplicity
 AliMultSelection *ams = (AliMultSelection*)aAOD->FindListObject("MultSelection");
 if(!ams){Fatal(sMethodName.Data(),"!ams");}
 Float_t eventCentrality = ams->GetMultiplicityPercentile(Form("%s", fCentralityEstimator.Data()));
 if(eventCentrality >= fMinCentrality && eventCentrality < fMaxCentrality)
 {
    fCentralityHist->Fill(eventCentrality);
 }
 else{ return;}
 // Apply other event cuts
 if(!this->CommonEventCuts(aAOD,eventCentrality)){return;}


 // b1) Start analysis over AODs:
 // Get the TExMap for GLOBAL tracks and the list of all tracks
 if(fGlobalTracksAODTEST[0]) fGlobalTracksAODTEST[0]->Delete();
 this->GlobalTracksAODTEST(aAOD,0);
 if(0 == fGlobalTracksAODTEST[0]->GetSize()){return;}
 *fAllTracksTEST[0] = *(aAOD->GetTracks());
 // Analyse track by track
 Int_t nTracks = aAOD->GetNumberOfTracks(); // number of all tracks in current event 
 fNumberOfTracksHist->Fill(nTracks);
 Int_t numberOfTracksAfterAllSelection = 0;
 for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
 {
  AliAODTrack *aTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack)); // getting a pointer to a track
  if(!aTrack){ continue;} // protection against NULL pointers
  if(aTrack->GetID()>=0 && aTrack->IsGlobalConstrained()){Fatal(sMethodName.Data(),"aTrack->GetID()>=0 && aTrack->IsGlobalConstrained()");} 
  if(fApplyCommonTrackCuts) {
   if(aTrack->TestFilterBit(128) && aTrack->IsGlobalConstrained()){Fatal(sMethodName.Data(),"aTrack->TestFiletrBit(128) && aTrack->IsGlobalConstrained()");}
  }
  if(!this->CommonTrackCuts(aTrack)){continue;}
  // get global track
  Int_t aTrackID = aTrack->GetID();
  AliAODTrack *gTrack = dynamic_cast<AliAODTrack*>(aTrackID>=0 ? fAllTracksTEST[0]->UncheckedAt(fGlobalTracksAODTEST[0]->GetValue(aTrackID)) : fAllTracksTEST[0]->UncheckedAt(fGlobalTracksAODTEST[0]->GetValue(-(aTrackID+1))));
  if(!gTrack){Fatal(sMethodName.Data(),"!gTrack");} 
  if(!this->GlobalTrackCuts(gTrack)){continue;}
  // check particle
  fTestTPCOnlyVsGlobal->Fill(gTrack->Pt()-aTrack->Pt());
  this->Pion(gTrack,1,kTRUE);
  this->Kaon(gTrack,1,kTRUE);
  this->Proton(gTrack,1,kTRUE);
  this->Pion(gTrack,-1,kTRUE);
  this->Kaon(gTrack,-1,kTRUE);
  this->Proton(gTrack,-1,kTRUE);
  Double_t pt = aTrack->Pt(); 
  fPtHist->Fill(pt); 
  numberOfTracksAfterAllSelection ++;
 } // for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
 fNumberOfTracksHistAfterAllCuts->Fill(numberOfTracksAfterAllSelection);
 
 // c) Reset event-by-event objects:
 // ...

 // d) PostData:
 PostData(1,fHistList);

} // void AliAnalysisTaskNBodyFemtoscopy::UserExec(Option_t *)




//================================================================================================================

Bool_t AliAnalysisTaskNBodyFemtoscopy::CommonEventCuts(AliVEvent *ave, Float_t centrality)
{
 // Apply event cuts.

 // a) Check which event type
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
  // a) Cuts on AliAODEvent:
  if(fRejectEventsNoPrimaryVertex && !aAOD->GetPrimaryVertex()) return kFALSE;

  // b) Cuts on AliAODVertex:
  AliAODVertex *avtx = (AliAODVertex*)aAOD->GetPrimaryVertex();
  if(fCutOnVertexZ)
  {
   if(avtx->GetZ() < fVertexZ[0]) return kFALSE;
   if(avtx->GetZ() > fVertexZ[1]) return kFALSE;
  }

 } // else if(aAOD)
 fCentralityHistVzCut->Fill(centrality);

 return kTRUE;

} // Bool_t AliAnalysisTaskNBodyFemtoscopy::CommonEventCuts(AliVEvent *ave)



//================================================================================================================
Bool_t AliAnalysisTaskNBodyFemtoscopy::CommonTrackCuts(AliAODTrack *atrack)
{
 // Apply track cuts.
 if(fApplyCommonTrackCuts) {

  TString sMethodName = "Bool_t AliAnalysisTaskNBodyFemtoscopy::CommonTrackCuts(AliAODTrack *atrack)";
  if(!atrack){Fatal(sMethodName.Data(),"!atrack");}

  if(!atrack->TestFilterBit(128)) return kFALSE; 

  Double_t pt = atrack->Pt(); // Pt
  if(atrack->Eta()<fEtaRange[0]) return kFALSE;
  if(atrack->Eta()>=fEtaRange[1]) return kFALSE;
  fPtHistEtaCut->Fill(pt); // filling pt distribution after eta cuts
  if(pt<fPtRange[0]) return kFALSE;
  if(pt>=fPtRange[1]) return kFALSE;
  fPtHistEtaCutPTCut->Fill(pt); // filling pt distribution after eta and pt cuts
  if(atrack->Phi()<fPhiRange[0]) return kFALSE;
  if(atrack->Phi()>=fPhiRange[1]) return kFALSE;
  fPtHistEtaCutPTCutPhiCut->Fill(pt); // filling pt distribution after eta, pt and phi cuts
  if(fRejectFakeTracks && atrack->GetLabel()<0) return kFALSE; // TBI is it more precise <=0 ?
  
 }
 return kTRUE;

} // Bool_t AliAnalysisTaskNBodyFemtoscopy::CommonTrackCuts(AliAODTrack *atrack)
//================================================================================================================
Bool_t AliAnalysisTaskNBodyFemtoscopy::GlobalTrackCuts(AliAODTrack *gtrack)
{
 // Check if the track passes common global track cuts (irrespectively of PID).

 TString sMethodName = "Bool_t AliAnalysisTaskNBodyFemtoscopy::GlobalTrackCuts(AliAODTrack *gtrack)";
 if(!gtrack){Fatal(sMethodName.Data(),"!gtrack");}

 if(gtrack->Pt()<fPtRange[0]) return kFALSE;
 if(gtrack->Pt()>=fPtRange[1]) return kFALSE;
 if(gtrack->Eta()<fEtaRange[0]) return kFALSE;
 if(gtrack->Eta()>=fEtaRange[1]) return kFALSE;
 if(gtrack->Phi()<fPhiRange[0]) return kFALSE;
 if(gtrack->Phi()>=fPhiRange[1]) return kFALSE;


 if(fRejectFakeTracks && gtrack->GetLabel()<0) return kFALSE; // TBI is it more precise <=0 ?

 return kTRUE;

} // Bool_t AliAnalysisTaskNBodyFemtoscopy::GlobalTrackCuts(AliAODTrack *gtrack)

//=======================================================================================================================

Bool_t AliAnalysisTaskNBodyFemtoscopy::Pion(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)
{
 // Is this primary/secondary (anti-)pion?

 // a) Insanity checks;
 // b) Trivial checks; // TBI
 // c) Track quality cuts;
 // d) PID;
 // e) PID MC (if available).

 // a) Insanity checks
 TString sMethodName = "Bool_t AliAnalysisTaskNBodyFemtoscopy::Pion(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)";
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

  if(!(dSigmaTPCPion < fInclusiveSigmaCuts[0] && dSigmaTPCProton > fExclusiveSigmaCuts[0][2] && dSigmaTPCKaon > fExclusiveSigmaCuts[0][1])) return kFALSE;
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
  if(charge < 0 && mcParticle->GetPdgCode() == -211) { fTestPIDTrueFalsePositive->Fill("Pion-", "Pion-",1); return kTRUE;}
  else if(charge > 0 && mcParticle->GetPdgCode() == 211) {fTestPIDTrueFalsePositive->Fill("Pion+","Pion+",1); return kTRUE;}
  else 
  {
   Int_t PID = mcParticle->GetPdgCode();
   if(charge<0)
   {
    if(PID==221) fTestPIDTrueFalsePositive->Fill("Pion+", "Pion-",1); 
    else if(PID==-321) fTestPIDTrueFalsePositive->Fill("Kaon-", "Pion-",1); 
    else if(PID==321) fTestPIDTrueFalsePositive->Fill("Kaon+", "Pion-",1); 
    else if(PID==-2212) fTestPIDTrueFalsePositive->Fill("Proton-", "Pion-",1); 
    else if(PID==2212) fTestPIDTrueFalsePositive->Fill("Proton+", "Pion-",1); 
    else fTestPIDTrueFalsePositive->Fill("Other", "Pion-",1); 
    return kFALSE;
   }
   else 
   {
   	if(PID==-221) fTestPIDTrueFalsePositive->Fill("Pion-", "Pion+",1); 
    else if(PID==-321) fTestPIDTrueFalsePositive->Fill("Kaon-", "Pion+",1); 
    else if(PID==321) fTestPIDTrueFalsePositive->Fill("Kaon+", "Pion+",1); 
    else if(PID==-2212) fTestPIDTrueFalsePositive->Fill("Proton-", "Pion+",1); 
    else if(PID==2212) fTestPIDTrueFalsePositive->Fill("Proton+", "Pion+",1); 
    else fTestPIDTrueFalsePositive->Fill("Other", "Pion+",1); 
    return kFALSE;
   }
  }
 } // if(fProcessBothKineAndReco)

 return kTRUE;

} // Bool_t AliAnalysisTaskNBodyFemtoscopy::Pion(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)


//=======================================================================================================================

Bool_t AliAnalysisTaskNBodyFemtoscopy::Kaon(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)
{
 // Is this primary/secondary (anti-)kaon?

 // a) Insanity checks;
 // b) Trivial checks;
 // c) Track quality cuts;
 // d) PID;
 // e) PID MC (if available).

 // a) Insanity checks:
 TString sMethodName = "Bool_t AliAnalysisTaskNBodyFemtoscopy::Kaon(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)";
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
  if(!(dSigmaTPCKaon < fInclusiveSigmaCuts[1] && dSigmaTPCProton > fExclusiveSigmaCuts[1][2] && dSigmaTPCPion > fExclusiveSigmaCuts[1][0])) return kFALSE;
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
  if(charge < 0 && mcParticle->GetPdgCode() == -321) {fTestPIDTrueFalsePositive->Fill("Kaon-","Kaon-",1); return kTRUE;}
  else if(charge > 0 && mcParticle->GetPdgCode() == 321) {fTestPIDTrueFalsePositive->Fill("Kaon+","Kaon+",1); return kTRUE;}
  else 
  {
   Int_t PID = mcParticle->GetPdgCode();
   if(charge<0)
   {
    if(PID==-221) fTestPIDTrueFalsePositive->Fill("Pion-", "Kaon-",1); 
    else if(PID==221) fTestPIDTrueFalsePositive->Fill("Pion+", "Kaon-",1); 
    else if(PID==321) fTestPIDTrueFalsePositive->Fill("Kaon+", "Kaon-",1); 
    else if(PID==-2212) fTestPIDTrueFalsePositive->Fill("Proton-", "Kaon-",1); 
    else if(PID==2212) fTestPIDTrueFalsePositive->Fill("Proton+", "Kaon-",1); 
    else fTestPIDTrueFalsePositive->Fill("Other", "Kaon-",1); 
    return kFALSE;
   }
   else 
   {
   	if(PID==-221) fTestPIDTrueFalsePositive->Fill("Pion-", "Kaon+",1); 
    else if(PID==221) fTestPIDTrueFalsePositive->Fill("Pion+", "Kaon+",1); 
    else if(PID==-321) fTestPIDTrueFalsePositive->Fill("Kaon-", "Kaon+",1); 
    else if(PID==-2212) fTestPIDTrueFalsePositive->Fill("Proton-", "Kaon+",1); 
    else if(PID==2212) fTestPIDTrueFalsePositive->Fill("Proton+", "Kaon+",1); 
    else fTestPIDTrueFalsePositive->Fill("Other", "Kaon+",1); 
    return kFALSE;
   }
  }
 } // if(fProcessBothKineAndReco)

 return kTRUE;

} // Bool_t AliAnalysisTaskNBodyFemtoscopy::Kaon(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)

//=======================================================================================================================

Bool_t AliAnalysisTaskNBodyFemtoscopy::Proton(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)
{
 // Is this primary/secondary (anti-)proton?

 // a) Insanity checks;
 // b) Trivial checks;
 // c) Track quality cuts;
 // d) PID;
 // e) PID MC (if available).

 // a) Insanity checks:
 TString sMethodName = "Bool_t AliAnalysisTaskNBodyFemtoscopy::Proton(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)";
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



  if(!(dSigmaTPCProton < fInclusiveSigmaCuts[2] && dSigmaTPCPion > fExclusiveSigmaCuts[2][0] && dSigmaTPCKaon > fExclusiveSigmaCuts[2][1])) return kFALSE;
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
  if(charge < 0 && mcParticle->GetPdgCode() == -2212)  {fTestPIDTrueFalsePositive->Fill("Proton-","Proton-",1); return kTRUE;}
  else if(charge > 0 && mcParticle->GetPdgCode() == 2212) {fTestPIDTrueFalsePositive->Fill("Proton+","Proton+",1); return kTRUE;}
  else 
  {
   Int_t PID = mcParticle->GetPdgCode();
   if(charge<0)
   {
    if(PID==-221) fTestPIDTrueFalsePositive->Fill("Pion-", "Proton-",1); 
    else if(PID==221) fTestPIDTrueFalsePositive->Fill("Pion+", "Proton-",1); 
    else if(PID==-321) fTestPIDTrueFalsePositive->Fill("Kaon-", "Proton-",1); 
    else if(PID==321) fTestPIDTrueFalsePositive->Fill("Kaon+", "Proton-",1); 
    else if(PID==2212) fTestPIDTrueFalsePositive->Fill("Proton+", "Proton-",1); 
    else fTestPIDTrueFalsePositive->Fill("Other", "Proton-",1); 
    return kFALSE;
   }
   else 
   {
   	if(PID==-221) fTestPIDTrueFalsePositive->Fill("Pion-", "Proton+",1); 
    else if(PID==221) fTestPIDTrueFalsePositive->Fill("Pion+", "Proton+",1); 
    else if(PID==321) fTestPIDTrueFalsePositive->Fill("Kaon+", "Proton+",1); 
    else if(PID==-321) fTestPIDTrueFalsePositive->Fill("Kaon-", "Proton+",1); 
    else if(PID==-2212) fTestPIDTrueFalsePositive->Fill("Proton-", "Proton+",1); 
    else fTestPIDTrueFalsePositive->Fill("Other", "Proton+",1); 
    return kFALSE;
   }
  }
 } // if(fProcessBothKineAndReco)


 return kTRUE;

} // Bool_t AliAnalysisTaskNBodyFemtoscopy::Proton(AliAODTrack *gtrack, Int_t charge, Bool_t bPrimary)

//=======================================================================================================================

void AliAnalysisTaskNBodyFemtoscopy::GlobalTracksAODTEST(AliAODEvent *aAOD, Int_t index)
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
 TString sMethodName = "void AliAnalysisTaskNBodyFemtoscopy::GlobalTracksAODTEST(AliAODEvent *aAOD, Int_t index)";

 if(!fGlobalTracksAODTEST[index]){Fatal(sMethodName.Data(),"fGlobalTracksAODTEST[%d]",index);}

 if(0 != fGlobalTracksAODTEST[index]->GetSize()){fGlobalTracksAODTEST[index]->Delete();} // yes, this method determines mapping from scratch each time

 // b) Determine the map:
 for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++)
 {
  AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack));
  if(aodTrack)
  {
   Int_t id = aodTrack->GetID();
   if(id>=0 && !aodTrack->IsGlobalConstrained()) // TBI rethink this, it seems that id>=0 is just enough, the second constraint is most likely just an overkill
   {
    fGlobalTracksAODTEST[index]->Add(id,iTrack); // "key" = id, "value" = iTrack
   } // if(id>=0 && !aodTrack->IsGlobalConstrained())
  } // if(aodTrack)
 } // for(Int_t iTrack=0;iTrack<aAOD->GetNumberOfTracks();iTrack++)

} // void AliAnalysisTaskNBodyFemtoscopy::GlobalTracksAODTEST(AliAODEvent *aAOD)

//=======================================================================================================================
void AliAnalysisTaskNBodyFemtoscopy::CreateRandomIndices(AliAODEvent *aAOD)
{
 // Fisher-Yates algorithm:
 Int_t nPrim = aAOD->GetNumberOfTracks();
 if(nPrim > 0)
 {
  fRandomIndices = new TArrayI(nPrim);
 }
 else
 {
  return;
 }

 for(Int_t i=0;i<nPrim;i++)
 {
  fRandomIndices->AddAt(i,i);
 }
 for(Int_t i=nPrim-1;i>=1;i--)
 {
  Int_t j = gRandom->Integer(i+1);
  Int_t temp = fRandomIndices->GetAt(j);
  fRandomIndices->AddAt(fRandomIndices->GetAt(i),j);
  fRandomIndices->AddAt(temp,i);
 } // end of for(Int_t i=nPrim-1;i>=1;i--) 

} // void AliAnalysisTaskNBodyFemtoscopy::CreateRandomIndices(AliAODEvent *aAOD)

//=======================================================================================================================


void AliAnalysisTaskNBodyFemtoscopy::Terminate(Option_t *)
{
 // Accessing the merged output list.

 fHistList = (TList*)GetOutputData(1);
 if(!fHistList){exit(1);}

 // Do some calculation in offline mode here:

 // ... your code for offline calculations ...

 // Update the output file with new results:
 TFile *f = new TFile("AnalysisResults.root","RECREATE");
 fHistList->Write(fHistList->GetName(),TObject::kSingleKey);

 delete f;

} // end of void AliAnalysisTaskNBodyFemtoscopy::Terminate(Option_t *)

//================================================================================================================

void AliAnalysisTaskNBodyFemtoscopy::InitializeArrays()
{
 // Initialize all data members which are arrays in this method.
  fPtRange[0]=0.2;
  fPtRange[1]=5.0;
  fEtaRange[0]=-0.8;
  fEtaRange[1]=0.8;
  fVertexZ[0]=-10.0;
  fVertexZ[1]=10.0;
  fPhiRange[0]=0.0;
  fPhiRange[1]=TMath::TwoPi();


  for(Int_t pf=0;pf<3;pf++) // PID function [0=Pion(...),1=Kaon(...),2=Proton(...)]
  {
   fInclusiveSigmaCuts[pf] = 0.;
  }

 for(Int_t pf=0;pf<3;pf++) //PID function [0=Pion(...),1=Kaon(...),2=Proton(...)]
 {
  for(Int_t pid=0;pid<3;pid++) //PID function [0=Pion(...),1=Kaon(...),2=Proton(...)]
  {
   fExclusiveSigmaCuts[pf][pid] = 0.;
  }
 }

 for(Int_t me=0;me<3;me++) // [0] is buffer for 1st event; [1] for 2nd, etc.
 {
  fGlobalTracksAODTEST[me] = NULL;
  fAllTracksTEST[me] = NULL;
 }


} // void AliAnalysisTaskNBodyFemtoscopy::InitializeArrays()

//=======================================================================================================================

void AliAnalysisTaskNBodyFemtoscopy::BookAndNestAllLists()
{
 // Book and nest all lists nested in the base list fHistList.

 // a) Book and nest lists for control histograms;
 // b) Book and nest lists for final results.

 TString sMethodName = "void AliAnalysisTaskNBodyFemtoscopy::BookAndNestAllLists()";
 if(!fHistList){Fatal(sMethodName.Data(),"fHistList is NULL");}

 // a) Book and nest lists for control histograms:
 fControlHistogramsList = new TList();
 fControlHistogramsList->SetName("ControlHistograms");
 fControlHistogramsList->SetOwner(kTRUE);
 fHistList->Add(fControlHistogramsList);

 // b) Book and nest lists for final results:
 fFinalResultsList = new TList();
 fFinalResultsList->SetName("FinalResults");
 fFinalResultsList->SetOwner(kTRUE);
 fHistList->Add(fFinalResultsList);

} // void AliAnalysisTaskNBodyFemtoscopy::BookAndNestAllLists()

//=======================================================================================================================

void AliAnalysisTaskNBodyFemtoscopy::BookControlHistograms()
{
 // Book all control histograms.

 // a) Book histogram to hold pt spectra;
 // b) Book histogram for centrality distribution
 // c) Book histogram for track number distribution
 // d) Book histogram for centrality distribution after event cuts
 // e) Book histogram for pt spectra after track cuts

 // a) Book histogram to hold pt spectra before cuts:
 fPtHist = new TH1F("fPtHist","atrack->Pt()",fNbinsPt,fMinBinPt,fMaxBinPt);
 fPtHist->SetFillColor(kBlue-10);
 fPtHist->GetXaxis()->SetTitle("p_{t}");
 fControlHistogramsList->Add(fPtHist);
 
 // b) Book histogram for centrality distribution before cuts:
 fCentralityHist = new TH1F("fCentralityHist","ams->GetMultiplicityPercentile()",fNbinsCentrality,fMinCentrality,fMaxCentrality);
 fCentralityHist->SetFillColor(kRed-10);
 fCentralityHist->GetXaxis()->SetTitle("Centrality percentile");
 //fCentralityHist->SetStats(kTRUE);
 fControlHistogramsList->Add(fCentralityHist);

 // c) Book histogram for track number distribution
 fNumberOfTracksHist = new TH1F("fNumberOfTracksHist","aAOD->GetNumberOfTracks()",fNbinsMultiplicity,fMinBinMultiplicity,fMaxBinMultiplicity);
 fNumberOfTracksHist->SetFillColor(kRed-10);
 fNumberOfTracksHist->GetXaxis()->SetTitle("Number of tracks");
 fControlHistogramsList->Add(fNumberOfTracksHist);

 // c) Book histogram for centrality distribution after event cuts
 fCentralityHistVzCut = new TH1F("fCentralityHistVzCut","ams->GetMultiplicityPercentile() after vz cut",fNbinsCentrality,fMinCentrality,fMaxCentrality);
 fCentralityHistVzCut->SetFillColor(kRed-10);
 fCentralityHistVzCut->GetXaxis()->SetTitle("Centrality percentile");
 fControlHistogramsList->Add(fCentralityHistVzCut);

 // d) Book histogram for pt spectra after track cuts
 fPtHistEtaCut = new TH1F("fPtHistEtaCut","atrack->Pt() after Eta cuts",fNbinsPt,fMinBinPt,fMaxBinPt);
 fPtHistEtaCut->SetFillColor(kBlue-10);
 fPtHistEtaCut->GetXaxis()->SetTitle("p_{t}");
 fControlHistogramsList->Add(fPtHistEtaCut);

 fPtHistEtaCutPTCut = new TH1F("fPtHistEtaCutPTCut","atrack->Pt() after Eta and pT cuts",fNbinsPt,fMinBinPt,fMaxBinPt);
 fPtHistEtaCutPTCut->SetFillColor(kBlue-10);
 fPtHistEtaCutPTCut->GetXaxis()->SetTitle("p_{t}");
 fControlHistogramsList->Add(fPtHistEtaCutPTCut);

 fPtHistEtaCutPTCutPhiCut = new TH1F("fPtHistEtaCutPTCutPhiCut","atrack->Pt() after Eta, pT and Phi cuts",fNbinsPt,fMinBinPt,fMaxBinPt);
 fPtHistEtaCutPTCutPhiCut->SetFillColor(kBlue-10);
 fPtHistEtaCutPTCutPhiCut->GetXaxis()->SetTitle("p_{t}");
 fControlHistogramsList->Add(fPtHistEtaCutPTCutPhiCut);

 fNumberOfTracksHistAfterAllCuts = new TH1F("fNumberOfTracksHistAfterAllCuts","Tracks that passed selection",fNbinsMultiplicity,fMinBinMultiplicity,fMaxBinMultiplicity);
 fNumberOfTracksHistAfterAllCuts->SetFillColor(kRed-10);
 fNumberOfTracksHistAfterAllCuts->GetXaxis()->SetTitle("Number of tracks");
 fControlHistogramsList->Add(fNumberOfTracksHistAfterAllCuts);
 const Int_t xBinsTemp = 7;
 const Int_t yBinsTemp = 6;
 fTestPIDTrueFalsePositive = new TH2F("fTestPIDTrueFalsePositive","True positives PID",xBinsTemp,0,6,yBinsTemp,0,5);
 const char *particleNameY[yBinsTemp] = {"Pion+","Pion-","Proton+","Proton-","Kaon+","Kaon-"};
 for (Int_t i=1;i<=yBinsTemp;i++) fTestPIDTrueFalsePositive->GetYaxis()->SetBinLabel(i,particleNameY[i-1]);
 const char *particleNameX[xBinsTemp] = {"Pion+","Pion-","Proton+","Proton-","Kaon+","Kaon-", "Other"};
 for (Int_t i=1;i<=xBinsTemp;i++) fTestPIDTrueFalsePositive->GetXaxis()->SetBinLabel(i,particleNameX[i-1]);
 fControlHistogramsList->Add(fTestPIDTrueFalsePositive);


 fTestTPCOnlyVsGlobal = new TH1F("fTestTPCOnlyVsGlobal","pt(tpc)-pt(global)",fNbinsTPCOnlyVsGlobal,fMinBinTPCOnlyVsGlobal,fMaxBinTPCOnlyVsGlobal);
 fTestTPCOnlyVsGlobal->SetFillColor(kRed-10);
 fTestTPCOnlyVsGlobal->GetXaxis()->SetTitle("Number of tracks");
 fControlHistogramsList->Add(fTestTPCOnlyVsGlobal);
} // void AliAnalysisTaskNBodyFemtoscopy::BookControlHistograms()

//=======================================================================================================================

void AliAnalysisTaskNBodyFemtoscopy::BookFinalResultsHistograms()
{
 // Book all histograms to hold the final results.

} // void AliAnalysisTaskNBodyFemtoscopy::BookFinalResultsHistograms()

//=======================================================================================================================

//=======================================================================================================================

void AliAnalysisTaskNBodyFemtoscopy::BookEverything()
{

 // General
 fPIDResponse = new AliPIDResponse();

 // PID sigma cuts: default
 if(fUseDefaultInclusiveSigmaCuts)
 {
  fInclusiveSigmaCuts[0] = 3.; // i.e. in function Pion(...) the inclusive cut for pions is 2 sigma
  fInclusiveSigmaCuts[1] = 3.; // i.e. in function Kaon(...) the inclusive cut for kaons is 2 sigma
  fInclusiveSigmaCuts[2] = 3.; // i.e. in function Proton(...) the inclusive cut for protons is 2 sigma
 }
 if(fUseDefaultExclusiveSigmaCuts)
 {
   // Pion(...)
  fExclusiveSigmaCuts[0][1] = 4.; // i.e. in function Pion(...) the exclusive cut for kaons is 4 sigma
  fExclusiveSigmaCuts[0][2] = 4.; // i.e. in function Pion(...) the exclusive cut for protons is 4 sigma
  // Kaon(...)
  fExclusiveSigmaCuts[1][0] = 4.; // i.e. in function Kaon(...) the exclusive cut for pions is 4 sigma
  fExclusiveSigmaCuts[1][2] = 4.; // i.e. in function Kaon(...) the exclusive cut for protons is 4 sigma
  // Proton(...)
  fExclusiveSigmaCuts[2][0] = 4.; // i.e. in function Proton(...) the exclusive cut for pions is 4 sigma
  fExclusiveSigmaCuts[2][1] = 4.; // i.e. in function Proton(...) the exclusive cut for kaons is 4 sigma
 } // if(fUseDefaultExclusiveSigmaCuts)

 // For global track map
 for(Int_t me=0;me<3;me++) // [0] is buffer for 1st event; [1] for 2nd, etc.
 {
  fAllTracksTEST[me] = new TClonesArray("AliAODTrack",10000);
  fGlobalTracksAODTEST[me] = new TExMap();
 }

} // void AliAnalysisTaskNBodyFemtoscopy::BookEverything()

//====