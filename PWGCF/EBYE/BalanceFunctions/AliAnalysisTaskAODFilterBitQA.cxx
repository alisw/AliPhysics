#include "TChain.h"
#include "TH2D.h"
#include "TH3D.h"


#include "AliAODEvent.h"
#include "AliAODRecoDecay.h"
#include "AliAODTrack.h"
#include "AliLog.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"
#include "AliInputEventHandler.h"
#include "AliMultSelection.h"



#include "AliAnalysisTaskAODFilterBitQA.h"

// Analysis task for the QA of AOD track filter bits
// Authors: m.weber@cern.ch

ClassImp(AliAnalysisTaskAODFilterBitQA)

//________________________________________________________________________
AliAnalysisTaskAODFilterBitQA::AliAnalysisTaskAODFilterBitQA(const char *name) 
  : AliAnalysisTaskSE(name),
  fArrayMC(0x0),
  fListQA(0x0),
  useCentrality(kFALSE),
  fUseMultSelectionFramework(kFALSE),
  fUseUncheckedCentrality(kFALSE),
  fillOnlySecondaries(kFALSE),
  fillHFVertexingTracks(kFALSE),
  fHFBranchName("D0toKpi"),
  fBitIgnore1(-1),
  fBitIgnore2(-1),
  fCentralityPercentileMin(0.),
  fCentralityPercentileMax(80.), 
  fPtMin(0),
  fPtMax(1000),
  fEtaMin(-10),
  fEtaMax(10),
  fHistEventStats(0),
  fHistTrackStats(0)
{

  for(Int_t iCharge = 0; iCharge < gNCharge; iCharge++){
    for(Int_t iTrackBit = 0; iTrackBit < gBitMax; iTrackBit++){
      fHistKinematics[iCharge][iTrackBit] = NULL;
      fHistDCAconstrained[iCharge][iTrackBit] = NULL;
      fHistDCAglobal[iCharge][iTrackBit]  = NULL;
      fHistChiClus[iCharge][iTrackBit]    = NULL;
    }
  }
  
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TList::Class());

}

//________________________________________________________________________
AliAnalysisTaskAODFilterBitQA::~AliAnalysisTaskAODFilterBitQA() {
  // Destructor
  // ... not implemented
}

//________________________________________________________________________
void AliAnalysisTaskAODFilterBitQA::UserCreateOutputObjects() {
  // Create histograms
  // Called once

   // global switch disabling the reference 
  // (to avoid "Replacing existing TH1" if several wagons are created in train)
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  // QA list
  fListQA = new TList();
  fListQA->SetName("listQA");
  fListQA->SetOwner();

  // QA histograms
  fHistEventStats = new TH2D("fHistEventStats","Event statistics;Centrality;EventTriggerBit;N_{events}",104,-2,102,32,0,32);
  fListQA->Add(fHistEventStats);
  fHistTrackStats = new TH2D("fHistTrackStats","Track statistics;Centrality;TrackFilterBit;N_{events}",100,0,100,gBitMax,0,gBitMax);
  fListQA->Add(fHistTrackStats);

  TString sCharge[gNCharge] = {"Plus","Minus"};
  
  for(Int_t iCharge = 0; iCharge < gNCharge; iCharge++){
    for(Int_t iTrackBit = 0; iTrackBit < gBitMax; iTrackBit++){
      fHistKinematics[iCharge][iTrackBit] = new TH3D(Form("Bit%d_%s_Kinematics",iTrackBit,sCharge[iCharge].Data()),Form("Bit%d_%s_Kinematics;#eta;#varphi (rad);p_{T} (GeV/c)",iTrackBit,sCharge[iCharge].Data()),100,-1.0,1.0,100,0,TMath::Pi()*2,100,0,10);
      fHistDCAconstrained[iCharge][iTrackBit] = new TH2D(Form("Bit%d_%s_DCAconstrained",iTrackBit,sCharge[iCharge].Data()),Form("Bit%d_%s_DCAconstrained;DCA XY [Constrained] (cm);DCA Z [Constrained] (cm)",iTrackBit,sCharge[iCharge].Data()),100,-5.0,5.0,100,-5.0,5.0);
      fHistDCAglobal[iCharge][iTrackBit]  = new TH3D(Form("Bit%d_%s_DCAglobal",iTrackBit,sCharge[iCharge].Data()),Form("Bit%d_%s_DCAglobal;DCA X [Global] (cm);DCA Y [Global] (cm);DCA Z [Global] (cm)",iTrackBit,sCharge[iCharge].Data()),100,-5.0,5.0,100,-5.0,5.0,100,-5.0,5.0);
fHistChiClus[iCharge][iTrackBit]    = new TH2D(Form("Bit%d_%s_ChiClus",iTrackBit,sCharge[iCharge].Data()),Form("Bit%d_%s_ChiClus;#chi^{2} [Fit];N_{clus} [TPC]",iTrackBit,sCharge[iCharge].Data()),100,-1.0,5.0,160,0,160.0);
      fListQA->Add(fHistKinematics[iCharge][iTrackBit]);
      fListQA->Add(fHistDCAconstrained[iCharge][iTrackBit]);
      fListQA->Add(fHistDCAglobal[iCharge][iTrackBit]);
      fListQA->Add(fHistChiClus[iCharge][iTrackBit]);
    }
  }

  // Post output data.
  PostData(1, fListQA);

  AliInfo("Finished setting up the Output");
  TH1::AddDirectory(oldStatus);
}

//________________________________________________________________________
void AliAnalysisTaskAODFilterBitQA::UserExec(Option_t *) {
  // Main loop
  // Called for each event

  AliVEvent* event = dynamic_cast<AliVEvent*>(InputEvent());     
  if(!event) {
    AliError("event not available");
    return;
  }

  // MC information (set if available)
  fArrayMC = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));   
    
  // check event cuts
  Double_t lMultiplicityVar = -1;
  lMultiplicityVar = IsEventAccepted(event);

  // event QA
  Bool_t isSelectedMB = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
  Bool_t isSelectedMu = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMUON);

  for(Int_t iTriggerBit = 0; iTriggerBit < 32; iTriggerBit++){
    Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & (1<<iTriggerBit));
    
    fHistEventStats->Fill(lMultiplicityVar,iTriggerBit,isSelected);
  }
   
  if(lMultiplicityVar < 0){ 
    return;
  }

  // fill HF vertexing tracks
  if(fillHFVertexingTracks){
    GetAcceptedHFVertexingTracks(event,lMultiplicityVar);
  }
  else{
    // get the accepted tracks in main event
    GetAcceptedTracks(event,lMultiplicityVar);
  }
}

//________________________________________________________________________
void  AliAnalysisTaskAODFilterBitQA::FinishTaskOutput(){
  // Finish task output
  // not implemented ...

}

//________________________________________________________________________
void AliAnalysisTaskAODFilterBitQA::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query
  // not implemented ...

}


//________________________________________________________________________
Double_t AliAnalysisTaskAODFilterBitQA::IsEventAccepted(AliVEvent *event){
  // Checks the Event cuts

  // still hard coded
  Double_t fVxMax = 0.5;
  Double_t fVyMax = 0.5;
  Double_t fVzMax = 10.0;
  TString fCentralityEstimator = "V0M";

  Float_t gCentrality      = -1.;
  const AliVVertex *vertex = event->GetPrimaryVertex();
  
  if(vertex) {
    Double32_t fCov[6];
    vertex->GetCovarianceMatrix(fCov);
    if(vertex->GetNContributors() > 0) {
      if(fCov[5] != 0) {
	if(TMath::Abs(vertex->GetX()) < fVxMax) {
	  if(TMath::Abs(vertex->GetY()) < fVyMax) {
	    if(TMath::Abs(vertex->GetZ()) < fVzMax) {

	      // get the reference multiplicty or centrality (if required)
	      if(useCentrality){

		// use AliMultSelection framework
		if (fUseMultSelectionFramework) {

		  AliMultSelection *multSelection = (AliMultSelection*) event->FindListObject("MultSelection");
		  if (!multSelection)
		    AliFatal("MultSelection not found in input event");
		  
		  if (fUseUncheckedCentrality)
		    gCentrality = multSelection->GetMultiplicityPercentile(fCentralityEstimator, kFALSE);
		  else
		    gCentrality = multSelection->GetMultiplicityPercentile(fCentralityEstimator, kTRUE);

		  // error handling
		  if (gCentrality > 100)
		    gCentrality = -1;
		}
		else{
		  AliAODHeader *header = (AliAODHeader*) event->GetHeader();
		  gCentrality = header->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());
		}

		if((gCentrality > fCentralityPercentileMin) && (gCentrality < fCentralityPercentileMax)){		  
		  return gCentrality;
		}
		else{
		  return -1;
		}//centrality range
	      }//use centrality

	      // if not using centrality/multiplicty, return 1
	      return 1;

	    }//Vz cut
	  }//Vy cut
	}//Vx cut
      }//proper vertex resolution
    }//proper number of contributors
  }//vertex object valid
  
  // in all other cases return -1 (event not accepted)
  return -1;
}

//________________________________________________________________________
void AliAnalysisTaskAODFilterBitQA::GetAcceptedTracks(AliVEvent *event, Double_t gCentrality){
  // Checks track cuts (filter bits)
  // Fills QA histograms


  Double_t vEta;
  Double_t vPhi;
  Double_t vPt;
  Double_t vDCAconstrainedxy;
  Double_t vDCAconstrainedz; 
  Double_t vDCAglobalx;
  Double_t vDCAglobaly;
  Double_t vDCAglobalz;
  Double_t vChi2;
  Double_t vClus;

  Double_t pos[3];
  Double_t v[3];

  const AliVVertex *vertex = event->GetPrimaryVertex();
  vertex->GetXYZ(v);

  // Loop over tracks in event
  for (Int_t iTracks = 0; iTracks < event->GetNumberOfTracks(); iTracks++) {
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack *>(event->GetTrack(iTracks));
    if (!aodTrack) {
      AliError(Form("Could not receive track %d", iTracks));
      continue;
    }

    // get MC information (if available)
    if(fArrayMC && fillOnlySecondaries){
      
      Int_t label = aodTrack->GetLabel();
      AliAODMCParticle *mcTrack = (AliAODMCParticle *)fArrayMC->At(TMath::Abs(label));

      if(mcTrack->IsPhysicalPrimary())
	continue;      
    }

    // track parameters
    vEta    = aodTrack->Eta();
    vPhi    = aodTrack->Phi();// * TMath::RadToDeg();
    vPt     = aodTrack->Pt();
    vDCAconstrainedxy  = aodTrack->DCA();     
    vDCAconstrainedz   = aodTrack->ZAtDCA();   
    vChi2   = aodTrack->Chi2perNDF(); 
    vClus   = aodTrack->GetTPCNcls(); 

    // kinematic cuts
    if( vPt > fPtMax || vPt < fPtMin )
      continue;
    if( vEta > fEtaMax || vEta < fEtaMin )
      continue;

    // if not constrained track the position is stored (primary vertex to be subtracted)
    aodTrack->GetXYZ(pos);
    vDCAglobalx  = pos[0] - v[0];
    vDCAglobaly  = pos[1] - v[1];
    vDCAglobalz  = pos[2] - v[2];
    
    // fill for separately for positive and negative charges
    Int_t iCharge = -1;
    // positive 
    if(aodTrack->Charge() > 0)
      iCharge = 0;
    else if(aodTrack->Charge() < 0)
      iCharge = 1;
    else{
      AliError("Charge==0?");
      iCharge = -1;
    }


    // AOD track cuts
    if(iCharge > -1){
      for(Int_t iTrackBit = 0; iTrackBit < gBitMax-1; iTrackBit++){
	fHistTrackStats->Fill(gCentrality,iTrackBit,aodTrack->TestFilterBit(1<<iTrackBit));
	
	if(aodTrack->TestFilterBit(1<<iTrackBit)){
	  fHistKinematics[iCharge][iTrackBit]->Fill(vEta,vPhi,vPt);
	  fHistDCAconstrained[iCharge][iTrackBit]->Fill(vDCAconstrainedxy,vDCAconstrainedz);
	  fHistDCAglobal[iCharge][iTrackBit]->Fill(vDCAglobalx,vDCAglobaly,vDCAglobalz);
	  fHistChiClus[iCharge][iTrackBit]->Fill(vChi2,vClus);
	} 
      }//bit loop
      
      // fill all tracks in last bit
      fHistTrackStats->Fill(gCentrality,gBitMax-1,1);
      fHistKinematics[iCharge][gBitMax-1]->Fill(vEta,vPhi,vPt);
      fHistDCAconstrained[iCharge][gBitMax-1]->Fill(vDCAconstrainedxy,vDCAconstrainedz);
      fHistDCAglobal[iCharge][gBitMax-1]->Fill(vDCAglobalx,vDCAglobaly,vDCAglobalz);
      fHistChiClus[iCharge][gBitMax-1]->Fill(vChi2,vClus);
      
    }//charge positive or negative
  }//track loop
  
  return;  
}


//________________________________________________________________________
void AliAnalysisTaskAODFilterBitQA::GetAcceptedHFVertexingTracks(AliVEvent *event, Double_t gCentrality){
  // Checks track cuts (filter bits)
  // from daughters of HF candidates
  // Fills QA histograms

  Double_t vEta;
  Double_t vPhi;
  Double_t vPt;
  Double_t vDCAconstrainedxy;
  Double_t vDCAconstrainedz; 
  Double_t vDCAglobalx;
  Double_t vDCAglobaly;
  Double_t vDCAglobalz;
  Double_t vChi2;
  Double_t vClus;

  Double_t pos[3];
  Double_t v[3];

  Int_t IDList[1000];
  Int_t IDListLength = 0;
  for(Int_t i = 0; i < 1000; i++){
    IDList[i] = -5;
  }

  const AliVVertex *vertex = event->GetPrimaryVertex();
  vertex->GetXYZ(v);

  // =================================================================================
  // HF part (taken from AliAnalysisTaskSEDmesonsFilterCJ)

  TClonesArray *arrayDStartoD0pi = (TClonesArray*)event->GetList()->FindObject(fHFBranchName.Data());
   
  if (!arrayDStartoD0pi) {
    AliInfo(Form("Could not find array %s, skipping the event",fHFBranchName.Data()));
    return;
  } else {
    AliDebug(2, Form("Found %d vertices",arrayDStartoD0pi->GetEntriesFast()));   
  }

  // loop over tracks
   const Int_t NVertices = arrayDStartoD0pi->GetEntriesFast();

   for (Int_t iVertex = 0; iVertex < NVertices; ++iVertex) {
     
     AliAODRecoDecay *HFvertex = static_cast<AliAODRecoDecay*>(arrayDStartoD0pi->At(iVertex));
   
     // Loop over tracks (daughters of D candidates)
     for (Int_t iProng = 0; iProng<HFvertex->GetNProngs(); iProng++) { 
       
       AliAODTrack *daughter = (AliAODTrack*)HFvertex->GetDaughter(iProng);
       if (!daughter){
	 AliError(Form("Could not receive track %d %d", iVertex, iProng));
	 continue;
       }


       // get MC information (if available)
       if(fArrayMC && fillOnlySecondaries){
	 
	 Int_t label = daughter->GetLabel();
	 AliAODMCParticle *mcTrack = (AliAODMCParticle *)fArrayMC->At(TMath::Abs(label));
	 
	 if(mcTrack->IsPhysicalPrimary())
	   continue;      
       }

       // track parameters
       vEta    = daughter->Eta();
       vPhi    = daughter->Phi();// * TMath::RadToDeg();
       vPt     = daughter->Pt();
       vDCAconstrainedxy  = daughter->DCA();     
       vDCAconstrainedz   = daughter->ZAtDCA();   
       vChi2   = daughter->Chi2perNDF(); 
       vClus   = daughter->GetTPCNcls(); 
       
       // kinematic cuts
       if( vPt > fPtMax || vPt < fPtMin )
	 continue;
       if( vEta > fEtaMax || vEta < fEtaMin )
	 continue;

       // avoid double counting (can be optimized)
       Bool_t doubleCount = kFALSE;
       Int_t  daughterID  = daughter->GetID();
       for(Int_t idx = 0; idx < IDListLength; idx++){
	 if(IDList[idx]==daughterID){
	   doubleCount = kTRUE;
	   break;
	 }
       }
       if(!doubleCount){
	 IDList[IDListLength] = daughterID;
	 IDListLength++;
       }
       else{
	 continue;
       }

       
       // if not constrained track the position is stored (primary vertex to be subtracted)
       daughter->GetXYZ(pos);
       vDCAglobalx  = pos[0] - v[0];
       vDCAglobaly  = pos[1] - v[1];
       vDCAglobalz  = pos[2] - v[2];
       
       // fill for separately for positive and negative charges
       Int_t iCharge = -1;
       // positive 
       if(daughter->Charge() > 0)
	 iCharge = 0;
       else if(daughter->Charge() < 0)
	 iCharge = 1;
       else{
	 AliError("Charge==0?");
	 iCharge = -1;
       }
       
       
       // AOD track cuts
       if(iCharge > -1){

	 // if some filter bits should be ignored, skip them here
	 if(fBitIgnore1 > -1 && daughter->TestFilterBit(1<<fBitIgnore1))
	   continue;
	 if(fBitIgnore2 > -1 && daughter->TestFilterBit(1<<fBitIgnore2))
	   continue;

	 for(Int_t iTrackBit = 0; iTrackBit < gBitMax; iTrackBit++){
	   fHistTrackStats->Fill(gCentrality,iTrackBit,daughter->TestFilterBit(1<<iTrackBit));
	   
	   if(daughter->TestFilterBit(1<<iTrackBit)){
	     fHistKinematics[iCharge][iTrackBit]->Fill(vEta,vPhi,vPt);
	     fHistDCAconstrained[iCharge][iTrackBit]->Fill(vDCAconstrainedxy,vDCAconstrainedz);
	     fHistDCAglobal[iCharge][iTrackBit]->Fill(vDCAglobalx,vDCAglobaly,vDCAglobalz);
	     fHistChiClus[iCharge][iTrackBit]->Fill(vChi2,vClus);
	   } 
	 }//bit loop

	 // fill all tracks in last bit
	 fHistTrackStats->Fill(gCentrality,gBitMax-1,1);
	 fHistKinematics[iCharge][gBitMax-1]->Fill(vEta,vPhi,vPt);
	 fHistDCAconstrained[iCharge][gBitMax-1]->Fill(vDCAconstrainedxy,vDCAconstrainedz);
	 fHistDCAglobal[iCharge][gBitMax-1]->Fill(vDCAglobalx,vDCAglobaly,vDCAglobalz);
	 fHistChiClus[iCharge][gBitMax-1]->Fill(vChi2,vClus);
	 
       }//charge positive or negative
       
     }//prong loop
   }//HF vertex loop
   
  return;  
}
