#include "TChain.h"
#include "TH2D.h"
#include "TH3D.h"


#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliLog.h"


#include "AliAnalysisTaskAODFilterBitQA.h"

// Analysis task for the QA of AOD track filter bits
// Authors: m.weber@cern.ch

ClassImp(AliAnalysisTaskAODFilterBitQA)

//________________________________________________________________________
AliAnalysisTaskAODFilterBitQA::AliAnalysisTaskAODFilterBitQA(const char *name) 
  : AliAnalysisTaskSE(name),
  fHistTrackStats(0)
{

  for(Int_t iTrackBit = 0; iTrackBit < gBitMax; iTrackBit++){
    fHistKinematics[iTrackBit] = NULL;
    fHistDCAconstrained[iTrackBit] = NULL;
    fHistDCAglobal[iTrackBit]  = NULL;
    fHistChiClus[iTrackBit]    = NULL;
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
  fHistTrackStats = new TH2D("fHistTrackStats","Track statistics;Centrality;TrackFilterBit;N_{events}",100,0,100,gBitMax,0,gBitMax);
  fListQA->Add(fHistTrackStats);

  for(Int_t iTrackBit = 0; iTrackBit < gBitMax; iTrackBit++){
    fHistKinematics[iTrackBit] = new TH3D(Form("Bit%d_Kinematics",iTrackBit),Form("Bit%d_Kinematics;#eta;#varphi (rad);p_{T} (GeV/c)",iTrackBit),100,-1.0,1.0,100,0,TMath::Pi()*2,100,0,10);
    fHistDCAconstrained[iTrackBit] = new TH2D(Form("Bit%d_DCAconstrained",iTrackBit),Form("Bit%d_DCAconstrained;DCA XY [Constrained] (cm);DCA Z [Constrained] (cm)",iTrackBit),100,-5.0,5.0,100,-5.0,5.0);
    fHistDCAglobal[iTrackBit]  = new TH3D(Form("Bit%d_DCAglobal",iTrackBit),Form("Bit%d_DCAglobal;DCA X [Global] (cm);DCA Y [Global] (cm);DCA Z [Global] (cm)",iTrackBit),100,-5.0,5.0,100,-5.0,5.0,100,-5.0,5.0);
    fHistChiClus[iTrackBit]    = new TH2D(Form("Bit%d_ChiClus",iTrackBit),Form("Bit%d_ChiClus;#chi^{2} [Fit];N_{clus} [TPC]",iTrackBit),100,-1.0,5.0,160,0,160.0);
    fListQA->Add(fHistKinematics[iTrackBit]);
    fListQA->Add(fHistDCAconstrained[iTrackBit]);
    fListQA->Add(fHistDCAglobal[iTrackBit]);
    fListQA->Add(fHistChiClus[iTrackBit]);
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

  

  // check event cuts
  Double_t lMultiplicityVar = -1;
  if((lMultiplicityVar = IsEventAccepted(event)) < 0){ 
    return;
  }

  
  // get the accepted tracks in main event
  GetAcceptedTracks(event,lMultiplicityVar);

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
  Double_t fCentralityPercentileMin = 0.;
  Double_t fCentralityPercentileMax = 80.;  
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
	      
	      // get the reference multiplicty or centrality
	      AliAODHeader *header = (AliAODHeader*) event->GetHeader();
	      gCentrality = header->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());
	      
	      if((gCentrality > fCentralityPercentileMin) && (gCentrality < fCentralityPercentileMax)){
		
		return gCentrality;
		
	      }//centrality range
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


  Short_t  vCharge;
  Double_t vEta;
  Double_t vY;
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

    // track parameters
    vCharge = aodTrack->Charge();
    vEta    = aodTrack->Eta();
    vY      = aodTrack->Y();
    vPhi    = aodTrack->Phi();// * TMath::RadToDeg();
    vPt     = aodTrack->Pt();
    vDCAconstrainedxy  = aodTrack->DCA();     
    vDCAconstrainedz   = aodTrack->ZAtDCA();   
    vChi2   = aodTrack->Chi2perNDF(); 
    vClus   = aodTrack->GetTPCNcls(); 

    // if not constrained track the position is stored (primary vertex to be subtracted)
    aodTrack->GetXYZ(pos);
    vDCAglobalx  = pos[0] - v[0];
    vDCAglobaly  = pos[1] - v[1];
    vDCAglobalz  = pos[2] - v[2];
    

    // AOD track cuts
    for(Int_t iTrackBit = 0; iTrackBit < gBitMax; iTrackBit++){
      fHistTrackStats->Fill(gCentrality,iTrackBit,aodTrack->TestFilterBit(1<<iTrackBit));
      
      if(aodTrack->TestFilterBit(1<<iTrackBit)){
	fHistKinematics[iTrackBit]->Fill(vEta,vPhi,vPt);
	fHistDCAconstrained[iTrackBit]->Fill(vDCAconstrainedxy,vDCAconstrainedz);
	fHistDCAglobal[iTrackBit]->Fill(vDCAglobalx,vDCAglobaly,vDCAglobalz);
	fHistChiClus[iTrackBit]->Fill(vChi2,vClus);
      }
      
    }//bit loop
  }//track loop

  return;  
}
