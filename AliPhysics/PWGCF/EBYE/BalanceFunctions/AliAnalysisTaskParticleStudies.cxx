#include "TChain.h"
#include "TH2D.h"
#include "TH3D.h"


#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliLog.h"
#include "AliAnalysisUtils.h"


#include "AliAnalysisTaskSE.h"
#include "AliCollisionGeometry.h"
#include "AliGenEventHeader.h"
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"



#include "AliAnalysisTaskParticleStudies.h"

// Analysis task for an analysis on rho mesons and other particles
// mainly to understand if BF behaviour could come from boosted rho mesons
// Authors: m.weber@cern.ch

ClassImp(AliAnalysisTaskParticleStudies)

//________________________________________________________________________
AliAnalysisTaskParticleStudies::AliAnalysisTaskParticleStudies(const char *name) 
  : AliAnalysisTaskSE(name),
  fListQA(0x0),
  fPdgCode(113),
  fMotherPdgCode(0),
  fPtMin(0),
  fPtMax(1000),
  fEtaMin(-10),
  fEtaMax(10),
  fHistTrackStats(0),
  fHistIPPt(0),
  fHistEtaPhiPt(0)
{

  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TList::Class());

}

//________________________________________________________________________
AliAnalysisTaskParticleStudies::~AliAnalysisTaskParticleStudies() {
  // Destructor
  // ... not implemented
}

//________________________________________________________________________
void AliAnalysisTaskParticleStudies::UserCreateOutputObjects() {
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
  fHistTrackStats = new TH2D("fHistTrackStats","Track statistics; particles; impact parameter (fm);N_{events}",500,0,200,500,0,20);
  fListQA->Add(fHistTrackStats);

  fHistIPPt = new TH2D("fHistIPPt","track distribution; impact parameter (fm);p_{T} (GeV/c)",200,0,20,200,0,20);
  fListQA->Add(fHistIPPt);

  fHistEtaPhiPt = new TH3D("fHistEtaPhiPt","track distribution;#eta;#varphi (rad);p_{T} (GeV/c)",100,-1.5,1.5,100,0,TMath::Pi(),100,0,10);
  fListQA->Add(fHistEtaPhiPt);

  // Post output data.
  PostData(1, fListQA);

  AliInfo("Finished setting up the Output");
  TH1::AddDirectory(oldStatus);
}

//________________________________________________________________________
void AliAnalysisTaskParticleStudies::UserExec(Option_t *) {
  // Main loop
  // Called for each event

  AliVEvent* event = MCEvent();     
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
  Int_t acceptedTracks = GetAcceptedTracks(event,lMultiplicityVar);

  fHistTrackStats->Fill(acceptedTracks,lMultiplicityVar);

}

//________________________________________________________________________
void  AliAnalysisTaskParticleStudies::FinishTaskOutput(){
  // Finish task output
  // not implemented ...

}

//________________________________________________________________________
void AliAnalysisTaskParticleStudies::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query
  // not implemented ...

}


//________________________________________________________________________
Double_t AliAnalysisTaskParticleStudies::IsEventAccepted(AliVEvent *event){
  // Checks the Event cuts

  Double_t centralityEvent = -1.;

  AliCollisionGeometry* headerH = dynamic_cast<AliCollisionGeometry*>(dynamic_cast<AliMCEvent*>(event)->GenEventHeader());      
  if(headerH){
    centralityEvent = headerH->ImpactParameter();
    return centralityEvent;
  }//MC header

  // Double_t centralityEvent = header->GetCentralityP()->GetCentralityPercentile(sCentralityEstimator[0].Data());
  // if(centralityEvent > fCentralityPercentileMin && centralityEvent < fCentralityPercentileMax ){
  //   return centralityEvent;
  // }

  // in all other cases return -1 (event not accepted)
  return -1.;
}

//________________________________________________________________________
Int_t AliAnalysisTaskParticleStudies::GetAcceptedTracks(AliVEvent *mcEvent, Double_t vIP){
  // Checks track cuts (filter bits)
  // Fills QA histograms

  Int_t acceptedTracks = 0;

  Double_t vEta;
  Double_t vPhi;
  Double_t vPt;

  
  // Loop over MC tracks in event
  for (Int_t iTracks = 0; iTracks < mcEvent->GetNumberOfTracks(); iTracks++) {
    AliMCParticle* mcTrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(iTracks));
    if (!mcTrack) {
      AliError(Form("Could not receive MC track %d", iTracks));
      continue;
    }

    // remove all non PDG particles
    if( TMath::Abs(mcTrack->PdgCode()) !=  fPdgCode ) 
      continue;

    // if fMotherPdgCode != 0
    // remove all particles that do not come from mother PDG
    if( fMotherPdgCode != 0 ) {
      Bool_t removeParticle = kTRUE;
      
      Int_t motherIndex = mcTrack->GetMother();
      if(motherIndex != -1) {
	AliMCParticle* motherParticle = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(motherIndex));
	if(motherParticle) {
	  Int_t pdgCodeOfMother = motherParticle->PdgCode();

	  if( TMath::Abs(pdgCodeOfMother) ==  fMotherPdgCode ) 
	    removeParticle = kFALSE;
	}
      }
      if(removeParticle)
	continue;
    }
      
    // some track parameters
    vEta    = mcTrack->Eta();
    vPhi    = mcTrack->Phi();
    vPt     = mcTrack->Pt();
    
    
    // kinematic cuts
    if( vPt > fPtMax || vPt < fPtMin )
      continue;
    if( vEta > fEtaMax || vEta < fEtaMin )
      continue;
    
    // fill QA
    fHistEtaPhiPt->Fill(vEta,vPhi,vPt);
    fHistIPPt->Fill(vIP,vPt);
    
    // count tracks
    acceptedTracks++;
    
  }//track loop
  
  return acceptedTracks;  
}
