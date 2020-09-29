#include "TChain.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliMCParticle.h"
#include "AliMCEventHandler.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliVVertex.h"

#include "AliAnalysisTaskLcToD0MC.h"


// Analysis task for the BF code
// Authors: Panos Cristakoglou@cern.ch

ClassImp(AliAnalysisTaskLcToD0MC)

//________________________________________________________________________
AliAnalysisTaskLcToD0MC::AliAnalysisTaskLcToD0MC(const char *name) 
: AliAnalysisTaskSE(name), 
  fList(0),
  fHistEventStats(0),
  fHistTrackStats(0),
  fHistVx(0),
  fHistVy(0),
  fHistVz(0),
  fHistMultiplicity(0),
  fHistPtVsNch(0),
  fHistProtonALICE(0), fHistPionALICE(0),
  fHistLambdaALICE(0), fHistK0sALICE(0),
  fHistLambdacALICE(0), fHistDZeroALICE(0),
  fHistProtonCMS(0), fHistPionCMS(0),
  fHistLambdaCMS(0), fHistK0sCMS(0),
  fHistLambdacCMS(0), fHistDZeroCMS(0),
  fHistProtonLHCb(0), fHistPionLHCb(0),
  fHistLambdaLHCb(0), fHistK0sLHCb(0),
  fHistLambdacLHCb(0), fHistDZeroLHCb(0),
  fVxMax(0.3),
  fVyMax(0.3),
  fVzMax(10.),
  fPtMin(0.1),
  fPtMax(10.5),
  fEtaMin(-0.8),
  fEtaMax(-0.8),
  fEtaMinCMS(-2.4),
  fEtaMaxCMS(2.4),
  fEtaMinLHCb(2.),
  fEtaMaxLHCb(4.5) {
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskLcToD0MC::~AliAnalysisTaskLcToD0MC() {
  //
}

//________________________________________________________________________
void AliAnalysisTaskLcToD0MC::UserCreateOutputObjects() {
  // Create histograms
  // Called once
  //QA list
  fList = new TList();
  fList->SetName("listQA");
  fList->SetOwner();

  //Event stats.
  TString gCutName[4] = {"Total","Offline trigger",
                         "Vertex","Analyzed"};
  fHistEventStats = new TH1F("fHistEventStats",
                             "Event statistics;;N_{events}",
                             4,0.5,4.5);
  for(Int_t i = 1; i <= 4; i++)
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());
  fList->Add(fHistEventStats);

  fHistTrackStats = new TH1F("fHistTrackStats","Event statistics;TriggerBit;N_{events}",130,0,130);
  fList->Add(fHistTrackStats);

  // Vertex distributions
  fHistVx = new TH1F("fHistVx","Primary vertex distribution - x coordinate;V_{x} (cm);Entries",100,-0.5,0.5);
  fList->Add(fHistVx);
  fHistVy = new TH1F("fHistVy","Primary vertex distribution - y coordinate;V_{y} (cm);Entries",100,-0.5,0.5);
  fList->Add(fHistVy);
  fHistVz = new TH1F("fHistVz","Primary vertex distribution - z coordinate;V_{z} (cm);Entries",100,-20.,20.);
  fList->Add(fHistVz);

  // QA histograms
  fHistMultiplicity = new TH1F("fHistMultiplicity",";N_{acc.};Counts",500,-0.5,499.5);
  fList->Add(fHistMultiplicity);

  fHistPtVsNch = new TH2F("fHistPtVsNch","pT spectrum;N_{ch.};p_{T} (GeV/c);Entries",500,-0.5,499.5,1000,0.,20.);
  fList->Add(fHistPtVsNch);

  fHistProtonALICE = new TH2F("fHistProtonALICE",";p_{T};Multiplicity;Counts",100,0,50,200,-0.5,199.5);
  fHistPionALICE = new TH2F("fHistPionALICE",";p_{T};Multiplicity;Counts",100,0,50,200,-0.5,199.5);
  fHistK0sALICE = new TH2F("fHistK0sALICE",";p_{T};Multiplicity;Counts",100,0,50,200,-0.5,199.5);
  fHistLambdaALICE = new TH2F("fHistLambdaALICE",";p_{T};Multiplicity;Counts",100,0,50,200,-0.5,199.5);
  fHistDZeroALICE = new TH2F("fHistDZeroALICE",";p_{T};Multiplicity;Counts",100,0,50,200,-0.5,199.5);
  fHistLambdacALICE = new TH2F("fHistLambdacALICE",";p_{T};Multiplicity;Counts",100,0,50,200,-0.5,199.5);

  fHistProtonCMS = new TH2F("fHistProtonCMS",";p_{T};Multiplicity;Counts",100,0,50,200,-0.5,199.5);
  fHistPionCMS = new TH2F("fHistPionCMS",";p_{T};Multiplicity;Counts",100,0,50,200,-0.5,199.5);
  fHistK0sCMS = new TH2F("fHistK0sCMS",";p_{T};Multiplicity;Counts",100,0,50,200,-0.5,199.5);
  fHistLambdaCMS = new TH2F("fHistLambdaCMS",";p_{T};Multiplicity;Counts",100,0,50,200,-0.5,199.5);
  fHistDZeroCMS = new TH2F("fHistDZeroCMS",";p_{T};Multiplicity;Counts",100,0,50,200,-0.5,199.5);
  fHistLambdacCMS = new TH2F("fHistLambdacCMS",";p_{T};Multiplicity;Counts",100,0,50,200,-0.5,199.5);

  fHistProtonLHCb = new TH2F("fHistProtonLHCb",";p_{T};Multiplicity;Counts",100,0,50,200,-0.5,199.5);
  fHistPionLHCb = new TH2F("fHistPionLHCb",";p_{T};Multiplicity;Counts",100,0,50,200,-0.5,199.5);
  fHistK0sLHCb = new TH2F("fHistK0sLHCb",";p_{T};Multiplicity;Counts",100,0,50,200,-0.5,199.5);
  fHistLambdaLHCb = new TH2F("fHistLambdaLHCb",";p_{T};Multiplicity;Counts",100,0,50,200,-0.5,199.5);
  fHistDZeroLHCb = new TH2F("fHistDZeroLHCb",";p_{T};Multiplicity;Counts",100,0,50,200,-0.5,199.5);
  fHistLambdacLHCb = new TH2F("fHistLambdacLHCb",";p_{T};Multiplicity;Counts",100,0,50,200,-0.5,199.5);

  fList->Add(fHistPionALICE);
  fList->Add(fHistProtonALICE);
  fList->Add(fHistK0sALICE);
  fList->Add(fHistLambdaALICE);
  fList->Add(fHistDZeroALICE);
  fList->Add(fHistLambdacALICE);

  fList->Add(fHistPionCMS);
  fList->Add(fHistProtonCMS);
  fList->Add(fHistK0sCMS);
  fList->Add(fHistLambdaCMS);
  fList->Add(fHistDZeroCMS);
  fList->Add(fHistLambdacCMS);

  fList->Add(fHistPionLHCb);
  fList->Add(fHistProtonLHCb);
  fList->Add(fHistK0sLHCb);
  fList->Add(fHistLambdaLHCb);
  fList->Add(fHistDZeroLHCb);
  fList->Add(fHistLambdacLHCb);

  // Post output data.
  PostData(1, fList);
  
}

//________________________________________________________________________
void AliAnalysisTaskLcToD0MC::UserExec(Option_t *) {
  // Main loop
  // Called for each event
  //AOD analysis (vertex and track cuts also here!!!!)
  //AliMCEvent* gMCEvent = dynamic_cast<AliMCEvent*>(InputEvent()); // from TaskSE
  AliMCEvent* gMCEvent = dynamic_cast<AliMCEvent*>(MCEvent());
  if(!gMCEvent) {
    Printf("ERROR: gMC not available");
    return;
  }

  Int_t nAcceptedParticles = 0; //counter for multiplicity in events with Lambdac
  Int_t nAcceptedParticlesControl = 0;  //counter for multiplicity in all events
  
  // event selection done in AliAnalysisTaskSE::Exec() --> this is not used
  fHistEventStats->Fill(1); //all events
  fHistEventStats->Fill(2); //triggered events
  
  const AliVVertex *vertex = gMCEvent->GetPrimaryVertex();
  if(vertex) {
    fHistEventStats->Fill(3); //events with a proper vertex
    if(TMath::Abs(vertex->GetX()) < fVxMax) {
      if(TMath::Abs(vertex->GetY()) < fVyMax) {
	if(TMath::Abs(vertex->GetZ()) < fVzMax) {
	  fHistEventStats->Fill(4); //analyzed events
	  fHistVx->Fill(vertex->GetX());
	  fHistVy->Fill(vertex->GetY());
	  fHistVz->Fill(vertex->GetZ());
	  	  
	  //Printf("There are %d tracks in this event", gMC->GetNumberOfTracks());
	  for (Int_t iParticle = 0; iParticle < gMCEvent->GetNumberOfTracks(); iParticle++) {
	    AliMCParticle* gParticle = dynamic_cast<AliMCParticle *>(gMCEvent->GetTrack(iParticle));
	    if (!gParticle) {
	      Printf(Form("Could not receive particle %d", iParticle));
	      continue;
	    }
	    
	    
	    //kinematic cuts cuts
	    Float_t pX  = gParticle->Px();
	    Float_t pY  = gParticle->Py();
	    Float_t pT  = gParticle->Pt();
	    Float_t eta = gParticle->Eta();
	    
	    // Kinematics cuts from ESD track cuts
	    if( pT < fPtMin || pT > fPtMax)      continue;
	    if( eta < fEtaMin || eta > fEtaMax)  continue;

	    //exclude non-stable particles for counting the multiplicity
	    if(gMCEvent->IsPhysicalPrimary(iParticle)) {
	      nAcceptedParticlesControl += 1;}
	  }
 	  fHistMultiplicity->Fill(nAcceptedParticlesControl);
	
	  for (Int_t iParticle = 0; iParticle < gMCEvent->GetNumberOfTracks(); iParticle++) {
	    AliMCParticle* gParticle = dynamic_cast<AliMCParticle *>(gMCEvent->GetTrack(iParticle));
	    if (!gParticle) {
	      Printf(Form("Could not receive particle %d", iParticle));
	      continue;
	    }
	    
	    
	    //kinematic cuts cuts
	    Float_t pX  = gParticle->Px();
	    Float_t pY  = gParticle->Py();
	    Float_t pT  = gParticle->Pt();
	    Float_t eta = gParticle->Eta();	  
       
	    if(eta > fEtaMin && eta < fEtaMax)
	      fHistPtVsNch->Fill(nAcceptedParticlesControl,pT);
	    
	    //getting epg code from each particle
	    Int_t epg = gParticle -> PdgCode();

	    //==============LIGHT====================//
	    //checking for PION particle
	    if(TMath::Abs(epg) == 211){
	      //getting multiplicity for events with Dzero particle
	      if(eta > fEtaMin && eta < fEtaMax){
		fHistPionALICE -> Fill(pT,nAcceptedParticlesControl);}
	      if(eta > fEtaMinCMS && eta < fEtaMaxCMS){
		fHistPionCMS -> Fill(pT,nAcceptedParticlesControl);}
	      if(eta > fEtaMinLHCb && eta < fEtaMaxLHCb){
		fHistPionLHCb -> Fill(pT,nAcceptedParticlesControl);}
	    }

	    //checking for proton particle
	    else if(epg == 2212){
	      //getting multiplicity for events with Lambdac particle      
	      if(eta > fEtaMin && eta < fEtaMax){
		fHistProtonALICE -> Fill(pT,nAcceptedParticlesControl);}
	      if(eta > fEtaMinCMS && eta < fEtaMaxCMS){
		fHistProtonCMS -> Fill(pT,nAcceptedParticlesControl);}
	      if(eta > fEtaMinLHCb && eta < fEtaMaxLHCb){
		fHistProtonLHCb -> Fill(pT,nAcceptedParticlesControl);}
	    }

	    //==============STRANGE====================//
	    //checking for K0s particle
	    if(TMath::Abs(epg) == 310){
	      //getting multiplicity for events with Dzero particle
	      if(eta > fEtaMin && eta < fEtaMax){
		fHistK0sALICE -> Fill(pT,nAcceptedParticlesControl);}
	      if(eta > fEtaMinCMS && eta < fEtaMaxCMS){
		fHistK0sCMS -> Fill(pT,nAcceptedParticlesControl);}
	      if(eta > fEtaMinLHCb && eta < fEtaMaxLHCb){
		fHistK0sLHCb -> Fill(pT,nAcceptedParticlesControl);}
	    }

	    //checking for Lambda particle
	    else if(epg == 3122){
	      //getting multiplicity for events with Lambdac particle      
	      if(eta > fEtaMin && eta < fEtaMax){
		fHistLambdaALICE -> Fill(pT,nAcceptedParticlesControl);}
	      if(eta > fEtaMinCMS && eta < fEtaMaxCMS){
		fHistLambdaCMS -> Fill(pT,nAcceptedParticlesControl);}
	      if(eta > fEtaMinLHCb && eta < fEtaMaxLHCb){
		fHistLambdaLHCb -> Fill(pT,nAcceptedParticlesControl);}
	    }

	    //==============CHARM====================//
	    //checking for Dzero particle
	    if(TMath::Abs(epg) == 421){
	      //getting multiplicity for events with Dzero particle
	      if(eta > fEtaMin && eta < fEtaMax){
		fHistDZeroALICE -> Fill(pT,nAcceptedParticlesControl);}
	      if(eta > fEtaMinCMS && eta < fEtaMaxCMS){
		fHistDZeroCMS -> Fill(pT,nAcceptedParticlesControl);}
	      if(eta > fEtaMinLHCb && eta < fEtaMaxLHCb){
		fHistDZeroLHCb -> Fill(pT,nAcceptedParticlesControl);}
	    }

	    //checking for Lambdac particle
	    else if(epg == 4122){
	      //getting multiplicity for events with Lambdac particle      
	      if(eta > fEtaMin && eta < fEtaMax){
		fHistLambdacALICE -> Fill(pT,nAcceptedParticlesControl);}
	      if(eta > fEtaMinCMS && eta < fEtaMaxCMS){
		fHistLambdacCMS -> Fill(pT,nAcceptedParticlesControl);}
	      if(eta > fEtaMinLHCb && eta < fEtaMaxLHCb){
		fHistLambdacLHCb -> Fill(pT,nAcceptedParticlesControl);}
	    }
	    
	  }//particle loop
	}//Vz cut
      }//Vy cut
    }//Vx cut
  }//vertex object valid
}      

//________________________________________________________________________
void  AliAnalysisTaskLcToD0MC::FinishTaskOutput(){
  //
}

//________________________________________________________________________
void AliAnalysisTaskLcToD0MC::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  // not implemented ...

}
