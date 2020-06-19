#include "TChain.h"
#include "TList.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TParticle.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TVector2.h"
#include "TVector3.h"

#include "AliAODEvent.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliAODTrack.h"
#include "AliAODInputHandler.h"

#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliEventplane.h"
#include "AliCollisionGeometry.h"

#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliQnCorrectionsManager.h"

#include "AliAnalysisTaskSignedBF.h"
#include <iostream>

//Analysis task for signed BF studies
//Author: Panos.Christakoglou@cern.ch

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskSignedBF)

//________________________________________________________________________
AliAnalysisTaskSignedBF::AliAnalysisTaskSignedBF(const char *name) 
: AliAnalysisTaskSE(name), 
  fListQA(0),
  fListBF(0),
  fHistEventStats(0),
  fHistVx(0), fHistVy(0), fHistVz(0),
  fHistMultiplicityPercentile(0),
  fHistMultiplicity(0),
  fHistMultiplicityvsPercentile(0),
  fHistEventPlane(0),
  fHistNumberOfAcceptedParticles(0),
  fHistNumberOfAcceptedPositiveParticles(0),
  fHistNumberOfAcceptedNegativeParticles(0),
  fHistPt(0), 
  fFlowQnVectorMgr(0), flowQnVectorTask(0),
  fEventPlaneDetector("VZEROA"),
  fHistP(0), fHistN(0),
  fHistPNRandomOut(0), fHistNPRandomOut(0), 
  fHistPPRandomOut(0), fHistNNRandomOut(0),
  fHistPNRandomIn(0), fHistNPRandomIn(0), 
  fHistPPRandomIn(0), fHistNNRandomIn(0),
  fHistDeltaBRandomOut(0), fHistDeltaBRandomIn(0),
  fHistPNLabOut(0), fHistNPLabOut(0), 
  fHistPPLabOut(0), fHistNNLabOut(0),
  fHistPNLabIn(0), fHistNPLabIn(0), 
  fHistPPLabIn(0), fHistNNLabIn(0),
  fHistDeltaBLabOut(0), fHistDeltaBLabIn(0),
  fHistPNRestOut(0), fHistNPRestOut(0), 
  fHistPPRestOut(0), fHistNNRestOut(0),
  fHistPNRestIn(0), fHistNPRestIn(0), 
  fHistPPRestIn(0), fHistNNRestIn(0),
  fHistDeltaBRestOut(0), fHistDeltaBRestIn(0),
  fLVParticle1(0), fLVParticle2(0),
  fLVParticlePair(0),
  fV2Particle1(0),fV2Particle2(0),
  fV3Particle1(0), fV3Particle2(0),
  fUtils(0x0),
  fAnalysisLevel("AOD"),
  fMultiplicityEstimator("V0M"),
  fUseOfflineTrigger(kTRUE),
  fVxMax(3.), fVyMax(3.), fVzMax(10.),
  fUseAdditionalVtxCuts(kFALSE),
  fCheckPileUp(kFALSE),
  fCentralityPercentileMin(0.), fCentralityPercentileMax(100.),
  fnAODtrackCutBit(768),
  fPtMin(0.2), fPtMax(2.0),
  fEtaMin(-0.8),fEtaMax(0.8) {
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskSignedBF::~AliAnalysisTaskSignedBF() {
  //
}

//________________________________________________________________________
void AliAnalysisTaskSignedBF::UserCreateOutputObjects() {
  // Create histograms
  // Called once

  // global switch disabling the reference 
  // (to avoid "Replacing existing TH1" if several wagons are created in train)
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  ///Output list
  fListQA = new TList();
  fListQA->SetName("listQA");
  fListQA->SetOwner();
  
  fListBF = new TList();
  fListBF->SetName("listBF");
  fListBF->SetOwner();

  //================Event level=================//
  //Event stats.
  TString gCutName[7] = {"Total","Offline trigger",
			 "Vertex","Add Vtx Cuts", "Pile up", "sel. Centrality",""};
  fHistEventStats = new TH1F("fHistEventStats",
                             "Event statistics;;N_{events}",
                             7,0.5,7.5);
  for(Int_t i = 1; i <= 7; i++)
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());
  fListQA->Add(fHistEventStats);  

  fHistVx = new TH1F("fHistVx","Primary vertex distribution - x coordinate;V_{x} (cm);Entries",100,-0.5,0.5);
  fListQA->Add(fHistVx);
  fHistVy = new TH1F("fHistVy","Primary vertex distribution - y coordinate;V_{y} (cm);Entries",100,-0.5,0.5);
  fListQA->Add(fHistVy);
  fHistVz = new TH1F("fHistVz","Primary vertex distribution - z coordinate;V_{z} (cm);Centrality percentile;Entries",100,-20.,20.);
  fListQA->Add(fHistVz);

  fHistMultiplicityPercentile = new TH1F("fHistMultiplicityPercentile","Multiplicity percentile; Event class (%);Entries",103,-1.5,100.5);
  fListQA->Add(fHistMultiplicityPercentile);

  fHistMultiplicity = new TH1F("fHistMultiplicity","Multiplicity distribution; Multiplicity;Entries",100002,-1.,20000.);
  fListQA->Add(fHistMultiplicity);

  fHistMultiplicityvsPercentile = new TH2F("fHistMultiplicityvsPercentile","Multiplicity vs percentile; Multiplicity;Event class (%)",100002,-1.,20000.,103,-1.5,100.5);
  fListQA->Add(fHistMultiplicityvsPercentile);

  //Event plane
  fHistEventPlane = new TH2F("fHistEventPlane",";#Psi_{2} [deg.];Centrality percentile;Counts",100,0,TMath::Pi(),220,-5,105);
  fListQA->Add(fHistEventPlane);

  fHistNumberOfAcceptedParticles = new TH1F("fHistNumberOfAcceptedParticles",";N_{acc.};Entries",10000,0,10000);
  fListQA->Add(fHistNumberOfAcceptedParticles);
  fHistNumberOfAcceptedPositiveParticles = new TH1F("fHistNumberOfAcceptedPositiveParticles",";N_{acc.}^{+};Entries",10000,0,10000);
  fListQA->Add(fHistNumberOfAcceptedPositiveParticles);
  fHistNumberOfAcceptedNegativeParticles = new TH1F("fHistNumberOfAcceptedNegativeParticles",";N_{acc.}^{-};Entries",10000,0,10000);
  fListQA->Add(fHistNumberOfAcceptedNegativeParticles);
  //============================================//

  //================Track level=================//
  fHistPt = new TH1F("fHistPt","pT spectrum; p_{T} (GeV/c);",100,0,10.);
  fListQA->Add(fHistPt);
  //============================================//

  //============================================//
  //AliQnCorrectionsManager *fFlowQnVectorMgr;
  flowQnVectorTask = dynamic_cast<AliAnalysisTaskFlowVectorCorrections *>(AliAnalysisManager::GetAnalysisManager()->GetTask("FlowQnVectorCorrections"));
  //AliAnalysisManager *mgr = (AliAnalysisManager::GetAnalysisManager());
  //flowQnVectorTask = dynamic_cast<AliAnalysisTaskFlowVectorCorrections *>(mgr->GetTask("FlowQnVectorCorrections"));
  if (flowQnVectorTask != NULL) {
    fFlowQnVectorMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
  }
  else {
    AliFatal("Flow Qn vector corrections framework needed but it is not present. ABORTING!!!");
  }
  //============================================//
  
  //============================================//
  fHistP = new TH1F("fHistP","P;p_{y} (GeV/c);Entries",1,-1.,1.);
  fHistN = new TH1F("fHistN","N;p_{y} (GeV/c);Entries",1,-1.,1.);
  //============================================//

  //============================================//
  //BF - random frame
  fHistPNRandomOut = new TH1F("fHistPNRandomOut","PN;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);
  fHistNPRandomOut = new TH1F("fHistNPRandomOut","NP;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);
  fHistPPRandomOut = new TH1F("fHistPPRandomOut","PP;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);
  fHistNNRandomOut = new TH1F("fHistNNRandomOut","NN;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);

  fHistPNRandomIn = new TH1F("fHistPNRandomIn","PN;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);
  fHistNPRandomIn = new TH1F("fHistNPRandomIn","NP;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);
  fHistPPRandomIn = new TH1F("fHistPPRandomIn","PP;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);
  fHistNNRandomIn = new TH1F("fHistNNRandomIn","NN;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);
  
  fHistDeltaBRandomOut = new TH2F("fHistDeltaBRandomOut","#Delta B out of plane (lab frame);#Delta B;Entries",1000,-99.999,99.999,102,-1.,101.);
  fListBF->Add(fHistDeltaBRandomOut);
  fHistDeltaBRandomIn = new TH2F("fHistDeltaBRandomIn","#Delta B in plane (lab frame);#Delta B;Entries",1000,-99.999,99.999,102,-1.,101.);
  fListBF->Add(fHistDeltaBRandomIn);
  //============================================//

  //============================================//
  //BF - lab frame
  fHistPNLabOut = new TH1F("fHistPNLabOut","PN;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);
  fHistNPLabOut = new TH1F("fHistNPLabOut","NP;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);
  fHistPPLabOut = new TH1F("fHistPPLabOut","PP;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);
  fHistNNLabOut = new TH1F("fHistNNLabOut","NN;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);

  fHistPNLabIn = new TH1F("fHistPNLabIn","PN;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);
  fHistNPLabIn = new TH1F("fHistNPLabIn","NP;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);
  fHistPPLabIn = new TH1F("fHistPPLabIn","PP;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);
  fHistNNLabIn = new TH1F("fHistNNLabIn","NN;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);
  
  fHistDeltaBLabOut = new TH2F("fHistDeltaBLabOut","#Delta B out of plane (lab frame);#Delta B;Entries",1000,-99.999,99.999,102,-1.,101.);
  fListBF->Add(fHistDeltaBLabOut);
  fHistDeltaBLabIn = new TH2F("fHistDeltaBLabIn","#Delta B in plane (lab frame);#Delta B;Entries",1000,-99.999,99.999,102,-1.,101.);
  fListBF->Add(fHistDeltaBLabIn);
  //============================================//

  //============================================//
  //BF - rest frame
  fHistPNRestOut = new TH1F("fHistPNRestOut","PN;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);
  fHistNPRestOut = new TH1F("fHistNPRestOut","NP;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);
  fHistPPRestOut = new TH1F("fHistPPRestOut","PP;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);
  fHistNNRestOut = new TH1F("fHistNNRestOut","NN;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);

  fHistPNRestIn = new TH1F("fHistPNRestIn","PN;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);
  fHistNPRestIn = new TH1F("fHistNPRestIn","NP;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);
  fHistPPRestIn = new TH1F("fHistPPRestIn","PP;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);
  fHistNNRestIn = new TH1F("fHistNNRestIn","NN;#Delta p_{y} (GeV/c);Entries",2,-1.,1.);

  fHistDeltaBRestOut = new TH2F("fHistDeltaBRestOut","#Delta B out of plane (rest frame);#Delta B;Entries",1000,-99.999,99.999,102,-1.,101.0);
  fListBF->Add(fHistDeltaBRestOut);
  fHistDeltaBRestIn = new TH2F("fHistDeltaBRestIn","#Delta B in plane (rest frame);#Delta B;Entries",1000,-99.999,99.999,102,-1.,101.);
  fListBF->Add(fHistDeltaBRestIn);
  //============================================//

  fLVParticle1 = new TLorentzVector();
  fLVParticle2 = new TLorentzVector();
  fLVParticlePair = new TLorentzVector();

  fV2Particle1 = new TVector2();
  fV2Particle2 = new TVector2();

  fV3Particle1 = new TVector3();
  fV3Particle2 = new TVector3();

  PostData(1, fListQA);
  PostData(2, fListBF);

  fUtils = new AliAnalysisUtils();
}

//________________________________________________________________________
void AliAnalysisTaskSignedBF::UserExec(Option_t *) {
  // Main loop
  // Called for each event
  //Int_t gNumberOfAcceptedTracks = 0;
  Double_t lMultiplicityVar     = -999.; //-1
  Double_t gReactionPlane       = -1.; 
  //Float_t bSign = 0.;

  // get the event (for generator level: MCEvent())
  AliVEvent* eventMain = NULL;
  if(fAnalysisLevel == "MC") {
    eventMain = dynamic_cast<AliVEvent*>(MCEvent());
  }
  else{
    eventMain = dynamic_cast<AliVEvent*>(InputEvent());
    // for HBT like cuts need magnetic field sign
    //bSign = (eventMain->GetMagneticField() > 0) ? 1 : -1;
  }
  if(!eventMain) {
    AliError("eventMain not available");
    return;
  }
  
  if((lMultiplicityVar = IsEventAccepted(eventMain)) < 0) { 
    return;
  }
    
  gReactionPlane = GetEventPlane(eventMain);
  fHistEventPlane->Fill(gReactionPlane,lMultiplicityVar);
  
  //Reset histograms
  fHistPNRandomOut->Reset(); fHistNPRandomOut->Reset();
  fHistPPRandomOut->Reset(); fHistNNRandomOut->Reset();
  fHistPNRandomIn->Reset(); fHistNPRandomIn->Reset();
  fHistPPRandomIn->Reset(); fHistNNRandomIn->Reset();

  fHistPNLabOut->Reset(); fHistNPLabOut->Reset();
  fHistPPLabOut->Reset(); fHistNNLabOut->Reset();
  fHistPNLabIn->Reset(); fHistNPLabIn->Reset();
  fHistPPLabIn->Reset(); fHistNNLabIn->Reset();

  fHistPNRestOut->Reset(); fHistNPRestOut->Reset();
  fHistPPRestOut->Reset(); fHistNNRestOut->Reset();
  fHistPNRestIn->Reset(); fHistNPRestIn->Reset();
  fHistPPRestIn->Reset(); fHistNNRestIn->Reset();
  
  fHistP->Reset(); fHistN->Reset();
  
  //process the event
  ProcessEvent(eventMain, lMultiplicityVar, gReactionPlane);

  PostData(1, fListQA);
  PostData(2, fListBF);
}

//________________________________________________________________________
void  AliAnalysisTaskSignedBF::FinishTaskOutput() {
  //
}

//________________________________________________________________________
void AliAnalysisTaskSignedBF::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  // not implemented ...

}

//________________________________________________________________________
void AliAnalysisTaskSignedBF::FillBFHistograms(TObjArray *tracksAccepted,
					       Double_t gCentrality,
					       Double_t gReactionPlane,
					       Int_t i1) {
  //Fill in the BF histograms
  if(!tracksAccepted) return;
  
  TParticle *particle1 = dynamic_cast<TParticle *>(tracksAccepted->At(i1));
  if(!particle1) return;
  fV3Particle1->SetXYZ(particle1->Px(),
		       particle1->Py(),
		       particle1->Pz());
  //Rotate wrt reaction plane
  fV3Particle1->RotateZ(gReactionPlane);
  fV2Particle1->Set(fV3Particle1->X(),
		    fV3Particle1->Z());
  
  for(Int_t i2 = 0; i2 < tracksAccepted->GetEntries(); i2++) {
    if(i1 == i2) continue;
    TParticle *particle2 = dynamic_cast<TParticle *>(tracksAccepted->At(i2));
    if(!particle2) continue;
    
    //======================================//
    //Rotate wrt reaction plane
    fV3Particle2->SetXYZ(particle2->Px(),
			 particle2->Py(),
			 particle2->Pz());
    fV3Particle2->RotateZ(gReactionPlane);
    fV2Particle2->Set(fV3Particle2->X(),
		      fV3Particle2->Z());
    //======================================//
    
    //======================================//
    //Boost to pair's rest frame from the Psi_RP rotated vectors
    fLVParticle1->SetPxPyPzE(fV3Particle1->X(),
			     fV3Particle1->Y(),
			     fV3Particle1->Z(),
			     particle1->Energy());
    
    fLVParticle2->SetPxPyPzE(fV3Particle2->X(),
			     fV3Particle2->Y(),
			     fV3Particle2->Z(),
			     particle2->Energy());
    
    fLVParticlePair->SetPxPyPzE(fV3Particle1->X() + fV3Particle2->X(),
				fV3Particle1->Y() + fV3Particle2->Y(),
				fV3Particle1->Z() + fV3Particle2->Z(),
				particle1->Energy() + particle2->Energy());
    
    TVector3 boostVector = fLVParticlePair->BoostVector();
    
    fLVParticle1->Boost(-boostVector);
    fLVParticle2->Boost(-boostVector);
    //======================================//
    
    Double_t deltaPx = particle1->Px() - particle2->Px();
    Double_t deltaPy = particle1->Py() - particle2->Py();
    
    Double_t deltaPxPrime = fV3Particle1->X() - fV3Particle2->X();
    Double_t deltaPyPrime = fV3Particle1->Y() - fV3Particle2->Y();
    
    Double_t deltaPxRest = fLVParticle1->Px() - fLVParticle2->Px();
    Double_t deltaPyRest = fLVParticle1->Py() - fLVParticle2->Py();
    
    if((particle1->GetPdgCode() > 0)&&(particle2->GetPdgCode() < 0)) {
      //================================================//
      //in-plane - random
      if(deltaPx < 0) {
	fHistPNRandomIn->Fill(-0.5); //fHistNPRandomIn->Fill(0.5);
      }
      else if(deltaPx >= 0) {
	fHistPNRandomIn->Fill(0.5); //fHistNPRandomIn->Fill(-0.5);
      }
      
      //out-of-plane - random
      if(deltaPy < 0) {
	fHistPNRandomOut->Fill(-0.5); //fHistNPRandomOut->Fill(0.5);
      }
      else if(deltaPy >= 0) {
	fHistPNRandomOut->Fill(0.5); //fHistNPRandomOut->Fill(-0.5);
      }

      //================================================//
      //in-plane - lab
      if(deltaPxPrime < 0) {
	fHistPNLabIn->Fill(-0.5); //fHistNPLabIn->Fill(0.5);
      }
      else if(deltaPxPrime >= 0) {
	fHistPNLabIn->Fill(0.5); //fHistNPLabIn->Fill(-0.5);
      }
      
      //out-of-plane - lab
      if(deltaPyPrime < 0) {
	fHistPNLabOut->Fill(-0.5); //fHistNPLabOut->Fill(0.5);
      }
      else if(deltaPyPrime >= 0) {
	fHistPNLabOut->Fill(0.5); //fHistNPLabOut->Fill(-0.5);
      }

      //================================================//
      //in-plane - pair rest frame
      if(deltaPxRest < 0) {
	fHistPNRestIn->Fill(-0.5); //fHistNPRestIn->Fill(0.5);
      }
      else if(deltaPxRest >= 0) {
	fHistPNRestIn->Fill(0.5); //fHistNPRestIn->Fill(-0.5);
      }
      
      //out-of-plane - pair rest frame
      if(deltaPyRest < 0) {
	fHistPNRestOut->Fill(-0.5); //fHistNPRestOut->Fill(0.5);
      }
      else if(deltaPyRest >= 0) {
	fHistPNRestOut->Fill(0.5); //fHistNPRestOut->Fill(-0.5);
      }
    }
    else if((particle1->GetPdgCode() < 0)&&(particle2->GetPdgCode() > 0)) {
      //================================================//
      //in-plane - random
      if(deltaPx < 0) {
	fHistNPRandomIn->Fill(-0.5); //fHistPNRandomIn->Fill(0.5);
      }
      else if(deltaPx >= 0) {
	fHistNPRandomIn->Fill(0.5); //fHistPNRandomIn->Fill(-0.5);
      }
      
      //out-of-plane - random
      if(deltaPy < 0) {
	fHistNPRandomOut->Fill(-0.5); //fHistPNRandomOut->Fill(0.5);
      }
      else if(deltaPy >= 0) {
	fHistNPRandomOut->Fill(0.5); //fHistPNRandomOut->Fill(-0.5);
      }

      //================================================//
      //in-plane - lab
      if(deltaPxPrime < 0) {
	fHistNPLabIn->Fill(-0.5); //fHistPNLabIn->Fill(0.5);
      }
      else if(deltaPxPrime >= 0) {
	fHistNPLabIn->Fill(0.5); //fHistPNLabIn->Fill(-0.5);
      }
      
      //out-of-plane - lab
      if(deltaPyPrime < 0) {
	fHistNPLabOut->Fill(-0.5); //fHistPNLabOut->Fill(0.5);
      }
      else if(deltaPyPrime >= 0) {
	fHistNPLabOut->Fill(0.5); //fHistPNLabOut->Fill(-0.5);
      }

      //================================================//
      //in-plane - pair rest frame
      if(deltaPxRest < 0) {
	fHistNPRestIn->Fill(-0.5); //fHistPNRestIn->Fill(0.5);
      }
      else if(deltaPxRest >= 0) {
	fHistNPRestIn->Fill(0.5); //fHistPNRestIn->Fill(-0.5);
      }
      
      //out-of-plane - pair rest frame
      if(deltaPyRest < 0) {
	fHistNPRestOut->Fill(-0.5); //fHistPNRestOut->Fill(0.5);
      }
      else if(deltaPyRest >= 0) {
	fHistNPRestOut->Fill(0.5); //fHistPNRestOut->Fill(-0.5);
      }
    }
    else if((particle1->GetPdgCode() > 0)&&(particle2->GetPdgCode() > 0)) {
      //================================================//
      //in-plane - random
      if(deltaPx < 0) {
	fHistPPRandomIn->Fill(-0.5); //fHistPPRandomIn->Fill(0.5);
      }
      else if(deltaPx >= 0) {
	fHistPPRandomIn->Fill(0.5); //fHistPPRandomIn->Fill(-0.5);
      }
      
      //out-of-plane - random
      if(deltaPy < 0) {
	fHistPPRandomOut->Fill(-0.5); //fHistPPRandomOut->Fill(0.5);
      }
      else if(deltaPy >= 0) {
	fHistPPRandomOut->Fill(0.5); //fHistPPRandomOut->Fill(-0.5);
      }

      //================================================//
      //in-plane - lab
      if(deltaPxPrime < 0) {
	fHistPPLabIn->Fill(-0.5); //fHistPPLabIn->Fill(0.5);
      }
      else if(deltaPxPrime >= 0) {
	fHistPPLabIn->Fill(0.5); //fHistPPLabIn->Fill(-0.5);
      }
      
      //out-of-plane - lab
      if(deltaPyPrime < 0) {
	fHistPPLabOut->Fill(-0.5); //fHistPPLabOut->Fill(0.5);
      }
      else if(deltaPyPrime >= 0) {
	fHistPPLabOut->Fill(0.5); //fHistPPLabOut->Fill(-0.5);
      }

      //================================================//
      //in-plane - pair rest frame
      if(deltaPxRest < 0) {
	fHistPPRestIn->Fill(-0.5); //fHistPPRestIn->Fill(0.5);
      }
      else if(deltaPxRest >= 0) {
	fHistPPRestIn->Fill(0.5); //fHistPPRestIn->Fill(-0.5);
      }
      
      //out-of-plane - pair rest frame
      if(deltaPyRest < 0) {
	fHistPPRestOut->Fill(-0.5); //fHistPPRestOut->Fill(0.5);
      }
      else if(deltaPyRest >= 0) {
	fHistPPRestOut->Fill(0.5); //fHistPPRestOut->Fill(-0.5);
      }
    }
    else if((particle1->GetPdgCode() < 0)&&(particle2->GetPdgCode() < 0)) {
      //================================================//
      //in-plane - random
      if(deltaPx < 0) {
	fHistNNRandomIn->Fill(-0.5); //fHistNNRandomIn->Fill(0.5);
      }
      else if(deltaPx >= 0) {
	fHistNNRandomIn->Fill(0.5); //fHistNNRandomIn->Fill(-0.5);
      }
      
      //out-of-plane - random
      if(deltaPy < 0) {
	fHistNNRandomOut->Fill(-0.5); //fHistNNRandomOut->Fill(0.5);
      }
      else if(deltaPy >= 0) {
	fHistNNRandomOut->Fill(0.5); //fHistNNRandomOut->Fill(-0.5);
      }

      //================================================//
      //in-plane - lab
      if(deltaPxPrime < 0) {
	fHistNNLabIn->Fill(-0.5); //fHistNNLabIn->Fill(0.5);
      }
      else if(deltaPxPrime >= 0) {
	fHistNNLabIn->Fill(0.5); //fHistNNLabIn->Fill(-0.5);
      }
      
      //out-of-plane - lab
      if(deltaPyPrime < 0) {
	fHistNNLabOut->Fill(-0.5); //fHistNNLabOut->Fill(0.5);
      }
      else if(deltaPyPrime >= 0) {
	fHistNNLabOut->Fill(0.5); //fHistNNLabOut->Fill(-0.5);
      }

      //================================================//
      //in-plane - pair rest frame
      if(deltaPxRest < 0) {
	fHistNNRestIn->Fill(-0.5); //fHistNNRestIn->Fill(0.5);
      }
      else if(deltaPxRest >= 0) {
	fHistNNRestIn->Fill(0.5); //fHistNNRestIn->Fill(-0.5);
      }
      
      //out-of-plane - pair rest frame
      if(deltaPyRest < 0) {
	fHistNNRestOut->Fill(-0.5); //fHistNNRestOut->Fill(0.5);
      }
      else if(deltaPyRest >= 0) {
	fHistNNRestOut->Fill(0.5); //fHistNNRestOut->Fill(-0.5);
      }
    }//(--)
  }//2nd particle loop
  if(particle1->GetPdgCode() > 0) fHistP->Fill(0);
  else if(particle1->GetPdgCode() < 0) fHistN->Fill(0);
}

//________________________________________________________________________
void AliAnalysisTaskSignedBF::ProcessEvent(AliVEvent *event,
					   Double_t gCentrality,
					   Double_t gReactionPlane) {
  // Loop over tracks in event
  Short_t vCharge = 0;
  Float_t vEta = 0.0;
  Float_t vP[3] = {0.,0.,0.};
  Float_t vPt = 0.0;
  Float_t vE = 0.0;

  Int_t gNumberOfAcceptedPositiveParticles = 0;
  Int_t gNumberOfAcceptedNegativeParticles = 0;
  Int_t gNumberOfAcceptedParticles = 0;
  Int_t i1 = 0;

  TParticle *pion = new TParticle();
  pion->SetPdgCode(211);
  Double_t gPionMass = pion->GetMass();

  Double_t gParticleMass = gPionMass;

  //output TObjArray holding all good tracks
  TObjArray* tracksAccepted = new TObjArray;
  tracksAccepted->SetOwner(kTRUE);

  for (Int_t iTracks = 0; iTracks < event->GetNumberOfTracks(); iTracks++) {
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack *>(event->GetTrack(iTracks));
    if (!aodTrack) {
      AliError(Form("Could not receive track %d", iTracks));
      continue;
    }
    
    // AOD track cuts
    if(!aodTrack->TestFilterBit(fnAODtrackCutBit)) continue;
    
    vEta = aodTrack->Eta();
    vPt = aodTrack->Pt();
    vCharge = aodTrack->Charge();

    //Acceptance 
    if((vEta < fEtaMin) || (vEta > fEtaMax)) continue;
    
    //pt coverage
    if((vPt < fPtMin) || (vPt > fPtMax)) continue;

    vP[0] = aodTrack->Px();
    vP[1] = aodTrack->Py();
    vP[2] = aodTrack->Pz();
    vE = TMath::Sqrt(TMath::Power(gParticleMass,2) +
		     TMath::Power(vP[0],2) +
		     TMath::Power(vP[1],2) +
		     TMath::Power(vP[2],2));

    fHistPt->Fill(aodTrack->Pt());

    if(vCharge > 0) {
      tracksAccepted->Add(new TParticle(211,0,-1,-1,-1,-1,vP[0],vP[1],vP[2],vE,0,0,0,0));
      gNumberOfAcceptedPositiveParticles += 1;
    }
    else {
      tracksAccepted->Add(new TParticle(-211,0,-1,-1,-1,-1,vP[0],vP[1],vP[2],vE,0,0,0,0));
      gNumberOfAcceptedNegativeParticles += 1;
    }

    //Fill the BF histograms
    FillBFHistograms(tracksAccepted, gCentrality, gReactionPlane, i1); 
    
    gNumberOfAcceptedParticles += 1;
    i1 += 1;
  }//loop over AOD tracks
  
  fHistNumberOfAcceptedParticles->Fill(gNumberOfAcceptedParticles);
  fHistNumberOfAcceptedPositiveParticles->Fill(gNumberOfAcceptedPositiveParticles);
  fHistNumberOfAcceptedNegativeParticles->Fill(gNumberOfAcceptedNegativeParticles);
  
  //=====================================//0
  //random frame
  Double_t bfPRandomIn = -999., bfNRandomIn = -999.;
  Double_t bfPRandomOut = -999., bfNRandomOut = -999.;
  if((fHistP->GetEntries() != 0)&&(fHistN->GetEntries() != 0)) {
    bfNRandomIn = (fHistPNRandomIn->GetBinContent(1) - fHistPPRandomIn->GetBinContent(1))/fHistP->GetEntries() - (fHistNPRandomIn->GetBinContent(1) - fHistNNRandomIn->GetBinContent(1))/fHistN->GetEntries();
    bfPRandomIn = (fHistPNRandomIn->GetBinContent(2) - fHistPPRandomIn->GetBinContent(2))/fHistP->GetEntries() - (fHistNPRandomIn->GetBinContent(2) - fHistNNRandomIn->GetBinContent(2))/fHistN->GetEntries();
    
    bfNRandomOut = (fHistPNRandomOut->GetBinContent(1) - fHistPPRandomOut->GetBinContent(1))/fHistP->GetEntries() - (fHistNPRandomOut->GetBinContent(1) - fHistNNRandomOut->GetBinContent(1))/fHistN->GetEntries();
    bfPRandomOut = (fHistPNRandomOut->GetBinContent(2) - fHistPPRandomOut->GetBinContent(2))/fHistP->GetEntries() - (fHistNPRandomOut->GetBinContent(2) - fHistNNRandomOut->GetBinContent(2))/fHistN->GetEntries();
  }
  
  Double_t deltaBRandomIn = bfPRandomIn - bfNRandomIn;
  fHistDeltaBRandomIn->Fill(deltaBRandomIn,gCentrality);
  
  Double_t deltaBRandomOut = bfPRandomOut - bfNRandomOut;
  fHistDeltaBRandomOut->Fill(deltaBRandomOut,gCentrality);

  //=====================================//0
  //lab frame
  Double_t bfPLabIn = -999., bfNLabIn = -999.;
  Double_t bfPLabOut = -999., bfNLabOut = -999.;
  if((fHistP->GetEntries() != 0)&&(fHistN->GetEntries() != 0)) {
    bfNLabIn = (fHistPNLabIn->GetBinContent(1) - fHistPPLabIn->GetBinContent(1))/fHistP->GetEntries() - (fHistNPLabIn->GetBinContent(1) - fHistNNLabIn->GetBinContent(1))/fHistN->GetEntries();
    bfPLabIn = (fHistPNLabIn->GetBinContent(2) - fHistPPLabIn->GetBinContent(2))/fHistP->GetEntries() - (fHistNPLabIn->GetBinContent(2) - fHistNNLabIn->GetBinContent(2))/fHistN->GetEntries();
    
    bfNLabOut = (fHistPNLabOut->GetBinContent(1) - fHistPPLabOut->GetBinContent(1))/fHistP->GetEntries() - (fHistNPLabOut->GetBinContent(1) - fHistNNLabOut->GetBinContent(1))/fHistN->GetEntries();
    bfPLabOut = (fHistPNLabOut->GetBinContent(2) - fHistPPLabOut->GetBinContent(2))/fHistP->GetEntries() - (fHistNPLabOut->GetBinContent(2) - fHistNNLabOut->GetBinContent(2))/fHistN->GetEntries();
  }
  
  Double_t deltaBLabIn = bfPLabIn - bfNLabIn;
  fHistDeltaBLabIn->Fill(deltaBLabIn,gCentrality);
  
  Double_t deltaBLabOut = bfPLabOut - bfNLabOut;
  fHistDeltaBLabOut->Fill(deltaBLabOut,gCentrality);

  //=====================================//
  //pair's rest frame
  Double_t bfPRestIn = -999., bfNRestIn = -999.;
  Double_t bfPRestOut = -999., bfNRestOut = -999.;
  if((fHistP->GetEntries() != 0)&&(fHistN->GetEntries() != 0)) {
    bfNRestIn = (fHistPNRestIn->GetBinContent(1) - fHistPPRestIn->GetBinContent(1))/fHistP->GetEntries() - (fHistNPRestIn->GetBinContent(1) - fHistNNRestIn->GetBinContent(1))/fHistN->GetEntries();
    bfPRestIn = (fHistPNRestIn->GetBinContent(2) - fHistPPRestIn->GetBinContent(2))/fHistP->GetEntries() - (fHistNPRestIn->GetBinContent(2) - fHistNNRestIn->GetBinContent(2))/fHistN->GetEntries();
    
    bfNRestOut = (fHistPNRestOut->GetBinContent(1) - fHistPPRestOut->GetBinContent(1))/fHistP->GetEntries() - (fHistNPRestOut->GetBinContent(1) - fHistNNRestOut->GetBinContent(1))/fHistN->GetEntries();
    bfPRestOut = (fHistPNRestOut->GetBinContent(2) - fHistPPRestOut->GetBinContent(2))/fHistP->GetEntries() - (fHistNPRestOut->GetBinContent(2) - fHistNNRestOut->GetBinContent(2))/fHistN->GetEntries();
  }
  
  Double_t deltaBRestIn = bfPRestIn - bfNRestIn;
  fHistDeltaBRestIn->Fill(deltaBRestIn,gCentrality);
  
  Double_t deltaBRestOut = bfPRestOut - bfNRestOut;
  fHistDeltaBRestOut->Fill(deltaBRestOut,gCentrality);
}

//________________________________________________________________________
void AliAnalysisTaskSignedBF::CalculateSignedBFEbyE(TObjArray *cObjAcceptedParticles,
						    Double_t gCentrality,
						    Double_t gReactionPlane) {
  //calculate the ebye balance function
  for(Int_t i1 = 0; i1 < cObjAcceptedParticles->GetEntries(); i1++) {
    TParticle *particle1 = dynamic_cast<TParticle *>(cObjAcceptedParticles->At(i1));
    if(!particle1) continue;
    fV3Particle1->SetXYZ(particle1->Px(),
			 particle1->Py(),
			 particle1->Pz());
    //Rotate wrt reaction plane
    fV3Particle1->RotateZ(gReactionPlane);
    fV2Particle1->Set(fV3Particle1->X(),
		      fV3Particle1->Z());
      
    for(Int_t i2 = 0; i2 < cObjAcceptedParticles->GetEntries(); i2++) {
      if(i1 == i2) continue;
      TParticle *particle2 = dynamic_cast<TParticle *>(cObjAcceptedParticles->At(i2));
      if(!particle2) continue;
      
      //======================================//
      //Rotate wrt reaction plane
      fV3Particle2->SetXYZ(particle2->Px(),
			   particle2->Py(),
			   particle2->Pz());
      fV3Particle2->RotateZ(gReactionPlane);
      fV2Particle2->Set(fV3Particle2->X(),
			fV3Particle2->Z());
      //======================================//
      
      //======================================//
      //Boost to pair's rest frame from the Psi_RP rotated vectors
      fLVParticle1->SetPxPyPzE(fV3Particle1->X(),
			       fV3Particle1->Y(),
			       fV3Particle1->Z(),
			       particle1->Energy());

      fLVParticle2->SetPxPyPzE(fV3Particle2->X(),
			       fV3Particle2->Y(),
			       fV3Particle2->Z(),
			       particle2->Energy());

      fLVParticlePair->SetPxPyPzE(fV3Particle1->X() + fV3Particle2->X(),
				  fV3Particle1->Y() + fV3Particle2->Y(),
				  fV3Particle1->Z() + fV3Particle2->Z(),
				  particle1->Energy() + particle2->Energy());

      TVector3 boostVector = fLVParticlePair->BoostVector();
	
      fLVParticle1->Boost(-boostVector);
      fLVParticle2->Boost(-boostVector);
      //======================================//
	
      //Double_t deltaPx = particle1->Px() - particle2->Px();
      //Double_t deltaPy = particle1->Py() - particle2->Py();
      
      Double_t deltaPxPrime = fV3Particle1->X() - fV3Particle2->X();
      Double_t deltaPyPrime = fV3Particle1->Y() - fV3Particle2->Y();
      
      Double_t deltaPxRest = fLVParticle1->Px() - fLVParticle2->Px();
      Double_t deltaPyRest = fLVParticle1->Py() - fLVParticle2->Py();
      
      if((particle1->GetPdgCode() > 0)&&(particle2->GetPdgCode() < 0)) {
	//in-plane - lab
	if(deltaPxPrime < 0) fHistPNLabIn->Fill(-0.5);
	else if(deltaPxPrime >= 0) fHistPNLabIn->Fill(0.5);
	
	//out-of-plane - lab
	if(deltaPyPrime < 0) fHistPNLabOut->Fill(-0.5);
	else if(deltaPyPrime >= 0) fHistPNLabOut->Fill(0.5);
	
	//in-plane - pair rest frame
	if(deltaPxRest < 0) fHistPNRestIn->Fill(-0.5);
	else if(deltaPxRest >= 0) fHistPNRestIn->Fill(0.5);
	
	//out-of-plane - pair rest frame
	if(deltaPyRest < 0) fHistPNRestOut->Fill(-0.5);
	else if(deltaPyRest >= 0) fHistPNRestOut->Fill(0.5);
      }
      else if((particle1->GetPdgCode() < 0)&&(particle2->GetPdgCode() > 0)) {
	//in-plane - lab
	if(deltaPxPrime < 0) fHistNPLabIn->Fill(-0.5);
	else if(deltaPxPrime >= 0) fHistNPLabIn->Fill(0.5);
	
	//out-of-plane - lab
	if(deltaPyPrime < 0) fHistNPLabOut->Fill(-0.5);
	else if(deltaPyPrime >= 0) fHistNPLabOut->Fill(0.5);
	
	//in-plane - pair rest frame
	if(deltaPxRest < 0) fHistNPRestIn->Fill(-0.5);
	else if(deltaPxRest >= 0) fHistNPRestIn->Fill(0.5);
	
	//out-of-plane - pair rest frame
	if(deltaPyRest < 0) fHistNPRestOut->Fill(-0.5);
	else if(deltaPyRest >= 0) fHistNPRestOut->Fill(0.5);
      }
      else if((particle1->GetPdgCode() > 0)&&(particle2->GetPdgCode() > 0)) {
	//in-plane - lab
	if(deltaPxPrime < 0) fHistPPLabIn->Fill(-0.5);
	else if(deltaPxPrime >= 0) fHistPPLabIn->Fill(0.5);
	
	//out-of-plane - lab
	if(deltaPyPrime < 0) fHistPPLabOut->Fill(-0.5);
	else if(deltaPyPrime >= 0) fHistPPLabOut->Fill(0.5);
	
	//in-plane - pair rest frame
	if(deltaPxRest < 0) fHistPPRestIn->Fill(-0.5);
	else if(deltaPxRest >= 0) fHistPPRestIn->Fill(0.5);
	
	//out-of-plane - pair rest frame
	if(deltaPyRest < 0) fHistPPRestOut->Fill(-0.5);
	else if(deltaPyRest >= 0) fHistPPRestOut->Fill(0.5);
      }
      else if((particle1->GetPdgCode() < 0)&&(particle2->GetPdgCode() < 0)) {
	//in-plane - lab
	if(deltaPxPrime < 0) fHistNNLabIn->Fill(-0.5);
	else if(deltaPxPrime >= 0) fHistNNLabIn->Fill(0.5);
	
	//out-of-plane - lab
	if(deltaPyPrime < 0) fHistNNLabOut->Fill(-0.5);
	else if(deltaPyPrime >= 0) fHistNNLabOut->Fill(0.5);
	
	//in-plane - pair rest frame
	if(deltaPxRest < 0) fHistNNRestIn->Fill(-0.5);
	else if(deltaPxRest >= 0) fHistNNRestIn->Fill(0.5);
	
	//out-of-plane - pair rest frame
	if(deltaPyRest < 0) fHistNNRestOut->Fill(-0.5);
	else if(deltaPyRest >= 0) fHistNNRestOut->Fill(0.5);
      }
    }//2nd particle loop
    
    if(particle1->GetPdgCode() > 0) fHistP->Fill(0);
    else if(particle1->GetPdgCode() < 0) fHistN->Fill(0);
  }//1st particle loop

  //=====================================//
  //lab frame
  Double_t bfPLabIn = -999., bfNLabIn = -999.;
  Double_t bfPLabOut = -999., bfNLabOut = -999.;
  if((fHistP->GetEntries() != 0)&&(fHistN->GetEntries() != 0)) {
    bfNLabIn = (fHistPNLabIn->GetBinContent(1) - fHistPPLabIn->GetBinContent(1))/fHistP->GetEntries() - (fHistNPLabIn->GetBinContent(1) - fHistNNLabIn->GetBinContent(1))/fHistN->GetEntries();
    bfPLabIn = (fHistPNLabIn->GetBinContent(2) - fHistPPLabIn->GetBinContent(2))/fHistP->GetEntries() - (fHistNPLabIn->GetBinContent(2) - fHistNNLabIn->GetBinContent(2))/fHistN->GetEntries();
    
    bfNLabOut = (fHistPNLabOut->GetBinContent(1) - fHistPPLabOut->GetBinContent(1))/fHistP->GetEntries() - (fHistNPLabOut->GetBinContent(1) - fHistNNLabOut->GetBinContent(1))/fHistN->GetEntries();
    bfPLabOut = (fHistPNLabOut->GetBinContent(2) - fHistPPLabOut->GetBinContent(2))/fHistP->GetEntries() - (fHistNPLabOut->GetBinContent(2) - fHistNNLabOut->GetBinContent(2))/fHistN->GetEntries();
  }
  
  Double_t deltaBLabIn = bfPLabIn - bfNLabIn;
  fHistDeltaBLabIn->Fill(deltaBLabIn,gCentrality);
  
  Double_t deltaBLabOut = bfPLabOut - bfNLabOut;
  fHistDeltaBLabOut->Fill(deltaBLabOut,gCentrality);
  
  //=====================================//
  //pair's rest frame
  Double_t bfPRestIn = -999., bfNRestIn = -999.;
  Double_t bfPRestOut = -999., bfNRestOut = -999.;
  if((fHistP->GetEntries() != 0)&&(fHistN->GetEntries() != 0)) {
    bfNRestIn = (fHistPNRestIn->GetBinContent(1) - fHistPPRestIn->GetBinContent(1))/fHistP->GetEntries() - (fHistNPRestIn->GetBinContent(1) - fHistNNRestIn->GetBinContent(1))/fHistN->GetEntries();
    bfPRestIn = (fHistPNRestIn->GetBinContent(2) - fHistPPRestIn->GetBinContent(2))/fHistP->GetEntries() - (fHistNPRestIn->GetBinContent(2) - fHistNNRestIn->GetBinContent(2))/fHistN->GetEntries();
    
    bfNRestOut = (fHistPNRestOut->GetBinContent(1) - fHistPPRestOut->GetBinContent(1))/fHistP->GetEntries() - (fHistNPRestOut->GetBinContent(1) - fHistNNRestOut->GetBinContent(1))/fHistN->GetEntries();
    bfPRestOut = (fHistPNRestOut->GetBinContent(2) - fHistPPRestOut->GetBinContent(2))/fHistP->GetEntries() - (fHistNPRestOut->GetBinContent(2) - fHistNNRestOut->GetBinContent(2))/fHistN->GetEntries();
  }
  
  Double_t deltaBRestIn = bfPRestIn - bfNRestIn;
  fHistDeltaBRestIn->Fill(deltaBRestIn,gCentrality);
  
  Double_t deltaBRestOut = bfPRestOut - bfNRestOut;
  fHistDeltaBRestOut->Fill(deltaBRestOut,gCentrality);
}

//________________________________________________________________________
Double_t AliAnalysisTaskSignedBF::GetRefMultiOrCentrality(AliVEvent *event) {
  // Checks the Event cuts
  // Fills Event statistics histograms
  Double_t gMultiplicityPercentile = -1.;
  Double_t gMultiplicity = -1.;
    
  AliMultSelection *multSelection = (AliMultSelection*) event->FindListObject("MultSelection");
  if (!multSelection)
    AliFatal("MultSelection not found in input event");
    
  gMultiplicityPercentile = multSelection->GetMultiplicityPercentile(fMultiplicityEstimator, 
								     kTRUE);
  gMultiplicity = multSelection->GetEstimator(fMultiplicityEstimator)->GetValue();
  
  // error handling
  if (gMultiplicityPercentile > 100)
    gMultiplicityPercentile = -1;
  
  // filling QA histograms
  fHistMultiplicityPercentile->Fill(gMultiplicityPercentile);
  fHistMultiplicity->Fill(gMultiplicity);
  fHistMultiplicityvsPercentile->Fill(gMultiplicity, gMultiplicityPercentile);

  // decide what should be returned only here
  Double_t lReturnVal = -100;
  lReturnVal = gMultiplicityPercentile;
  
  return lReturnVal;
}

//________________________________________________________________________
Double_t AliAnalysisTaskSignedBF::IsEventAccepted(AliVEvent *event){
  // Checks the Event cuts
  // Fills Event statistics histograms
  Bool_t isSelectedMain = kTRUE;
  Float_t gRefMultiplicity = -1.;

  //AliMCEvent *mcevent = dynamic_cast<AliMCEvent*>(event);  
  fHistEventStats->Fill(1); //all events

  if(fUseOfflineTrigger)
    isSelectedMain = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  
  if(isSelectedMain) {
    fHistEventStats->Fill(2); //triggered events
 
    const AliVVertex *vertex = event->GetPrimaryVertex();  
    if(vertex) {
      Double32_t fCov[6];
      vertex->GetCovarianceMatrix(fCov);
      if(vertex->GetNContributors() > 0) {
	if(fCov[5] != 0) {
	  fHistEventStats->Fill(3); //proper vertex
	  if(TMath::Abs(vertex->GetX()) < fVxMax) {
	    if(TMath::Abs(vertex->GetY()) < fVyMax) {
	      if( TMath::Abs(vertex->GetZ()) < fVzMax) {
		
		fHistEventStats->Fill(4);//analyzed events
		
		if (fUseAdditionalVtxCuts){
		  if(vertex->GetNContributors()<1) {
		    Printf("No Primary Vertex");
		    return -1;
		  }
		  
		  const AliVVertex *vSPD = ((AliAODEvent*)event)->GetPrimaryVertexSPD();
		  if (!vSPD) {
		    Printf("No vertex SPD");
		    return -1;
		  }
		  Double_t dz = vSPD->GetZ()- vertex->GetZ();
		  double covTrc[6],covSPD[6];
		  vertex->GetCovarianceMatrix(covTrc);
		  vSPD->GetCovarianceMatrix(covSPD);
		  double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
		  double errTrc = TMath::Sqrt(covTrc[5]);
		  double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
		  if (TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrc>20)
		    return -1;
		     
		  fHistEventStats->Fill(5); 
		}//additional vertex cuts
		
		// check for pile-up event
		if(fCheckPileUp){
		  fUtils->SetUseMVPlpSelection(kTRUE);
		  if(fUtils->IsPileUpEvent(event))
		    return -1.;
		  fHistEventStats->Fill(6); 
		}
			        
		gRefMultiplicity = GetRefMultiOrCentrality(event);
		
		fHistVx->Fill(vertex->GetX());
		fHistVy->Fill(vertex->GetY());
		fHistVz->Fill(vertex->GetZ());
		
				   
		// take only events inside centrality class
		// if(fUseCentrality) {
		if((gRefMultiplicity > fCentralityPercentileMin) && (gRefMultiplicity < fCentralityPercentileMax)){
		  fHistEventStats->Fill(7); //events with correct centrality

		  return gRefMultiplicity;		
		}//centrality class
	      }//Vz cut
	    }//Vy cut
	  }//Vx cut
	}//proper vertex resolution
      }//proper number of contributors
    }//vertex object valid
  }//triggered event 
    
  // in all other cases return -1 (event not accepted)
  return -1;
}

//________________________________________________________________________
Double_t AliAnalysisTaskSignedBF::GetEventPlane(AliVEvent *event){
  // Get the event plane
  Float_t gVZEROEventPlane    = -10.;
  Float_t gReactionPlane      = -10.;
  Double_t qxTot = 0.0, qyTot = 0.0;
  Int_t gHarmonic = 2;

  //MC: from reaction plane
  /*if(fAnalysisLevel == "MC"){
    if(!event) {
      AliError("mcEvent not available");
      return 0x0;
    }

    AliMCEvent *gMCEvent = dynamic_cast<AliMCEvent*>(event);
    if(gMCEvent){
      TString genName;
      AliCollisionGeometry* headerH;
      TList *ltgen = (TList*)gMCEvent->GetCocktailList();
      if (ltgen) {
	for(auto&& listObject: *ltgen){
	    genName = Form("%s",listObject->GetName());
	    if (genName.Contains("Hijing")) {
		headerH = dynamic_cast<AliCollisionGeometry*>(listObject);
		break;
	      }
	  }
      }
      else 
	headerH = dynamic_cast<AliCollisionGeometry*>(gMCEvent->GenEventHeader());  
      if (headerH) {
	gReactionPlane = headerH->ReactionPlaneAngle();
	//gReactionPlane *= TMath::RadToDeg();
      }//MC header
    }//MC event cast
    }*///MC
  // AOD,ESD,ESDMC: from VZERO Event Plane
  //else {
  const AliQnCorrectionsQnVector *gQnVector;
  Double_t gEventPlane = 0.0;
  /* get the fully corrected Qn vector from VZEROA sub-detector */
  gQnVector = fFlowQnVectorMgr->GetDetectorQnVector(fEventPlaneDetector.Data());
  if (gQnVector != NULL)
    gEventPlane = gQnVector->EventPlane(gHarmonic);

  gReactionPlane = gVZEROEventPlane;
  //}//AOD,ESD,ESDMC
  
  return gReactionPlane;
}
