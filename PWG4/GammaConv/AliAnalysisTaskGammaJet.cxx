#include <iostream>
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskGammaJet.h"

#include "AliESDEvent.h"
#include "AliESDCaloCluster.h"
#include "AliESDInputHandler.h"

#include "AliAODPWG4ParticleCorrelation.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODCaloCluster.h"
#include "AliGammaConversionAODObject.h"
#include "AliAODConversionParticle.h"
#include "AliAODJet.h"

#include "AliAODInputHandler.h"

// Gamma - jet correlation analysis task
// Authors: Svein Lindal


using namespace std;

ClassImp(AliAnalysisTaskGammaJet)

//________________________________________________________________________
AliAnalysisTaskGammaJet::AliAnalysisTaskGammaJet() : AliAnalysisTaskSE(), 
  fOutputList(NULL), 
  fHistPt(NULL),
  fHistPtPhos(NULL),
  fHistPtEmcal(NULL),

  fHistPhotPhi(NULL), 
  fHistHadPhi(NULL), 
  fHistJetPhi(NULL), 

  fHistPtJets(NULL),
  fHistGammaJets(NULL),
  fHistGammaJetsIso(NULL),
  fHistMaxdPhi(NULL),
  fHistMaxdPhiIso(NULL),
  fHistMaxdPhiIsoPt(NULL),

  fHadHistPt(NULL), //! Pt spectrum
  fHadHistdPhi(NULL), //!Phi correlations
  fHadHistdPhiIso(NULL), //!Phi correlations
  fHadHistMaxdPhi(NULL), //!Phi correlations
  fHadHistMaxdPhiIso(NULL), //!Phi correlations
  fHadHistMaxdPhiIsoPt(NULL), //!Phi correlations

  fMinPt(2.0),
  fConeSize(0.9),
  fPtThreshold(2.0),
  fDeltaAODFileName(""),
  fPhotons(NULL)
{
  // Dummy Constructor
}

//________________________________________________________________________________
AliAnalysisTaskGammaJet::~AliAnalysisTaskGammaJet() {

  if(fOutputList)
    delete fOutputList; 
  fOutputList = NULL;

  if(fHistPt)
    fHistPt = NULL;
  delete fHistPt;
 
  if(fHistPtPhos)
    fHistPtPhos = NULL;
  delete fHistPtPhos;

  if(fHistPtEmcal)
    fHistPtEmcal = NULL;
  delete fHistPtEmcal;

  if(fHistPtJets)
    fHistPtJets= NULL;
  delete fHistPtJets;

    if(fHistGammaJets)
      fHistGammaJets = NULL;
  delete fHistGammaJets;

  if(fHistGammaJetsIso)
    fHistGammaJetsIso = NULL;
  delete fHistGammaJetsIso;
 
  if(fHistMaxdPhi)
    fHistMaxdPhi = NULL;
  delete fHistMaxdPhi;
 
  if(fHistMaxdPhiIso)
    fHistMaxdPhiIso = NULL;
  delete fHistMaxdPhiIso;
 
  if(fHistMaxdPhiIsoPt)
    fHistMaxdPhiIsoPt = NULL;
  delete fHistMaxdPhiIsoPt;
   
  if(fPhotons)
    fPhotons = NULL;
  delete fPhotons;


}



//________________________________________________________________________
AliAnalysisTaskGammaJet::AliAnalysisTaskGammaJet(const char *name) : 
  AliAnalysisTaskSE(name), 
  fOutputList(0), 
  fHistPt(0),
  fHistPtPhos(0),
  fHistPtEmcal(0),

  fHistPhotPhi(NULL), 
  fHistHadPhi(NULL), 
  fHistJetPhi(NULL), 

  fHistPtJets(0),
  fHistGammaJets(NULL),
  fHistGammaJetsIso(NULL),
  fHistMaxdPhi(NULL),
  fHistMaxdPhiIso(NULL),
  fHistMaxdPhiIsoPt(NULL),

  fHadHistPt(NULL), //! Pt spectrum
  fHadHistdPhi(NULL), //!Phi correlations
  fHadHistdPhiIso(NULL), //!Phi correlations
  fHadHistMaxdPhi(NULL), //!Phi correlations
  fHadHistMaxdPhiIso(NULL), //!Phi correlations
  fHadHistMaxdPhiIsoPt(NULL), //!Phi correlations

  fMinPt(0.0),
  fConeSize(0.0),
  fPtThreshold(0.0),
  fDeltaAODFileName(""),
  fPhotons(NULL)
{
  // Constructor
  // Define input and output slots here
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD

  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskGammaJet::UserCreateOutputObjects()
{
  //Create histograms add, to outputlist
  fOutputList = new TList();

  fHistPt = new TH1F("fHistPt", "P_{T} distribution", 150, 0.1, 50);
  fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPt->SetMarkerStyle(kFullCircle);
  fOutputList->Add(fHistPt);
  
  fHistPtPhos = new TH1F("fHistPtPhos", "P_{T} distribution", 150, 0.1, 50);
  fHistPtPhos->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPtPhos->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPtPhos->SetMarkerStyle(kFullCircle);
  fOutputList->Add(fHistPtPhos);
  
  fHistPtEmcal = new TH1F("fHistPtEmcal", "P_{T} distribution", 150, 0.1, 50);
  fHistPtEmcal->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPtEmcal->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPtEmcal->SetMarkerStyle(kFullCircle);
  fOutputList->Add(fHistPtEmcal);

  fHistPhotPhi = new TH1F("phi gamma", "phi gamma", 120, -6.3, 6.3);
  fOutputList->Add(fHistPhotPhi);
  fHistHadPhi = new TH1F("phi track", "phi track", 120, -6.3, 6.3);
  fOutputList->Add(fHistHadPhi);
  fHistJetPhi = new TH1F("phi jet", "phi jet", 120, -6.3, 6.3);
  fOutputList->Add(fHistJetPhi);


  fHistPtJets = new TH1F("fHistPtJets", "P_{T} distribution", 150, 0.1, 50);
  fHistPtJets->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPtJets->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPtJets->SetMarkerStyle(kFullCircle);
  fOutputList->Add(fHistPtJets);

  fHistGammaJets = new TH1F("fHistGammaJets", "fHistGammaJets", 200, -2*TMath::Pi(), 2*TMath::Pi());
  fOutputList->Add(fHistGammaJets);
  
  fHistGammaJetsIso = new TH1F("fHistGammaJetsIso", "fHistGammaJetsIso", 200, -2*TMath::Pi(), 2*TMath::Pi());
  fOutputList->Add(fHistGammaJetsIso);


  fHistMaxdPhi = new TH1F("fHistMaxdPhi", "fHistMaxdPhi", 200, -2*TMath::Pi(), 2*TMath::Pi());
  fOutputList->Add(fHistMaxdPhi);
  
  fHistMaxdPhiIso = new TH1F("fHistMaxdPhiIso", "fHistMaxdPhiIso", 200, -2*TMath::Pi(), 2*TMath::Pi());
  fOutputList->Add(fHistMaxdPhiIso);

  fHistMaxdPhiIsoPt = new TH1F("fHistMaxdPhiIsoPt", "fHistMaxdPhiIsoPt", 200, -2*TMath::Pi(), 2*TMath::Pi());
  fOutputList->Add(fHistMaxdPhiIsoPt);




  fHadHistPt = new TH1F("fHadHistPt", "P_{T} distribution", 150, 0.1, 50);
  fHadHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHadHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHadHistPt->SetMarkerStyle(kFullCircle);
  fOutputList->Add(fHadHistPt);
  

  fHadHistdPhi = new TH1F("fHadHistdPhi", "fHadHistdPhi", 200, -2*TMath::Pi(), 2*TMath::Pi());
  fOutputList->Add(fHadHistdPhi);
  
  fHadHistdPhiIso = new TH1F("fHadHistdPhiIso", "fHadHistdPhiIso", 200, -2*TMath::Pi(), 2*TMath::Pi());
  fOutputList->Add(fHadHistdPhiIso);


  fHadHistMaxdPhi = new TH1F("fHadHistMaxdPhi", "fHadHistMaxdPhi", 200, -2*TMath::Pi(), 2*TMath::Pi());
  fOutputList->Add(fHadHistMaxdPhi);
  
  fHadHistMaxdPhiIso = new TH1F("fHadHistMaxdPhiIso", "fHadHistMaxdPhiIso", 200, -2*TMath::Pi(), 2*TMath::Pi());
  fOutputList->Add(fHadHistMaxdPhiIso);

  fHadHistMaxdPhiIsoPt = new TH1F("fHadHistMaxdPhiIsoPt", "fHadHistMaxdPhiIsoPt", 200, -2*TMath::Pi(), 2*TMath::Pi());
  fOutputList->Add(fHadHistMaxdPhiIsoPt);




  
  //TNtuple * tuple = new TNtuple("fNtuple", "fNtuple", dPhi, 


  ///Create AOD branch
  fPhotons = new TClonesArray("AliAODPWG4ParticleCorrelation", 0);
  fPhotons->SetName("ConversionGamma");
  AddAODBranch("TClonesArray", &fPhotons);


}

//________________________________________________________________________
void AliAnalysisTaskGammaJet::UserExec(Option_t *) 
{


  //BALLE BALLE not do always
  AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kTRUE);


  //Clear stuff for new event
  CleanUp();


  


  ///Get AOD event
  AliAODEvent * aodEvent = GetAODEvent();
  if(!aodEvent) {
    AliError("No AOD event!!");
    return;
  }
  
  ProcessConvGamma(aodEvent);
  //ProcessCalorimeters(aodEvent);
    

  PostData(1, fOutputList);
        
}
//_____________________________________________________________________
void AliAnalysisTaskGammaJet::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query
}

//_____________________________________________________________________
AliAODEvent * AliAnalysisTaskGammaJet::GetAODEvent() {
  //Get the AOD event from whereever it might be
  AliAODEvent * aodEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!aodEvent) {
    aodEvent = AODEvent();
  }
  
  return aodEvent;

}


//_____________________________________________________________________
TClonesArray * AliAnalysisTaskGammaJet::GetConversionGammas(const AliAODEvent * aodEvent) {

  //Get Conversion gamma branch of AOD. First try standard AOD
  TClonesArray * convGamma = dynamic_cast<TClonesArray*>(aodEvent->FindListObject("GammaConv_gamma"));
  
  //If it's there, send it back
  if(convGamma)  return convGamma;


  //If AOD not in standard file have to locate it in delta AOD
  if( !(fDeltaAODFileName.Length() > 0)  ) return NULL;
  
  AliAODHandler * aodHandler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()); 
  if(aodHandler) {
    AliAODExtension * gExt = dynamic_cast<AliAODExtension*>(aodHandler->GetExtensions()->FindObject(fDeltaAODFileName));
    if(gExt) {
      AliAODEvent * gcEvent = gExt->GetAOD();
      return dynamic_cast<TClonesArray*>(gcEvent->FindListObject("GammaConv_gamma"));
    }
  }  
  return NULL;
}


// //_____________________________________________________________________
// void AliAnalysisTaskGammaJet::FillPWG4PartCorrBranch( TClonesArray * gcBranch, TClonesArray * partCorrBranch , TString detector ) {
  
//   for(int i = 0; i < gcBranch->GetEntriesFast(); i++) {
//     AliGammaConversionAODObject * gcObject = dynamic_cast<AliGammaConversionAODObject*>(gcBranch->At(i));
//     if ( gcObject ) {
//       AliAODPWG4ParticleCorrelation pc(gcObject->Px(), gcObject->Py(), gcObject->Pz(), gcObject->E()); 
//       pc.SetTagged(gcObject->IsTagged());
//       pc.SetTrackLabel(gcObject->GetLabel1(), gcObject->GetLabel2());
//       pc.SetDetector(detector);
//       new((*partCorrBranch)[i]) AliAODPWG4ParticleCorrelation(pc);
    
//     } else {
//       AliError(Form("Couldn't get gamma conversion aod object"));
//     }
   
//   }
// }


// //_____________________________________________________________________
// AliAODPWG4ParticleCorrelation * AliAnalysisTaskGammaJet::PWG4PartFromGammaConvAODObject(AliGammaConversionAODObject * gcObject, TString detector ) {

//   AliAODPWG4ParticleCorrelation * pc = new AliAODPWG4ParticleCorrelation(gcObject->Px(), gcObject->Py(), gcObject->Pz(), gcObject->E());
//   pc->SetTagged(gcObject->IsTagged());
//   pc->SetTrackLabel(gcObject->GetLabel1(), gcObject->GetLabel2());
//   pc->SetDetector(detector);
//   return pc;
// }


//_________________________________________________________________________
void AliAnalysisTaskGammaJet::CleanUp() {
  fPhotons->Delete();
}

//_________________________________________________________________________
Bool_t AliAnalysisTaskGammaJet::IsIsolated( AliAODPWG4Particle * particle, TClonesArray * tracks, Float_t coneSize, Float_t ptThreshold ) {
  //See header file for documentation
  for(int it = 0; it < tracks->GetEntriesFast(); it++) {
    if ( (it == particle->GetTrackLabel(0)) || it == particle->GetTrackLabel(1) ) 
      continue;

    //BALLE Svein:How are you checking the calorimeters for whether they are decay particles ?

    AliAODTrack * track = dynamic_cast<AliAODTrack*>(tracks->At(it));
    if (track) {
      if ( IsInCone(particle->Eta() - track->Eta(), particle->Phi() - track->Phi(), coneSize) ) {
	if (track->Pt() > ptThreshold) {
	  return kFALSE;
	}
      }
    } else {
      AliError(Form("Bad track!!!! "));
    }
  }
  
  //No particle above threshold, it's isolated
  return kTRUE;
}


//______________________________________________________________________________________________
void AliAnalysisTaskGammaJet::ProcessCalorimeters( const AliAODEvent * const aodEvent ) {
  
  TClonesArray * clusters = aodEvent->GetCaloClusters();
  

  for(int ic = 0; ic < clusters->GetEntriesFast(); ic++) {
    AliAODCaloCluster * cluster = dynamic_cast<AliAODCaloCluster*>(clusters->At(ic));
    if (!cluster) { 
      AliError(Form("Error getting cluster"));
      continue;
    }


    if (cluster->GetNCells() < 6) continue;
    if (cluster->GetEmcCpvDistance() < 15) continue;

    TLorentzVector tlvec;
    
    AliAODVertex * vertex = aodEvent->GetPrimaryVertex();
    Double_t vertexPosition[3];
    vertex->GetXYZ(vertexPosition);
    cluster->GetMomentum(tlvec, vertexPosition);
    if (tlvec.Pt() < GetMinPt()) continue; 
    
  }
  
}

// // ///____________________________________________________________________________________
// // void AddToAOD(AliAODConversionParticle * photon, TClonesArray * branch) {


// // }


// //__________________________________________________________________________________________
// void AddToOutputAOD(AliGammaConversionAODObject * photon, TClonesArray * branch) {
  
//   cout <<"BALLE BALLE BALLE AddToOutputAOD"<<endl;

//   AliAODPWG4ParticleCorrelation particle = AliAODPWG4ParticleCorrelation(photon->Px(), photon->Py(), photon->Pz(), photon->E());
//   particle.SetTrackLabel(photon->GetLabel1(), photon->GetLabel2());
//   particle.SetDetector("CTS");
//   //particle.SetPdg(AliPID::kElecon);

//   Int_t i = branch->GetEntriesFast();
//   // if(! (branch.GetClass("Correlation")) ) {
//   //   new((*branch)[i])  AliAODPWG4Particle(particle);
//   // } else {
//   TList * list = particle.GetObjArrayList();
//   if(list)
//     cout <<"BALLE BALLE we have the list"<<endl;
//   else
//     cout <<"BALLE BALLE we don't"<<endl;

//  particle.GetObjArrayList()->Dump();
//   new((*branch)[i])  AliAODPWG4ParticleCorrelation(particle);
//     //}
// }


//___________________________________________________________________________________________
void AliAnalysisTaskGammaJet::ProcessConvGamma( const AliAODEvent * const aodEvent ) {
  
  TClonesArray * tracks = aodEvent->GetTracks();
  if(!tracks) {
    cout << "No tracks!!!"<<endl;
    return;
  }

 
  TClonesArray * convGamma = dynamic_cast<TClonesArray*>(aodEvent->FindListObject("GammaConversionTask_900356200010031"));
  

  if(!convGamma) {

    convGamma = GetConversionGammas(aodEvent);
    if(!convGamma) {
      AliError(Form("No convgamma"));
      return;
    }
  }
  


  for (Int_t iPhot = 0; iPhot < convGamma->GetEntriesFast(); iPhot++) {
    AliAODPWG4Particle * photon = dynamic_cast<AliAODPWG4Particle*>(convGamma->At(iPhot));
    
    
    if(!photon) {
      AliGammaConversionAODObject * aodO = dynamic_cast<AliGammaConversionAODObject*>(convGamma->At(iPhot));
      if (!aodO) {
	AliError(Form("ERROR: Could not receive ga %d\n", iPhot));
	continue;
      }
      
      if(aodO->Pt() < GetMinPt()) continue;
      photon = AddToAOD(aodO, fPhotons, "ConvGamma");
    } 
  
    if(photon) {
      Bool_t isolated = IsIsolated(photon, tracks, GetConeSize(), GetPtThreshold() );
      
      
      CorrelateWithJets(photon, aodEvent->GetJets(), isolated);
      CorrelateWithHadrons(photon, aodEvent->GetTracks(), isolated);

      fHistPt->Fill(photon->Pt());
      fHistPhotPhi->Fill(photon->Phi());
    }
  }
}


AliAODPWG4ParticleCorrelation * AliAnalysisTaskGammaJet::AddToAOD(AliGammaConversionAODObject * aodO, TClonesArray * branch, TString detector) {
  new((*branch)[branch->GetEntriesFast()]) AliAODPWG4ParticleCorrelation(aodO->Px(), aodO->Py(), aodO->Pz(), aodO->E());
  AliAODPWG4ParticleCorrelation * photon = dynamic_cast<AliAODPWG4ParticleCorrelation*>(branch->Last());
  if(photon) {
    photon->SetTagged(aodO->IsTagged());
    photon->SetTrackLabel(aodO->GetLabel1(), aodO->GetLabel2());
    photon->SetDetector(detector);
  }
  
  return photon;
}

AliAODPWG4ParticleCorrelation * AliAnalysisTaskGammaJet::AddToAOD(AliAODConversionParticle * aodO, TClonesArray * branch, TString detector) {
  new((*branch)[branch->GetEntriesFast()]) AliAODPWG4ParticleCorrelation(aodO->Px(), aodO->Py(), aodO->Pz(), aodO->E());
  AliAODPWG4ParticleCorrelation * photon = dynamic_cast<AliAODPWG4ParticleCorrelation*>(branch->Last());
  if(photon) {
    photon->SetTrackLabel(aodO->GetLabel1(), aodO->GetLabel2());
    photon->SetDetector(detector);
  }

  return photon;
}

///____________________________________________________________________________________________________
void AliAnalysisTaskGammaJet::CorrelateWithHadrons(AliAODPWG4Particle * photon, const TClonesArray * tracks, Bool_t const isolated) {

  //See header file for documentation
  if (tracks) {
 
    Float_t maxdPhi = 0.0;
    Float_t maxdPhiPt = 0.0;

    for(int ij = 0; ij < tracks->GetEntriesFast(); ij++) {
      AliAODTrack * track = dynamic_cast<AliAODTrack*>(tracks->At(ij));
      if(track) {
	if (track->Pt() < 2.0 ) continue;

	fHadHistPt->Fill(track->Pt());
	fHistHadPhi->Fill(track->Phi());
	Float_t dPhi = TMath::Abs(photon->Phi() - track->Phi());
	
	fHadHistdPhi->Fill(dPhi);
	if (isolated) {
	  fHadHistdPhiIso->Fill(dPhi, track->Pt()/photon->Pt());
	}


	if(photon->Phi() < 0)
	  cout << "BALLE"<<endl;
	if(track->Phi() < 0 )
	  cout << "KUKKK"<<endl;
	  
	//cout << dPhi << " " << maxdPhi << endl;
	//cout << TMath::Abs( dPhi - TMath::Pi() ) << " " << TMath::Abs( maxdPhi - TMath::Pi()) <<endl;;

	if ( TMath::Abs( dPhi - TMath::Pi() ) <   TMath::Abs( maxdPhi - TMath::Pi()) ){
	  maxdPhi = dPhi;

	  maxdPhiPt = track->Pt();
	  //cout << dPhi << " " << maxdPhi << endl;
	
	}

      }
    }


    if(tracks->GetEntriesFast() > 0) {
      //cout << maxdPhi << endl;
      fHadHistMaxdPhi->Fill(maxdPhi);
      if(isolated) {
	fHadHistMaxdPhiIso->Fill(maxdPhi);
	fHadHistMaxdPhiIsoPt->Fill(maxdPhi, maxdPhiPt/photon->Pt());
      }
    }
  }



}


///________________________________________________________________________________________________________________
void AliAnalysisTaskGammaJet::CorrelateWithJets(AliAODPWG4ParticleCorrelation * photon, const TClonesArray * const jets) {
  //See header file for documentation
  if (jets) {
    for(int ij = 0; ij < jets->GetEntriesFast(); ij++) {
      AliAODJet * jet = dynamic_cast<AliAODJet*>(jets->At(ij));
      if(jet) {
	fHistPtJets->Fill(jet->Pt());
	
	Float_t dPhi = TMath::Abs(photon->Phi() - jet->Phi());
	if (photon->IsIsolated())
	  fHistGammaJetsIso->Fill(dPhi, jet->Pt()/photon->Pt());
	else
	  fHistGammaJets->Fill(dPhi);
	
      }
    }
  }
}



///________________________________________________________________________________________________________________
void AliAnalysisTaskGammaJet::CorrelateWithJets(AliAODPWG4Particle * photon, const TClonesArray * const jets, Bool_t const isolated ) {
  //See header file for documentation
  if (jets) {
 
    Float_t maxdPhi = 0.0;
    Float_t maxdPhiPt = 0.0;

    for(int ij = 0; ij < jets->GetEntriesFast(); ij++) {
      AliAODJet * jet = dynamic_cast<AliAODJet*>(jets->At(ij));
      if(jet) {
	fHistPtJets->Fill(jet->Pt());
	fHistJetPhi->Fill(jet->Phi());
	Float_t dPhi = TMath::Abs(photon->Phi() - jet->Phi());

	fHistGammaJets->Fill(dPhi);
	if (isolated) {
	  fHistGammaJetsIso->Fill(dPhi, jet->Pt()/photon->Pt());
	}



	if ( TMath::Abs( dPhi - TMath::Pi() ) <   TMath::Abs( maxdPhi - TMath::Pi()) ){
	  maxdPhi = dPhi;

	  maxdPhiPt = jet->Pt();
	  //cout << dPhi << " " << maxdPhi << endl;
	
	}
	 
      }
    }


    if(jets->GetEntriesFast() > 0) {
      //cout << maxdPhi << endl;
      fHistMaxdPhi->Fill(maxdPhi);
      if(isolated) {
	fHistMaxdPhiIso->Fill(maxdPhi);
	fHistMaxdPhiIsoPt->Fill(maxdPhi, maxdPhiPt/photon->Pt());
      }
    }
  }
}
