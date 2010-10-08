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
#include "AliAODJet.h"

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
  fHistPtJets(NULL),
  fHistGammaJets(NULL),
  fHistGammaJetsIso(NULL),
  fMinPt(5.0),
  fConeSize(0.3),
  fPtThreshold(2.0),
  fDeltaAODFileName(""),
  fPhotons(NULL)
{
  // Dummy Constructor
}


//________________________________________________________________________
AliAnalysisTaskGammaJet::AliAnalysisTaskGammaJet(const char *name) : 
  AliAnalysisTaskSE(name), 
  fOutputList(0), 
  fHistPt(0),
  fHistPtPhos(0),
  fHistPtEmcal(0),
  fHistPtJets(0),
  fHistGammaJets(NULL),
  fHistGammaJetsIso(NULL),
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


  fHistPtJets = new TH1F("fHistPtJets", "P_{T} distribution", 150, 0.1, 50);
  fHistPtJets->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPtJets->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPtJets->SetMarkerStyle(kFullCircle);
  fOutputList->Add(fHistPtJets);

  fHistGammaJets = new TH1F("fHistGammaJets", "fHistGammaJets", 200, -TMath::Pi(), 2*TMath::Pi());
  fOutputList->Add(fHistGammaJets);
  
  fHistGammaJetsIso = new TH1F("fHistGammaJetsIso", "fHistGammaJetsIso", 200, -TMath::Pi(), 2*TMath::Pi());
  fOutputList->Add(fHistGammaJetsIso);
  
  //TNtuple * tuple = new TNtuple("fNtuple", "fNtuple", dPhi, 


  ///Create AOD branch
  fPhotons = new TClonesArray("AliAODPWG4ParticleCorrelation", 0);
  fPhotons->SetName("fPhotons");
  AddAODBranch("TClonesArray", &fPhotons);

  ///Isolation class
  // fIsolation = new AliAnaParticleIsolation();
  // fIsolation->SetInputAODName("fPhotons");


}

//________________________________________________________________________
void AliAnalysisTaskGammaJet::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  

  //Clear stuff for new event
  CleanUp();


  ///Get AOD event
  AliAODEvent * aodEvent = GetAODEvent();
  if(!aodEvent) {
    AliError("No AOD event!!");
    return;
  }
  
  //FillPWG4PartCorrBranch(convGamma, fPhotons, "ConvGamma");
  //fIsolation->MakeAnalysisFillAOD();
  
  ProcessConvGamma(aodEvent);
  ProcessCalorimeters(aodEvent);
    

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


//_____________________________________________________________________
void AliAnalysisTaskGammaJet::FillPWG4PartCorrBranch( TClonesArray * gcBranch, TClonesArray * partCorrBranch , TString detector ) {
  
  for(int i = 0; i < gcBranch->GetEntriesFast(); i++) {
    AliGammaConversionAODObject * gcObject = dynamic_cast<AliGammaConversionAODObject*>(gcBranch->At(i));
    if ( gcObject ) {
      AliAODPWG4ParticleCorrelation pc(gcObject->Px(), gcObject->Py(), gcObject->Pz(), gcObject->E()); 
      pc.SetTagged(gcObject->IsTagged());
      pc.SetTrackLabel(gcObject->GetLabel1(), gcObject->GetLabel2());
      pc.SetDetector(detector);
      new((*partCorrBranch)[i]) AliAODPWG4ParticleCorrelation(pc);
    
    } else {
      AliError(Form("Couldn't get gamma conversion aod object"));
    }
   
  }
}


//_____________________________________________________________________
AliAODPWG4ParticleCorrelation * AliAnalysisTaskGammaJet::PWG4PartFromGammaConvAODObject(AliGammaConversionAODObject * gcObject, TString detector ) {

  AliAODPWG4ParticleCorrelation * pc = new AliAODPWG4ParticleCorrelation(gcObject->Px(), gcObject->Py(), gcObject->Pz(), gcObject->E());
  pc->SetTagged(gcObject->IsTagged());
  pc->SetTrackLabel(gcObject->GetLabel1(), gcObject->GetLabel2());
  pc->SetDetector(detector);
  return pc;
}


//_________________________________________________________________________
void AliAnalysisTaskGammaJet::CleanUp() {
  fPhotons->Delete();
}

//_________________________________________________________________________
Bool_t AliAnalysisTaskGammaJet::IsIsolated( AliAODPWG4ParticleCorrelation * particle, TClonesArray * tracks, Float_t coneSize, Float_t ptThreshold ) {
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
    
    AliAODPWG4ParticleCorrelation * photon = new AliAODPWG4ParticleCorrelation(tlvec);
    
    photon->SetIsolated( IsIsolated(photon, aodEvent->GetTracks(), GetConeSize(), GetPtThreshold()) );
    CorrelateWithJets(photon, aodEvent->GetJets());
  }
  
}
//___________________________________________________________________________________________
void AliAnalysisTaskGammaJet::ProcessConvGamma( const AliAODEvent * const aodEvent ) {

  TClonesArray * convGamma = GetConversionGammas(aodEvent);
  if(!convGamma) {
    AliError(Form("No convgamma"));
    return;
  }

  for (Int_t iPhot = 0; iPhot < convGamma->GetEntriesFast(); iPhot++) {
    AliGammaConversionAODObject * aodO = dynamic_cast<AliGammaConversionAODObject*>(convGamma->At(iPhot));
    
    if (!aodO) {
      AliError(Form("ERROR: Could not receive ga %d\n", iPhot));
      continue;
    }
    
    if(aodO->Pt() < GetMinPt()) continue;
    

    //Use the AODPWG4PartCorr shit!
    AliAODPWG4ParticleCorrelation * photon = PWG4PartFromGammaConvAODObject(aodO, "ConvGamma");
    photon->SetIsolated( IsIsolated(photon, aodEvent->GetTracks(), GetConeSize(), GetPtThreshold()) );
    
    
    // if ( (aodO->Phi()) < 0 )
    //   cout << aodO->Phi() << endl;
    
    CorrelateWithJets(photon, aodEvent->GetJets());

    fHistPt->Fill(photon->Pt());
    delete photon;
    
  }
}

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
