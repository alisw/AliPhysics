
/* $Id$ */

#include "TChain.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <TLorentzVector.h>
#include <TVector3.h>

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAnalysisTaskAODvsESD.h"

ClassImp(AliAnalysisTaskAODvsESD)

//________________________________________________________________________
AliAnalysisTaskAODvsESD::AliAnalysisTaskAODvsESD(const char *name) 
: AliAnalysisTask(name, ""), fESD(0), fAOD(0), fList(0), fMuonNtuple(0), fMuonNtupleAOD(0), fInvMass(0), fInvMassAOD(0)
{
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes NTuple/histos into a TList
  DefineOutput(0, TList::Class());  
}

//________________________________________________________________________
void AliAnalysisTaskAODvsESD::ConnectInputData(Option_t *) 
{
  // Connect ESD here
  TTree* esdTree = dynamic_cast<TTree*> (GetInputData(0));
  if (!esdTree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());   
    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else
      fESD = esdH->GetEvent();
  }  

  // Connect AOD here
  AliAODHandler *aodH = dynamic_cast<AliAODHandler*> (AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());  
  if (!aodH) {
    Printf("ERROR: Could not get AODHandler");
  } else
    fAOD = aodH->GetAOD();
}

//________________________________________________________________________
void AliAnalysisTaskAODvsESD::CreateOutputObjects()
{
  // This method has to be called INSIDE the user redefined CreateOutputObjects
  // method, before creating each object corresponding to the output containers
  // that are to be written to a file. This need to be done in general for the big output
  // objects that may not fit memory during processing. 
  OpenFile(0); 

  // Create Ntuples 
  // For single muons: thetaX, thetaY, ptInv, eta, phi, theta, px, py, pz,ch
  fMuonNtuple = new TNtuple("fMuonNtuple","Muon information","thX:thY:ptI:eta:phi:theta:px:py:pz");
  fMuonNtupleAOD = new TNtuple("fMuonNtupleAOD","Muon information","etaAOD:phiAOD:thetaAOD:pxAOD:pyAOD:pzAOD");
  // Create histos for inv mass
  fInvMass = new TH1F("fInvMass","Inv. mass from ESDs",140,0,7);
  fInvMassAOD = new TH1F("fInvMassAOD","Inv. mass from AOD",140,0,7);

  // Add Ntuples to the list
  fList = new TList();
  fList->Add(fMuonNtuple);
  fList->Add(fMuonNtupleAOD);
  fList->Add(fInvMass);
  fList->Add(fInvMassAOD);
}

//________________________________________________________________________
void AliAnalysisTaskAODvsESD::Exec(Option_t *) 
{
  // Main loop, called for each event
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }
  if (!fAOD) {
    Printf("ERROR: fAOD not available");
    return;
  }
  
  // ESDs
  // for muons
  Float_t muonMass = 0.105658369;
  Float_t thetaX1=0;    Float_t thetaY1=0;
  Float_t ptInv1=0;     Float_t pYZ1=0;
  Float_t pxM1=0;       Float_t pyM1=0;        Float_t pzM1=0;
  Float_t etaM1=0;      Float_t phiM1=0;       Float_t thetaM1=0;
  Int_t chargeM1=0;      Float_t energyM1=0;    TLorentzVector lvM1;
  Float_t thetaX2=0;    Float_t thetaY2=0;
  Float_t ptInv2=0;     Float_t pYZ2=0;
  Float_t pxM2=0;       Float_t pyM2=0;        Float_t pzM2=0;
  Int_t chargeM2=0;      Float_t energyM2=0;    TLorentzVector lvM2;
  //for J/psi
  Float_t invM=0;
  TLorentzVector lvJpsi;
 
  int nmt = fESD->GetNumberOfMuonTracks();
  // first loop over muon tracks in event
  for (Int_t mTrack = 0; mTrack < nmt; mTrack++) {
    thetaX1=0; thetaY1=0; ptInv1=0; pYZ1=0; chargeM1=0;
    pxM1=0;    pyM1=0;    pzM1=0;   energyM1=0;
    AliESDMuonTrack* muonTrack1 = fESD->GetMuonTrack(mTrack);
    thetaX1 = muonTrack1->GetThetaX();
    thetaY1 = muonTrack1->GetThetaY();
    ptInv1 = TMath::Abs(muonTrack1->GetInverseBendingMomentum());
    pYZ1 = 1/ptInv1;
    pzM1  = - pYZ1 / TMath::Sqrt(1.0 + TMath::Tan(thetaY1)*TMath::Tan(thetaY1));
    pxM1  = pzM1 * TMath::Tan(thetaX1);
    pyM1  = pzM1 * TMath::Tan(thetaY1);
    energyM1 = TMath::Sqrt(muonMass*muonMass + pxM1*pxM1 + pyM1*pyM1 + pzM1*pzM1);
    lvM1.SetPxPyPzE(pxM1,pyM1,pzM1,energyM1);
    chargeM1 = Int_t(TMath::Sign(1.,muonTrack1->GetInverseBendingMomentum()));      
    etaM1 = muonTrack1->Eta();
    phiM1 = muonTrack1->Phi();
    thetaM1 = muonTrack1->Theta();
    fMuonNtuple->Fill(thetaX1,thetaY1,ptInv1,etaM1,phiM1,thetaM1,pxM1,pyM1,pzM1);
    // second loop over muon tracks in event
    for (Int_t mTrack2 = mTrack+1; mTrack2 < nmt; mTrack2++) {
      thetaX2=0; thetaY2=0; ptInv2=0; pYZ2=0; chargeM2=0;
      pxM2=0; pyM2=0; pzM2=0; energyM2=0;
      AliESDMuonTrack* muonTrack2 = fESD->GetMuonTrack(mTrack2);
      thetaX2 = muonTrack2->GetThetaX();
      thetaY2 = muonTrack2->GetThetaY();
      ptInv2 = TMath::Abs(muonTrack2->GetInverseBendingMomentum());
      pYZ2 = 1/ptInv2;
      pzM2  = - pYZ2 / TMath::Sqrt(1.0 + TMath::Tan(thetaY2)*TMath::Tan(thetaY2));
      pxM2  = pzM2 * TMath::Tan(thetaX2);
      pyM2  = pzM2 * TMath::Tan(thetaY2);
      chargeM2 = Int_t(TMath::Sign(1.,muonTrack2->GetInverseBendingMomentum()));      
      energyM2 = TMath::Sqrt(muonMass*muonMass + pxM2*pxM2 + pyM2*pyM2 + pzM2*pzM2);
      // if muons have opposite charge
      if(chargeM1*chargeM2 == -1){
	lvM2.SetPxPyPzE(pxM2, pyM2, pzM2,energyM2);
	lvJpsi = lvM1 + lvM2;
	invM = lvJpsi.M();
	fInvMass->Fill(invM);
      } // end if muons with opposite charge
    } // end second loop over muon tracks in event
  } // end first loop over muon tracks in event


  // Created AOD
  Float_t pxAodM1=0;   Float_t pyAodM1=0;     Float_t pzAodM1=0;
  Float_t etaAodM1=0;  Float_t phiAodM1=0;    Float_t thetaAodM1=0;
  Int_t chargeAodM1=0;  Float_t energyAodM1=0; TLorentzVector lvAodM1;
  Float_t pxAodM2=0;   Float_t pyAodM2=0;     Float_t pzAodM2=0;
  Float_t etaAodM2=0;  Float_t phiAodM2=0;    Float_t thetaAodM2=0;
  Int_t chargeAodM2=0;  Float_t energyAodM2=0; TLorentzVector lvAodM2;
  //for J/psi
  Float_t invMAOD=0;
  TLorentzVector lvJpsiAOD;

  int nmtAOD = fAOD->GetNumberOfTracks();
  // first loop over tracks
  for (Int_t mTrack = 0; mTrack < nmtAOD; mTrack++) {
    pxM1=0; pyM1=0; pzM1=0; 
    AliAODTrack* muonTrack1 = fAOD->GetTrack(mTrack);
    if(muonTrack1->IsMuonTrack()){
      etaAodM1 = muonTrack1->Eta();
      phiAodM1 = muonTrack1->Phi();
      thetaAodM1 = muonTrack1->Theta();
      pxAodM1 = muonTrack1->Px();
      pyAodM1 = muonTrack1->Py();
      pzAodM1 = muonTrack1->Pz();
      chargeAodM1 = muonTrack1->Charge();
      energyAodM1 = TMath::Sqrt(muonMass*muonMass + pxAodM1*pxAodM1 + pyAodM1*pyAodM1 + pzAodM1*pzAodM1);
      lvAodM1.SetPxPyPzE(pxAodM1, pyAodM1, pzAodM1,energyAodM1);
      fMuonNtupleAOD->Fill(etaAodM1,phiAodM1,thetaAodM1,pxAodM1,pyAodM1,pzAodM1);
    }
    // second loop over tracks in event
    for (Int_t mTrack2 = mTrack+1; mTrack2 < nmtAOD; mTrack2++) {
      chargeAodM2=0;
      pxAodM2=0; pyAodM2=0; pzAodM2=0; energyAodM2=0;
      AliAODTrack* muonTrack2 = fAOD->GetTrack(mTrack2);
      if(muonTrack2->IsMuonTrack()){
	etaAodM2 = muonTrack2->Eta();
	phiAodM2 = muonTrack2->Phi();
	thetaAodM2 = muonTrack2->Theta();
	pxAodM2 = muonTrack2->Px();
	pyAodM2 = muonTrack2->Py();
	pzAodM2 = muonTrack2->Pz();
 	energyAodM2 = TMath::Sqrt(muonMass*muonMass + pxAodM2*pxAodM2 + pyAodM2*pyAodM2 + pzAodM2*pzAodM2);
	chargeAodM2 = muonTrack2->Charge();
	// if muons of opposite charge
	if(chargeAodM1*chargeAodM2 == -1){
	  lvAodM2.SetPxPyPzE(pxAodM2, pyAodM2, pzAodM2,energyAodM2);
	  lvJpsiAOD = lvAodM1 + lvAodM2;
	  invMAOD = lvJpsiAOD.M();
	  fInvMassAOD->Fill(invMAOD);
	}// end if muons with opposite charge
      }// end if muon track
    }// end second loop over tracks in event  
  }// end first loop over tracks in event
  
  // Post final data. Write histo list to a file with option "RECREATE"
  PostData(0,fList);
}      

//________________________________________________________________________
void AliAnalysisTaskAODvsESD::Terminate(const Option_t*)
{
  // check if major differences between the two Ntuples (more comparisons can be added)
  int n1 = fMuonNtuple->GetEntries();
  int n2 = fMuonNtupleAOD->GetEntries();
 
  if(n1!=n2){
    printf("ERROR: Different number of entries in single muon Ntuples\n");
    return;
  }
  else
    printf("Same number of entries in single muon Ntuples\n");
  
//   TCanvas* cv1 = new TCanvas("cvn1","cvn1",500,350);
//   cv1->cd(1);
//   fInvMass->SetMarkerStyle(29);
//   fInvMass->SetMarkerColor(4);
//   fInvMass->Draw();
//   fInvMassAOD->SetMarkerStyle(28);
//   fInvMassAOD->SetMarkerColor(2);
//   fInvMassAOD->Draw("same"); 
}

