/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//
// class to benchmark the V0 pid capabilties
// runs over reconstructed V0 candidates and uses MC information for 
// further analysis [purity, PID perfortmance]
//
// authors:
//    Matus Kalisky <matus.kalisky@cern.ch>
//

#include "TIterator.h"

#include "AliVParticle.h"
#include "AliESDtrack.h"
#include "AliMCParticle.h"
#include "AliMCEvent.h"

#include "AliHFEcollection.h"

#include "AliHFEV0pidMC.h"
ClassImp(AliHFEV0pidMC)

//____________________________________________________________
  AliHFEV0pidMC::AliHFEV0pidMC():
    fMC(0x0)
    , fColl(NULL)
    , fV0cuts(NULL)
{
  //
  // default constructor
  //

}
//____________________________________________________________
AliHFEV0pidMC::~AliHFEV0pidMC(){
  //
  // destructor
  //
  if (fColl) delete fColl;
  //if (fV0cuts) delete fV0cuts;
}
//____________________________________________________________
void AliHFEV0pidMC::Init(){
  //
  // initialize objects
  //

  fColl = new AliHFEcollection("V0pidMC", "MC based V0 benchmarking");
  // QA
  fColl->CreateTH1F("h_QA_nParticles", "QA on track processing", 10, -0.5, 9.5);

  const Int_t nBins = 20;
  const Float_t pMin = 0.1;
  const Float_t pMax = 10.;

  // before PID
  fColl->CreateTH1F("h_Electron", "all electron candidates (no MC); p (GeV/c); counts", nBins, pMin, pMax, 0);
  fColl->CreateTH1F("h_PionK0", "all K0 pion candidates (no MC); p (GeV/c); counts",  nBins, pMin, pMax, 0);
  fColl->CreateTH1F("h_PionL", "all Lambda pion candidates (no MC); p (GeV/c); counts", nBins, pMin, pMax, 0);
  fColl->CreateTH1F("h_Kaon", "all Kaon candidates (no MC); p (GeV/c); counts",  nBins, pMin, pMax, 0);
  fColl->CreateTH1F("h_Proton", "all Lambda proton candidates (no MC); p (GeV/c); counts",  nBins, pMin, pMax, 0);
  
  fColl->CreateTH1F("h_mis_Electron", "all NON electron candidates MC tagged; p (GeV/c); counts",  nBins, pMin, pMax, 0);
  fColl->CreateTH1F("h_mis_PionK0", "all NON K0 pion candidates MC tagged; p (GeV/c); counts",  nBins, pMin, pMax, 0);
  fColl->CreateTH1F("h_mis_PionL", "all NON Lambda pion candidates MC tagged ; p (GeV/c); counts", nBins, pMin, pMax, 0);
  fColl->CreateTH1F("h_mis_Kaon", "all NON Kaon candidates MC tagged; p (GeV/c); counts", nBins, pMin, pMax, 0);
  fColl->CreateTH1F("h_mis_Proton", "all NON Lambda proton candidates MC tagged; p (GeV/c); counts", nBins, pMin, pMax, 0);  

  fColl->CreateTH1Fvector1(5, "h_tag_Electron", "electron candidate MC tagged; p (GeV/c); counts", nBins, pMin, pMax, 0);
  fColl->CreateTH1Fvector1(5, "h_tag_PionK0", "K0 pion candidate MC tagged; p (GeV/c); counts", nBins, pMin, pMax, 0);
  fColl->CreateTH1Fvector1(5, "h_tag_PionL", "Lambda pion candidate MC tagged; p (GeV/c); counts", nBins, pMin, pMax, 0);
  fColl->CreateTH1Fvector1(5, "h_tag_Kaon", "kaon candidate MC tagged; p (GeV/c); counts", nBins, pMin, pMax, 0);
  fColl->CreateTH1Fvector1(5, "h_tag_Proton", "Lambda proton candidate MC tagged; p (GeV/c); counts", nBins, pMin, pMax, 0);

}
//____________________________________________________________
Bool_t  AliHFEV0pidMC::Process(TObjArray * const particles, Int_t type){
  //
  // process the selected V0 daughter tracks
  //
  
  TString hname;
  const TString typeName[5] = {"Electron", "PionK0", "PionL", "Kaon", "Proton"};
  const Int_t  typePID[5] = {0, 2, 2, 3, 4};

  if(!fMC) return kFALSE;
  if(!particles) return kFALSE;
  
  AliVParticle *recTrack = NULL;
  TIterator *trackIter = particles->MakeIterator(); 
  while((recTrack = dynamic_cast<AliVParticle *>(trackIter->Next()))){
    fColl->Fill("h_QA_nParticles", 0);
    // only ESD for now
    AliESDtrack *track = dynamic_cast<AliESDtrack *>(recTrack);
    if(!track) continue;
    const AliExternalTrackParam *ext = track->GetOuterParam();
    if(!ext) continue;
    // MC label
    Int_t label = track->GetLabel();
    if(label <0){
       fColl->Fill("h_QA_nParticles", 1);
      continue;
    }
    AliMCParticle *mcpD = dynamic_cast<AliMCParticle*>(fMC->GetTrack(label));
    if(!mcpD){
      fColl->Fill("h_QA_nParticles", 2);
      continue;
    }

    Float_t p = ext->P();
    //Short_t charge = ext->Charge();
    Int_t pdgD = mcpD->PdgCode();
    AliMCParticle *mcpM = dynamic_cast<AliMCParticle*>(fMC->GetTrack(mcpD->GetMother()));
    if(!mcpM){
      fColl->Fill("h_QA_nParticles", 3);
      continue;
    }
    //Int_t pdgM = mcpM->PdgCode();

    // all candidates
    hname = "h_" + typeName[type];
    fColl->Fill(hname, p);
    Int_t pidD = PDGtoPIDdaughter(pdgD);
    
    // all misidentified candidates
    hname = "h_mis_" + typeName[type];
    if(typePID[type] != pidD){
      fColl->Fill(hname, p);
    }

    // for every particle fill detailed information about the missidentified particles
    hname = "h_tag_" + typeName[type];
    Int_t aliPID = PDGtoAliPID(pdgD);
    fColl->Fill(hname, aliPID, p);
    
  }// .. loop over array
  
  
  return kTRUE;
}  
//____________________________________________________________
Int_t AliHFEV0pidMC::PDGtoPIDdaughter(Int_t pdg) const {
  //
  // convert PDG to local pid 
  //
  switch (TMath::Abs(pdg)){
  case 11:
    return 0;  // electron gamma
  case 211:
    return 2; // pion K0 or pion Lambda
  case 321:
    return 3; //kaon Phi
  case 2212:
    return 4; // proton Lambda
  default:
    return -1;
  };
  
  return -1;
}
//____________________________________________________________
Int_t AliHFEV0pidMC::PDGtoPIDmother(Int_t pdg) const {
  //
  // convert PDG to local pid
  //
  switch (TMath::Abs(pdg)){
  case 22:
      return 0; // gamma
    case 310: 
    return 1; // K0s
  case 333:
    return 2; // Phi
  case 3122:
    return 3; // Lambda
  default:
    return -1;
  };
  
  return -1;
}
//____________________________________________________________
Int_t AliHFEV0pidMC::PDGtoAliPID(Int_t pdg) const {
  //
  // convert PDG to AliPID
  //
  
  switch (TMath::Abs(pdg)){
  case 11:
    return 0; // electron 
  case 13:
    return 1; // muon
  case 211:
    return 2; // pion 
  case 321:
    return 3; // kaon 
  case 2212:
    return 4; // proton 
  default:
    return -1;
  };
  return -1;
}

