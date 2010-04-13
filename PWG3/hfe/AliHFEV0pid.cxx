/*************************************************************************
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
//
// Utility class for V0 PID
// Identifies Electrons, Pions and Protons using gamma conversions and
// the decays of K0s and Lambdas
// Containers with samples of Electrons, Pions and Protons can be accessed
// via GetListOfElectrons() etc.
//
// Authors:
//    Matus Kalisky <matus.kalisky@cern.ch>
//    Markus Heide <mheide@uni-muenster.de>
//    Markus Fasel <M.Fasel@gsi.de>
//
#include <TDatabasePDG.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TString.h>

#include <TDatabasePDG.h>

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODVertex.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliESDVertex.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliVEvent.h"
#include "AliVTrack.h"

#include "AliHFEV0pid.h"
ClassImp(AliHFEV0pid)

//____________________________________________________________
AliHFEV0pid::AliHFEV0pid():
  TObject()
  , fInputEvent(NULL)
  , fPrimaryVertex(NULL)
  , fElectrons(NULL)
  , fPionsK0(NULL)
  , fPionsL(NULL)
  , fKaons(NULL)
  , fProtons(NULL)
  , fIndices(NULL)
  , fQA(NULL)
{
  //
  // Default constructor
  //
  fElectrons = new TObjArray();
  fPionsK0 = new TObjArray();
  fPionsL = new TObjArray();
  fKaons = new TObjArray();
  fProtons = new TObjArray();
  fIndices = new AliHFEV0pidTrackIndex;
}

//____________________________________________________________
AliHFEV0pid::~AliHFEV0pid(){
  //
  // Destructor
  // Remove Containers
  //
  if(fInputEvent) delete fInputEvent;
  //if(fPrimaryVertex) delete fPrimaryVertex;
  if(fElectrons) delete fElectrons;
  if(fPionsK0) delete fPionsK0;
  if(fPionsL) delete fPionsL;
  if(fKaons) delete fKaons;
  if(fProtons) delete fProtons;
  if(fIndices) delete fIndices;
  if(fQA) delete fQA;
}

//____________________________________________________________
void AliHFEV0pid::InitQA(){
  //
  // Initialize QA histograms
  //
  if(!fQA){
    fQA = new AliHFEcollection("v0pidQA", "QA histograms for V0 PID");

    // QA histograms for cut statistics
    fQA->CreateTH1F("h_cutEfficiencyGamma", "Cut Efficiency for Gammas", 10, 0, 10);
    fQA->CreateTH1F("h_cutEfficiencyK0s", "Cut Efficiency for K0s", 10, 0, 10);
    fQA->CreateTH1F("h_cutEfficiencyPhi", "Cut Efficiency for Phi", 10, 0, 10);
    fQA->CreateTH1F("h_cutEfficiencyLambda", "Cut Efficiency for Lambdas", 10, 0, 10);

    // QA histograms for invariant mass
    fQA->CreateTH1F("h_InvMassGamma", "Gamma invariant mass; inv mass [GeV/c^{2}]; counts", 100, 0, 0.25);
    fQA->CreateTH1F("h_InvMassK0s", "K0s invariant mass; inv mass [GeV/c^{2}]; counts", 100, 0.4, 0.65);
    fQA->CreateTH1F("h_InvMassPhi", "Phi invariant mass; inv mass [GeV/c^{2}]; counts", 100, 0.4, 0.65);
    fQA->CreateTH1F("h_InvMassLambda", "Lambda invariant mass; inv mass [GeV/c^{2}]; counts", 100, 1.05, 1.15);
    
    // QA histograms for p distribution (of the daughters)
    fQA->CreateTH1F("h_P_electron", "P distribution of the gamma electrons; p (GeV/c); counts", 100, 0.1, 10);
    fQA->CreateTH1F("h_P_K0pion", "P distribution of the K0 pions; p (GeV/c); counts", 100, 0.1, 10);
    fQA->CreateTH1F("h_P_Lpion", "P distribution of the Lambda pions; p (GeV/c); counts", 100, 0.1, 10);
    fQA->CreateTH1F("h_P_Lproton", "P distribution of the Lambda protons; p (GeV/c); counts", 100, 0.1, 10);

    // QA invariant mass as a functin of pt
    fQA->CreateTH1Fvector1(20, "h_InvMassGamma_pt", "Gamma invarinat mass in pt bins; inv mass [GeV/c^{2}]; counts", 250, 0, 2);
    fQA->CreateTH1Fvector1(20, "h_InvMassK0_pt", "K0 invarinat mass in pt bins; inv mass [GeV/c^{2}]; counts", 250, 0, 2);
    fQA->CreateTH1Fvector1(20, "h_InvMassPhi_pt", "Phi invarinat mass in pt bins; inv mass [GeV/c^{2}]; counts", 250, 0, 2);
    fQA->CreateTH1Fvector1(20, "h_InvMassLambda_pt", "Lambda invarinat mass in pt bins; inv mass [GeV/c^{2}]; counts", 250, 0, 2);

    // QA pt of the V0
    fQA->CreateTH1F("h_Pt_Gamma", "Pt of the gamma conversion; p_{T} (GeV/c); counts", 100, 0, 10);
    fQA->CreateTH1F("h_Pt_K0", "Pt of the K0; p_{T} (GeV/c); counts", 100, 0, 10);
    fQA->CreateTH1F("h_Pt_Phi", "Pt of the Phi; p_{T} (GeV/c); counts", 100, 0, 10);
    fQA->CreateTH1F("h_Pt_Lambda", "Pt of the Lambda; p_{T} (GeV/c); counts", 100, 0, 10);
    //fQA->CreateTH1F("h_Pt_electrons", "Pt of the conversion electrons; p_{T} (GeV/c); counts");
    //fQA->CreateTH1F("h_Pt_pionsK0", "Pt of the K0 pions; p_{T} (GeV/c); counts");
    //fQA->CreateTH1F("h_Pt_pionsL", "Pt of the Lambda pions; p_{T} (GeV/c); counts");
    //fQA->CreateTH1F("h_Pt_protons", "Pt of the Lambda protons; p_{T} (GeV/c); counts");


    // QA histogram for both Lambda candidate combinations - 
    fQA->CreateTH2F("h_L0_dca_v_dMass", "L0 dca verus dMass; dMass [GeV/c^{2}]; dDCA [cm]; ", 100, -1., 1., 100, 0., 5.);
    
    // Chi2 histograms
    fQA->CreateTH1F("h_chi2_gamma", "Chi2 for gammas", 10000, 0, 1000);
    fQA->CreateTH1F("h_chi2_K0s", "Chi2 for K0s", 10000, 0, 500);
    fQA->CreateTH1F("h_chi2_Phi", "Chi2 for K0s", 10000, 0, 500);
    fQA->CreateTH1F("h_chi2_Lambda", "Chi2 for Lambdas", 10000, 0, 1000);
  }
}

//____________________________________________________________
void AliHFEV0pid::Process(AliVEvent * const inputEvent){

  //
  // Find protons, pions and electrons using V0 decays and 
  // store the pointers in the TObjArray
  //
  Int_t nGamma = 0, nK0s = 0, nLambda = 0, nPhi = 0;
  fInputEvent = inputEvent;
  
  fIndices->Init(fInputEvent->GetNumberOfV0s() * 2);
  fPrimaryVertex = new AliKFVertex(*(fInputEvent->GetPrimaryVertex()));
  Int_t v0status = 0;
  for(Int_t iv0 = 0; iv0 < fInputEvent->GetNumberOfV0s(); iv0++){
    if(!TString(fInputEvent->IsA()->GetName()).CompareTo("AliESDEvent")){
      // case ESD
      SetESDanalysis();
      AliESDv0 *esdV0 = (dynamic_cast<AliESDEvent *>(fInputEvent))->GetV0(iv0);
      if(!esdV0->GetOnFlyStatus()) continue; // Take only V0s from the On-the-fly v0 finder
      v0status = ProcessV0(esdV0);
    } else {
      // case AOD
      SetAODanalysis();
      AliAODv0 *aodV0 = (dynamic_cast<AliAODEvent *>(fInputEvent))->GetV0(iv0);
      if(aodV0->GetOnFlyStatus()) continue; // Take only V0s from the On-the-fly v0 finder
      v0status = ProcessV0(aodV0);
    }
    switch(v0status){
    case kRecoGamma: nGamma++; break;
    case kRecoK0s: nK0s++; break;
    case kRecoPhi: nPhi++; break;  
    case kRecoLambda: nLambda++; break;
    };
  }

  AliDebug(1, Form("Number of gammas  : %d", nGamma));
  AliDebug(1, Form("Number of K0s     : %d", nK0s));
  AliDebug(1, Form("Number of Phis    : %d", nPhi));
  AliDebug(1, Form("Number of Lambdas : %d", nLambda));

  AliDebug(1, "Number of stored tracks:");
  AliDebug(1, Form("Number of electrons      : %d", fElectrons->GetEntries()));
  AliDebug(1, Form("Number of K0 pions       : %d", fPionsK0->GetEntries()));
  AliDebug(1, Form("Number of Lambda pions   : %d", fPionsL->GetEntries()));
  AliDebug(1, Form("Number of Phi kaons      : %d", fKaons->GetEntries()));
  AliDebug(1, Form("Number of protons        : %d", fProtons->GetEntries()));
  
  delete  fPrimaryVertex;
}

//____________________________________________________________
Int_t AliHFEV0pid::ProcessV0(TObject *v0){
  //
  // Process single V0
  // Apply general cut and special cuts for gamma, K0s, Lambda
  //
  AliVTrack* daughter[2];
  daughter[0] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack((dynamic_cast<AliESDv0 *>(v0))->GetPindex()));
  daughter[1] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack((dynamic_cast<AliESDv0 *>(v0))->GetPindex()));
  if(!daughter[0] || !daughter[1]) return kUndef;

  if(IsESDanalysis()){
    for(Int_t i=0; i<2; ++i){
      if(!CutESDtrack(dynamic_cast<AliESDtrack*>(daughter[i]))) return kUndef;
    }    
  }

  if(IsGammaConv(v0)) return kRecoGamma;
  else if(IsK0s(v0)) return kRecoK0s;
  else if(IsLambda(v0)) return kRecoLambda;
  else return kUndef;
  
  
}
//____________________________________________________________
void AliHFEV0pid::Flush(){
  //
  // Clear the Lists
  //
  AliDebug(1, "Flushing containers");
  fProtons->Clear();
  fPionsK0->Clear();
  fPionsL->Clear();
  fElectrons->Clear();
  fIndices->Flush();
}

//____________________________________________________________
Bool_t AliHFEV0pid::IsGammaConv(TObject *v0){
  //
  // Identify Gamma
  //
  AliVTrack* daughter[2];
  Int_t pIndex = 0, nIndex = 0;
  Double_t invMass = 0.;
  if(IsESDanalysis()){
    // ESD - cut V0
    AliESDv0 *esdV0 = dynamic_cast<AliESDv0 *>(v0);
    if(!CutV0(esdV0, kRecoGamma)) return kFALSE; 
    if(LooseRejectK0(esdV0) || LooseRejectLambda(esdV0)) return kFALSE;
    // DEBUG
    //invMass = esdV0->GetEffMass(AliPID::kElectron, AliPID::kElectron);
    invMass = GetEffMass(esdV0, AliPID::kElectron, AliPID::kElectron);
    //..

    pIndex = esdV0->GetPindex();
    nIndex = esdV0->GetNindex();
  } else {
    // AOD Analysis - not possible to cut
    AliAODv0 *aodV0 = dynamic_cast<AliAODv0 *>(v0);
    pIndex = aodV0->GetPosID();
    nIndex = aodV0->GetNegID();
    invMass = aodV0->InvMass2Prongs(0, 1, kElectron, kElectron);
  }
  daughter[0] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(pIndex));
  daughter[1] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(nIndex));
  if(!daughter[0] || !daughter[1]) return kFALSE;

  // Get Invariant mass and chi2/ndf 
  AliKFParticle *kfMother = CreateMotherParticle(daughter[0], daughter[1], TMath::Abs(kElectron), TMath::Abs(kElectron));
  kfMother->SetProductionVertex(*fPrimaryVertex);
  kfMother->SetMassConstraint(0, 1.);
  Double_t ndf = kfMother->GetNDF();
  Double_t chi2 = kfMother->GetChi2();
  delete kfMother; 

  if(fQA) fQA->Fill("h_chi2_gamma", chi2/ndf);
  if(chi2/ndf > 7) return kFALSE;
  Double_t mPt = kfMother->GetPt();
  fQA->Fill("h_Pt_Gamma", mPt);
  Int_t ptBin = (int)(mPt*10.0);

  if(fQA) fQA->Fill("h_InvMassGamma", invMass);
  if(invMass > 0.05) return kFALSE; 
  fQA->Fill("h_InvMassGamma_pt", ptBin+1, invMass);
  
  AliDebug(1, Form("Gamma identified, daughter IDs: %d,%d", daughter[0]->GetID(), daughter[1]->GetID()));
  // Identified gamma - store tracks in the electron containers
  if(!fIndices->Find(daughter[0]->GetID())){
    AliDebug(1, Form("Gamma identified, daughter IDs: %d,%d", daughter[0]->GetID(), daughter[1]->GetID()));    
    fElectrons->Add(daughter[0]);
    fIndices->Add(daughter[0]->GetID(), AliHFEV0pid::kRecoElectron);
  }
  if(!fIndices->Find(daughter[1]->GetID())){
    AliDebug(1, Form("Gamma identified, daughter IDs: %d,%d", daughter[1]->GetID(), daughter[1]->GetID()));
    fElectrons->Add(daughter[1]);
    fIndices->Add(daughter[1]->GetID(), AliHFEV0pid::kRecoElectron);
  }
  return kTRUE;
}



//____________________________________________________________
Bool_t AliHFEV0pid::IsK0s(TObject *v0){
  //
  // Identify K0s
  //
  const Double_t kK0smass=TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();  // PDG K0s mass
  AliVTrack* daughter[2];
  Int_t pIndex = 0, nIndex = 0;
  Double_t invMass = 0.;
  if(IsESDanalysis()){
    // ESD - cut V0
    AliESDv0 *esdV0 = dynamic_cast<AliESDv0 *>(v0);
    if(!CutV0(esdV0, kRecoK0s)) return kFALSE; 
    invMass = esdV0->GetEffMass(AliPID::kPion, AliPID::kPion);
    pIndex = esdV0->GetPindex();
    nIndex = esdV0->GetNindex();
  } else {
    // AOD Analysis - not possible to cut
    AliAODv0 *aodV0 = dynamic_cast<AliAODv0 *>(v0);
    pIndex = aodV0->GetPosID();
    nIndex = aodV0->GetNegID();
    invMass = aodV0->MassK0Short();
  }
  daughter[0] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(pIndex));
  daughter[1] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(nIndex));
  if(!daughter[0] || !daughter[1]) return kFALSE;

  // Get Invariant mass and chi2/ndf 
  AliKFParticle *kfMother = CreateMotherParticle(daughter[0], daughter[1], TMath::Abs(kPiPlus), TMath::Abs(kPiPlus));
  kfMother->SetProductionVertex(*fPrimaryVertex);
  kfMother->SetMassConstraint(kK0smass, 0.);
  Double_t ndf = kfMother->GetNDF();
  Double_t chi2 = kfMother->GetChi2();
  delete kfMother; 

  if(fQA) fQA->Fill("h_chi2_K0s", chi2/ndf);
  if(chi2/ndf > 5) return kFALSE;
  Double_t mPt = kfMother->GetPt();
  fQA->Fill("h_Pt_K0", mPt);
  Int_t ptBin = (int)(mPt*10.0);
  fQA->Fill("h_InvMassK0_pt", ptBin+1, invMass);

  if(fQA) fQA->Fill("h_InvMassK0s", invMass);

  if(invMass < 0.485 || invMass > 0.51) return kFALSE; 
  AliDebug(1, Form("K0 identified, daughter IDs: %d,%d", daughter[0]->GetID(), daughter[1]->GetID()));

  // Identified gamma - store tracks in the electron containers
  if(!fIndices->Find(daughter[0]->GetID())){
    AliDebug(1, Form("Adding K0 Pion track with ID %d", daughter[0]->GetID()));
    fPionsK0->Add(daughter[0]);
    fIndices->Add(daughter[0]->GetID(), AliHFEV0pid::kRecoPionK0);
  }
  if(!fIndices->Find(daughter[1]->GetID())){
    AliDebug(1, Form("Adding K0 Pion track with ID %d", daughter[1]->GetID()));
    fPionsK0->Add(daughter[1]);
    fIndices->Add(daughter[1]->GetID(), AliHFEV0pid::kRecoPionK0);
  }
  return kTRUE; 
}

//____________________________________________________________
Bool_t AliHFEV0pid::IsPhi(TObject *v0){
  //
  // Identify Phi - very preliminary - requires diffrent approach (V0 fnder is not effective)
  //

  //const Double_t kPhiMass=TDatabasePDG::Instance()->GetParticle(333)->Mass();  // PDG phi mass
  //AliVTrack* daughter[2];
  //AliKFParticle *mother = NULL;
  //Double_t invMass = 0.;
 
  Int_t pIndex = 0, nIndex = 0;
  if(IsESDanalysis()){
    // ESD - cut V0
    AliESDv0 *esdV0 = dynamic_cast<AliESDv0 *>(v0);
    if(!CutV0(esdV0, kRecoPhi)) return kFALSE; 
    if(LooseRejectGamma(esdV0) || LooseRejectK0(esdV0)) return kFALSE;
    pIndex = esdV0->GetPindex();
    nIndex = esdV0->GetNindex();
  } else {
    // PRELIMINARY - !!!
    // AOD Analysis - not possible to cut
  } 

  return kTRUE;
}

//____________________________________________________________
Bool_t AliHFEV0pid::IsLambda(TObject *v0){
  //
  // Identify Lambda
  //
  const Double_t kL0mass=TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();  // PDG lambda mass
  AliVTrack* daughter[2];
  Int_t pIndex = 0, nIndex = 0;
  AliKFParticle *mother[2] = {NULL, NULL};
  Float_t dMass[2]; // lambda mass difference for the two hypoteses
  //Int_t   lambda = 0; // [1] for lambda and [-1] for anti-lambda
  Double_t chi2 = 0.;

  Double_t invMass = 0.;
  Double_t invMassEOD = 0;
  if(IsESDanalysis()){
    // ESD - cut V0
    AliESDv0 *esdV0 = dynamic_cast<AliESDv0 *>(v0);
    if(!CutV0(esdV0, kRecoLambda)) return kFALSE; 
    if(LooseRejectK0(esdV0) || LooseRejectGamma(esdV0)) return kFALSE;
    pIndex = esdV0->GetPindex();
    nIndex = esdV0->GetNindex();
    
    //invMass = esdV0->GetEffMass(AliPID::kPion, AliPID::kPion);

  } else {
    // PRELIMINARY - !!!
    // AOD Analysis - not possible to cut
    
    // again - two cases as above
    AliAODv0 *aodV0 = dynamic_cast<AliAODv0 *>(v0);
    pIndex = aodV0->GetPosID();
    nIndex = aodV0->GetNegID();
    invMassEOD = aodV0->MassLambda();
  }

  daughter[0] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(pIndex));
  daughter[1] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(nIndex));
  if(!daughter[0] || !daughter[1]) return kFALSE;

  //
  // now - go over two case - Lambda and AntiLambda
  // choose the on ewhere the resulting lambda mass is close to the
  // expected value and at the same time points closer to the primary vertex in XY plane
  //
  // A)      lambda -> proton + negative-pion
  // B) anti-lambda -> anti-proton + positive-pion
  //    

  //
  // A) :: proton
  //
  mother[0] = CreateMotherParticle(daughter[0], daughter[1], TMath::Abs(kProton), TMath::Abs(kPiPlus));
  AliHFELambdaInf lambda0(mother[0], fPrimaryVertex);

  // Undo changes
  *fPrimaryVertex -= *mother[0];
  AddTrackToKFVertex(daughter[0], TMath::Abs(kProton));
  AddTrackToKFVertex(daughter[1], TMath::Abs(kPiPlus));

  //
  // B) :: anti-proton
  //
  mother[1] = CreateMotherParticle(daughter[0], daughter[1], TMath::Abs(kPiPlus), TMath::Abs(kProton));
  AliHFELambdaInf lambda1(mother[1], fPrimaryVertex);

  // qa histograms
  dMass[0] = lambda0.GetInvariantMass() - kL0mass;
  dMass[1] = lambda1.GetInvariantMass() - kL0mass;
  fQA->Fill("h_L0_dca_v_dMass", dMass[0], lambda0.GetDistanceFromPrimaryVertex());
  fQA->Fill("h_L0_dca_v_dMass", dMass[1], lambda1.GetDistanceFromPrimaryVertex());

  //
  // decide on the true Lambda candidate bsaed on
  // - mass difference
  // - vertex dca  (testing needed first)
  //
  AliHFELambdaInf *lambdaInf[2] = {&lambda0, &lambda1};
  Bool_t isLambda = kFALSE;
  Int_t index = -1;
  index = (dMass[0] <dMass[1]) ? 0 : 1;
  invMass = lambdaInf[index]->GetInvariantMass();
  chi2 = lambdaInf[index]->GetChi2NDF();
  //if(!IsESDanalysis()){
  //  invMass = invMassEOD;
  //}
  //else{
  //  invMass = mother[index]->GetMass();
  //}
  if(fQA) fQA->Fill("h_chi2_Lambda", chi2);
  if(chi2 < 3){
    if(fQA) fQA->Fill("h_InvMassLambda", invMass);
    Double_t mPt = mother[index]->GetPt();
    fQA->Fill("h_Pt_Lambda", mPt);
    Int_t ptBin = (int)(mPt*10.0);
    fQA->Fill("h_InvMassLambda_pt", ptBin+1, invMass);
    
    // cut on the invariant mass for the proton and pion selection
    if(invMass > 1.11 || invMass < 1.12){
      // Identified lambdas - store the protons and pions and update primary vertex
      *fPrimaryVertex += *mother[index];
    
      // lambda
      if(0 == index){
        if(!fIndices->Find(daughter[0]->GetID())){
	        fProtons->Add(daughter[0]);
	        fIndices->Add(daughter[0]->GetID(), AliHFEV0pid::kRecoProton);
        }
        if(!fIndices->Find(daughter[1]->GetID())){
	        fPionsL->Add(daughter[1]);
	        fIndices->Add(daughter[1]->GetID(), AliHFEV0pid::kRecoPionL);
        }
      }
      // antilambda
      if(1 == index){
        if(!fIndices->Find(daughter [1]->GetID())){
	        fProtons->Add(daughter[1]);
	        fIndices->Add(daughter[1]->GetID(), AliHFEV0pid::kRecoProton);
        }
        if(!fIndices->Find(daughter [0]->GetID())){
	        fPionsL->Add(daughter[0]);
	        fIndices->Add(daughter [0]->GetID(), AliHFEV0pid::kRecoPionL);
        }
      }
      isLambda = kTRUE;
    }
  }

  if(!isLambda){
    // Add the daughters again to the primary vertex
    AddTrackToKFVertex(daughter[0], TMath::Abs(kPiPlus));
    AddTrackToKFVertex(daughter[1], TMath::Abs(kProton));
  }
  // remove the objects
  for(Int_t i=0; i<2; ++i){
    if (mother[i]) delete mother[i];
  }
 
  return isLambda;
}

//____________________________________________________________
AliKFParticle *AliHFEV0pid::CreateMotherParticle(AliVTrack *pdaughter, AliVTrack *ndaughter, Int_t pspec, Int_t nspec){
  //
  // Creates a mother particle
  //
  AliKFParticle pkfdaughter(*pdaughter, pspec);
  AliKFParticle nkfdaughter(*ndaughter, nspec);

  // check if the daughter particles are coming from the primary vertex
  if(IsESDanalysis()){
    // ESD Analyis
    const AliESDVertex *esdvertex = dynamic_cast<const AliESDVertex *>(fInputEvent->GetPrimaryVertex());
    UShort_t *contrib = esdvertex->GetIndices();

    Int_t nfound = 0;
    for(Int_t id = 0; id < esdvertex->GetNIndices(); id++){
      if(contrib[id] == pdaughter->GetID()){
        *fPrimaryVertex -= pkfdaughter;
        nfound++;
      }
      if(contrib[id] == ndaughter->GetID()){
        *fPrimaryVertex -= nkfdaughter;
        nfound++;
      }
      if(nfound == 2) break;
    }
  } else {
    // AOD Analysis: AOD Vertex 
    const AliAODVertex *aodvertex = dynamic_cast<const AliAODVertex *>(fInputEvent->GetPrimaryVertex());
    if(aodvertex->HasDaughter(pdaughter))
      *fPrimaryVertex -= pkfdaughter;
    if(aodvertex->HasDaughter(ndaughter))
      *fPrimaryVertex -= nkfdaughter;
  }

  // Create the mother particle and add them to the primary vertex
  AliKFParticle *mother = new AliKFParticle(pkfdaughter, nkfdaughter);
  *fPrimaryVertex += *mother;

  return mother;
}

//____________________________________________________________
void AliHFEV0pid::AddTrackToKFVertex(AliVTrack *track, Int_t species){
  //
  // Add track to the primary vertex (if it was used in the vertex
  // calculation)
  //
  Bool_t isAdd = kFALSE;
  if(IsESDanalysis()){
     const AliESDVertex *esdvertex = dynamic_cast<const AliESDVertex *>(fInputEvent->GetPrimaryVertex());
    UShort_t *contrib = esdvertex->GetIndices();
    for(Int_t id = 0; id < esdvertex->GetNIndices(); id++){
      if(contrib[id] == track->GetID()){
        isAdd = kTRUE;
        break;
      }
    }
  } else {
    const AliAODVertex *aodvertex = dynamic_cast<const AliAODVertex *>(fInputEvent->GetPrimaryVertex());
    if(aodvertex->HasDaughter(track)) isAdd = kTRUE;
  }
  if(isAdd){
    AliKFParticle kftrack(*track, species);
    *fPrimaryVertex += kftrack;
  }
}

//____________________________________________________________
Bool_t AliHFEV0pid::CutV0(AliESDv0 *v0, Int_t V0species){
  //
  // Cut the V0
  //
  Int_t cutRequired = 0;
  // For cut always take min and max
  Double_t  cutCosPoint[2] = {0., 0.}, 
  cutDCA[2] = {0., 0.}, 
  cutProdVtx[2] = {0., 0.},
	cutOpAng[2] = {0., 0.},
	cutPsiPair[2] = {0., 0.};
	  
  switch(V0species){
  case kRecoGamma:
    cutCosPoint[1] = 0.03;
    SETBIT(cutRequired, 1);
    cutDCA[1] = 0.25;
    SETBIT(cutRequired, 3);
    cutProdVtx[0] = 6;
    SETBIT(cutRequired, 4);
    cutOpAng[1] = 0.1;
    SETBIT(cutRequired, 7);
    cutPsiPair[1] = 0.05;
    SETBIT(cutRequired, 9);
    break;
  case kRecoK0s:
    cutCosPoint[1] = 0.03;
    SETBIT(cutRequired, 1);
    cutDCA[1] = 0.1;
    SETBIT(cutRequired, 3);
    cutProdVtx[1] = 8.1;
    SETBIT(cutRequired, 5);
    break;
  case kRecoPhi:
    break;
  case kRecoLambda:
    cutCosPoint[1] = 0.03;
    SETBIT(cutRequired, 1);
    cutDCA[1] = 0.2;
    SETBIT(cutRequired, 3);
    cutProdVtx[1] = 24;
    SETBIT(cutRequired, 5);
    break;
  default:
    // unidentified, return
    return kFALSE;
  };
  
  Char_t hname[256];
  const Char_t *specname[4] = {"Gamma", "K0s", ", Phi", "Lambda"};
  sprintf(hname, "h_cutEfficiency%s", specname[V0species-1]);

  // Cut on pointing angle
  Double_t cosPoint = v0->GetV0CosineOfPointingAngle();
  if(TESTBIT(cutRequired, 0) && TMath::ACos(cosPoint) < cutCosPoint[0]) return kFALSE; 
  if(fQA) fQA->Fill(hname, 0);
  if(TESTBIT(cutRequired, 1) && TMath::ACos(cosPoint) > cutCosPoint[1]) return kFALSE;
  if(fQA) fQA->Fill(hname, 1); 

  // Cut on DCA between daughters
  Double_t dca = v0->GetDcaV0Daughters();
  if(TESTBIT(cutRequired, 2) && dca < cutDCA[0]) return kFALSE;
  if(fQA) fQA->Fill(hname, 2);
  if(TESTBIT(cutRequired, 3) && dca > cutDCA[1]) return kFALSE;
  if(fQA) fQA->Fill(hname, 3);

  // Cut on reconstructed verted position
  Double_t x, y, z; 
  v0->GetXYZ(x,y,z);
  Double_t r = TMath::Sqrt(x*x + y*y);
  if(TESTBIT(cutRequired, 4) && r < cutProdVtx[0]) return kFALSE;
  if(fQA) fQA->Fill(hname, 4);
  if(TESTBIT(cutRequired, 5) && r > cutProdVtx[1]) return kFALSE;
  if(fQA) fQA->Fill(hname, 5);

  //Cut on Opening angle (conversions only)
  if(TESTBIT(cutRequired, 6) && OpenAngle(v0) < cutOpAng[0]) return kFALSE;
  if(fQA) fQA->Fill(hname, 6);
  if(TESTBIT(cutRequired, 7) && OpenAngle(v0) > cutOpAng[1]) return kFALSE;
  if(fQA) fQA->Fill(hname, 7);

  //Cut on PsiPair angle (conversons only)
  if(TESTBIT(cutRequired, 8) && PsiPair(v0) < cutPsiPair[0]) return kFALSE;
  if(fQA) fQA->Fill(hname, 8);
  if(TESTBIT(cutRequired, 9) && PsiPair(v0) > cutPsiPair[1]) return kFALSE;
  if(fQA) fQA->Fill(hname, 9);
  return kTRUE;
}

//____________________________________________________________
Bool_t AliHFEV0pid::CutESDtrack(AliESDtrack *track){
  // 
  // Single track cuts
  //
  // Hard coaded cut values for the beginning 
  //

  if(!track) return kFALSE;
  
  // status word
  ULong_t status = track->GetStatus();

  // DCA - to Vertex R & Z
  Float_t dcaR = -1.;
  Float_t dcaZ = -1.;
  track->GetImpactParameters(dcaR, dcaZ);
  //if(dcaR > 4.0) return kFALSE;
  //if(dcaZ > 10.0) return kFALSE;

  // No. of TPC clusters
  if(track->GetTPCNcls() < 80) return kFALSE;

  // TPC refit
  if(!(status & AliESDtrack::kTPCrefit)) return kFALSE;

  // Chi2 per TPC cluster
  Int_t nTPCclusters = track->GetTPCclusters(0);
  Float_t chi2perTPCcluster = track->GetTPCchi2()/Float_t(nTPCclusters);
  if(chi2perTPCcluster > 3.5) return kFALSE;

  // TPC cluster ratio
  Float_t cRatioTPC = track->GetTPCNclsF() > 0. ? static_cast<Float_t>(track->GetTPCNcls())/static_cast<Float_t> (track->GetTPCNclsF()) : 1.;
  if(cRatioTPC < 0.6) return kFALSE;

  // kinks
  if(track->GetKinkIndex(0) != 0) return kFALSE;

  // Covariance matrix - TO BE RECONSIDERED
  Double_t extCov[15];
  track->GetExternalCovariance(extCov);
  //if(extCov[0]  > 2. ) return kFALSE;
  //if(extCov[2]  > 2. ) return kFALSE;
  //if(extCov[5]  > 0.5) return kFALSE;
  //if(extCov[9]  > 0.5) return kFALSE;
  //if(extCov[14] > 2. ) return kFALSE;
  
  // pt
  if(track->Pt() < 0.1 || track->Pt() > 100) return kFALSE;

  // eta
  if(TMath::Abs(track->Eta()) > 0.9) return kFALSE;

  // the track made it through! :-)
  return kTRUE;
}

//_________________________________________________
Bool_t AliHFEV0pid::LooseRejectK0(AliESDv0 * const v0) const {
  //
  // Reject K0 based on loose cuts
  //
  Double_t mass = v0->GetEffMass(AliPID::kPion, AliPID::kPion);
  if(mass > 0.494 && mass < 0.501) return kTRUE;
  return kFALSE;
}

//_________________________________________________
Bool_t AliHFEV0pid::LooseRejectLambda(AliESDv0 * const v0) const {
  //
  // Reject Lambda based on loose cuts
  //
  Double_t mass1 = v0->GetEffMass(AliPID::kPion, AliPID::kProton);
  Double_t mass2 = v0->GetEffMass(AliPID::kProton, AliPID::kPion);
  
  if(mass1 > 1.1 && mass1 < 1.12) return kTRUE;
  if(mass2 > 1.1 && mass2 < 1.12) return kTRUE;
  return kFALSE;
}

//_________________________________________________
Bool_t AliHFEV0pid::LooseRejectGamma(AliESDv0 * const v0) const {
  //
  // Reject Lambda based on loose cuts
  //
 
  //Double_t mass = v0->GetEffMass(AliPID::kElectron, AliPID::kElectron);
  // DEBUG temporary solution, see the comment in GetEffMass
  Double_t mass = GetEffMass(v0, AliPID::kElectron, AliPID::kElectron);
  //..
  if(mass < 0.02) return kTRUE;
  return kFALSE;
}

//_________________________________________________
Float_t AliHFEV0pid::OpenAngle(AliESDv0 *v0) const {
  //
  // Opening angle between two daughter tracks
  //
  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0};
    

  v0->GetNPxPyPz(mn[0],mn[1],mn[2]);//reconstructed cartesian momentum components of negative daughter;
  v0->GetPPxPyPz(mp[0],mp[1],mp[2]);//reconstructed cartesian momentum components of positive daughter;

  
  Float_t openAngle = TMath::ACos((mp[0]*mn[0] + mp[1]*mn[1] + mp[2]*mn[2])/(TMath::Sqrt(mp[0]*mp[0] + mp[1]*mp[1] + mp[2]*mp[2])*TMath::Sqrt(mn[0]*mn[0] + mn[1]*mn[1] + mn[2]*mn[2])));
  
  return TMath::Abs(openAngle);
}

//_________________________________________________
Float_t AliHFEV0pid::PsiPair(AliESDv0 *v0) {
  //
  // Angle between daughter momentum plane and plane 
  // 
  Float_t magField = fInputEvent->GetMagneticField();

  Int_t pIndex = v0->GetPindex();
  Int_t nIndex = v0->GetNindex();

  AliESDtrack* daughter[2];

  daughter[0] = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(pIndex));
  daughter[1] = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(nIndex));

  Double_t x, y, z;
  v0->GetXYZ(x,y,z);//Reconstructed coordinates of V0; to be replaced by Markus Rammler's method in case of conversions!
  
  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0};
  

  v0->GetNPxPyPz(mn[0],mn[1],mn[2]);//reconstructed cartesian momentum components of negative daughter;
  v0->GetPPxPyPz(mp[0],mp[1],mp[2]);//reconstructed cartesian momentum components of positive daughter; 


  Double_t deltat = 1.;
  deltat = TMath::ATan(mp[2]/(TMath::Sqrt(mp[0]*mp[0] + mp[1]*mp[1])+1.e-13)) -  TMath::ATan(mn[2]/(TMath::Sqrt(mn[0]*mn[0] + mn[1]*mn[1])+1.e-13));//difference of angles of the two daughter tracks with z-axis

  Double_t radiussum = TMath::Sqrt(x*x + y*y) + 50;//radius to which tracks shall be propagated

  Double_t momPosProp[3];
  Double_t momNegProp[3];
    
  AliExternalTrackParam pt(*daughter[0]), nt(*daughter[1]);
    
  Float_t psiPair = 4.;

  if(nt.PropagateTo(radiussum,magField) == 0)//propagate tracks to the outside
    psiPair =  -5.;
  if(pt.PropagateTo(radiussum,magField) == 0)
    psiPair = -5.;
  pt.GetPxPyPz(momPosProp);//Get momentum vectors of tracks after propagation
  nt.GetPxPyPz(momNegProp);
  
  Double_t pEle =
    TMath::Sqrt(momNegProp[0]*momNegProp[0]+momNegProp[1]*momNegProp[1]+momNegProp[2]*momNegProp[2]);//absolute momentum value of negative daughter
  Double_t pPos =
    TMath::Sqrt(momPosProp[0]*momPosProp[0]+momPosProp[1]*momPosProp[1]+momPosProp[2]*momPosProp[2]);//absolute momentum value of positive daughter
    
  Double_t scalarproduct =
    momPosProp[0]*momNegProp[0]+momPosProp[1]*momNegProp[1]+momPosProp[2]*momNegProp[2];//scalar product of propagated positive and negative daughters' momenta
    
  Double_t chipair = TMath::ACos(scalarproduct/(pEle*pPos));//Angle between propagated daughter tracks

  psiPair =  TMath::Abs(TMath::ASin(deltat/chipair));  

  return psiPair; 
}

//____________________________________________________________
AliHFEV0pid::AliHFEV0pidTrackIndex::AliHFEV0pidTrackIndex():
    fNElectrons(0)
  , fNPionsK0(0)
  , fNPionsL(0)
  , fNKaons(0)
  , fNProtons(0)
  , fIndexElectron(NULL)
  , fIndexPionK0(NULL)
  , fIndexPionL(NULL)
  , fIndexKaon(NULL)
  , fIndexProton(NULL)
{
  //
  // Default Constructor
  //
}

//____________________________________________________________
AliHFEV0pid::AliHFEV0pidTrackIndex::~AliHFEV0pidTrackIndex(){
  //
  // Destructor
  //
  if(fIndexElectron) delete[] fIndexElectron;
  if(fIndexPionK0) delete[] fIndexPionK0;
  if(fIndexPionL) delete[] fIndexPionL;
  if(fIndexProton) delete[] fIndexProton;
}

//____________________________________________________________
void AliHFEV0pid::AliHFEV0pidTrackIndex::Flush(){
  //
  // Reset containers
  //
  
  if(fIndexElectron) delete[] fIndexElectron;
  fIndexElectron = NULL;
  if(fIndexPionK0) delete[] fIndexPionK0;
  fIndexPionK0 = NULL;
  if(fIndexPionL) delete[] fIndexPionL;
  fIndexPionL = NULL;
  if(fIndexKaon) delete[] fIndexKaon;
  fIndexKaon = NULL;
  if(fIndexProton) delete[] fIndexProton;
  fIndexProton = NULL;

  fNElectrons = 0;
  fNPionsK0 = 0;
  fNPionsL = 0;
  fNKaons = 0;
  fNProtons = 0;
}

//____________________________________________________________
void AliHFEV0pid::AliHFEV0pidTrackIndex::Init(Int_t capacity){
  //
  // Initialize container
  //
  fIndexElectron = new Int_t[capacity];
  fIndexPionK0 = new Int_t[capacity];
  fIndexPionL = new Int_t[capacity];
  fIndexProton = new Int_t[capacity];
}

//____________________________________________________________
void AliHFEV0pid::AliHFEV0pidTrackIndex::Add(Int_t index, Int_t species){
  //
  // Add new index to the list of identified particles
  //
  switch(species){
    case AliHFEV0pid::kRecoElectron:
      fIndexElectron[fNElectrons++] = index;
      break;
    case AliHFEV0pid::kRecoPionK0:
      fIndexPionK0[fNPionsK0++] = index;
      break;
    case AliHFEV0pid::kRecoPionL:
      fIndexPionL[fNPionsL++] = index;
      break;
    case AliHFEV0pid::kRecoProton:
      fIndexProton[fNProtons++] = index;
      break;
  };
}

//____________________________________________________________
Bool_t AliHFEV0pid::AliHFEV0pidTrackIndex::Find(Int_t index, Int_t species) const {
  //
  // Find track index in the specific sample of particles
  //

  Int_t *container = NULL; Int_t n = 0;
  switch(species){
  case AliHFEV0pid::kRecoElectron:
    container = fIndexElectron;
    n = fNElectrons;
    break;
  case AliHFEV0pid::kRecoPionK0:
    container = fIndexPionK0;
    n = fNPionsK0;
    break;
  case AliHFEV0pid::kRecoPionL:
    container = fIndexPionL;
    n = fNPionsL;
    break;
  case AliHFEV0pid::kRecoProton:
    container = fIndexProton;
    n = fNProtons;
    break;
  }
  if(!container) return kFALSE;
  if(n == 0) return kFALSE;
  Bool_t found = kFALSE;
  for(Int_t i = 0; i < n; i++){
    if(container[i] == index){
      found = kTRUE;
      break;
    }
  }
  return found;
}

//____________________________________________________________
Bool_t AliHFEV0pid::AliHFEV0pidTrackIndex::Find(Int_t index) const {
  // 
  // Find index in all samples
  //
  if(Find(index, AliHFEV0pid::kRecoElectron)) return kTRUE;
  else if(Find(index, AliHFEV0pid::kRecoPionK0)) return kTRUE;
  else if(Find(index, AliHFEV0pid::kRecoPionL)) return kTRUE;
  else return Find(index, AliHFEV0pid::kRecoProton);
}

//____________________________________________________________
AliHFEV0pid::AliHFELambdaInf::AliHFELambdaInf(AliKFParticle *mother, AliKFVertex * const primaryVertex):
  fInvariantMass(0.),
  fChi2NDF(0.),
  fDistVtx(0.)
{
  //
  // Constructor
  // Fill infos
  //
  const Double_t kL0mass=TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();  // PDG lambda mass
  fInvariantMass = mother->GetMass();
  // distance frm the privary vertex and mass difference for the (A) case
  fDistVtx = mother->GetDistanceFromVertex(*primaryVertex);
  // apply constraints for the fit
  mother->SetMassConstraint(kL0mass, 0.);
  mother->SetProductionVertex(*primaryVertex);
  fChi2NDF = mother->GetChi2()/mother->GetNDF();
}
//____________________________________________________________
AliHFEV0pid::AliHFELambdaInf::~AliHFELambdaInf(){
  //
  // destructor
  //
}
//____________________________________________________________
// DEBUG
//____________________________________________________________
Double_t AliHFEV0pid::GetEffMass(AliESDv0 *v0, UInt_t p1, UInt_t p2) const{
  //
  // TEMPORARY - this function should become obsolete with v4-18-Rev-10 or 11
  // calculate effective mass
  //
  const Double_t kpmass[5] = {TDatabasePDG::Instance()->GetParticle(kElectron)->Mass(),
			     TDatabasePDG::Instance()->GetParticle(kMuonMinus)->Mass(),
			     TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass(),
			     TDatabasePDG::Instance()->GetParticle(kKPlus)->Mass(),
			      TDatabasePDG::Instance()->GetParticle(kProton)->Mass()}; // float
  if (p1>4) return -1;
  if (p2>4) return -1;
  Double_t mass1 = kpmass[p1]; // float
  Double_t mass2 = kpmass[p2]; // float  

  Double_t pMom[3];
  Double_t nMom[3];
  
  v0->GetPPxPyPz(pMom[0], pMom[1], pMom[2]);
  v0->GetNPxPyPz(nMom[0], nMom[1], nMom[2]);
  

  const Double_t *m1 = pMom;
  const Double_t *m2 = nMom;
  //
  //if (fRP[p1]+fRM[p2]<fRP[p2]+fRM[p1]){
  //  m1 = fPM;
  //  m2 = fPP;
  //}
  //
  Double_t e1    = TMath::Sqrt(mass1*mass1+
                              m1[0]*m1[0]+
                              m1[1]*m1[1]+
                              m1[2]*m1[2]); // float
  Double_t e2    = TMath::Sqrt(mass2*mass2+
                              m2[0]*m2[0]+
                              m2[1]*m2[1]+
                              m2[2]*m2[2]); // float
  Double_t mass =  
    (m2[0]+m1[0])*(m2[0]+m1[0])+
    (m2[1]+m1[1])*(m2[1]+m1[1])+
    (m2[2]+m1[2])*(m2[2]+m1[2]); // float

  
  mass = (e1+e2)*(e1+e2)-mass;
  //if(mass < 0.00001){
  //  printf("-D: mass: %f\n", mass);
  // }
  mass = TMath::Sqrt(mass);
  return mass;
}
