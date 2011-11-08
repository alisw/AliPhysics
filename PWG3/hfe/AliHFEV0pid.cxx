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
#include <TObjArray.h>

#include "AliAODEvent.h"
#include "AliAODv0.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliKFVertex.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliMCEvent.h"

#include "AliHFEV0cuts.h"
#include "AliHFEV0info.h"
#include "AliHFEcollection.h"

#include "AliHFEV0pid.h"
ClassImp(AliHFEV0pid)

AliHFEV0pid::AliHFEV0pid():
  TNamed()
  , fInputEvent(NULL)
  , fNtracks(0)
  , fMCEvent(NULL)
  , fMCon(kFALSE)
  , fPrimaryVertex(NULL)
  , fElectrons(NULL)
  , fPionsK0(NULL)
  , fPionsL(NULL)
  , fKaons(NULL)
  , fProtons(NULL)
  , fGammas(NULL)
  , fK0s(NULL)
  , fLambdas(NULL)
  , fAntiLambdas(NULL)
  , fIndices(NULL)
  , fQA(NULL)
  , fV0cuts(NULL)
  , fOutput(NULL)
  , fDestBits(0)

{
  //
  // Default constructor
  //
}
//____________________________________________________________
AliHFEV0pid::AliHFEV0pid(const char *name):
  TNamed(name, "")
  , fInputEvent(NULL)
  , fNtracks(0)
  , fMCEvent(NULL)
  , fMCon(kFALSE)
  , fPrimaryVertex(NULL)
  , fElectrons(NULL)
  , fPionsK0(NULL)
  , fPionsL(NULL)
  , fKaons(NULL)
  , fProtons(NULL)
  , fGammas(NULL)
  , fK0s(NULL)
  , fLambdas(NULL)
  , fAntiLambdas(NULL)
  , fIndices(NULL)
  , fQA(NULL)
  , fV0cuts(NULL)
  , fOutput(NULL)
  , fDestBits(0)
{
  //
  // Standard constructor
  //
  fElectrons = new TObjArray();
  fPionsK0 = new TObjArray();
  fPionsL = new TObjArray();
  fKaons = new TObjArray();
  fProtons = new TObjArray();

  fElectrons->SetOwner();
  fPionsK0->SetOwner();
  fProtons->SetOwner();
  fPionsL->SetOwner();
  fKaons->SetOwner();
  
  fGammas = new TObjArray();
  fK0s = new TObjArray();
  fLambdas = new TObjArray();
  fAntiLambdas = new TObjArray();

  fIndices = new AliHFEV0pidTrackIndex();
  
}
//____________________________________________________________

AliHFEV0pid::~AliHFEV0pid(){
  //
  // Destructor
  // Remove Containers
  //
  if(fElectrons) delete fElectrons;
  if(fPionsK0) delete fPionsK0;
  if(fPionsL) delete fPionsL;
  if(fKaons) delete fKaons;
  if(fProtons) delete fProtons;

  if(fGammas) delete fGammas;
  if(fK0s) delete fK0s;
  if(fLambdas) delete fLambdas;
  if(fAntiLambdas) delete fAntiLambdas;

  if(fIndices) delete fIndices;
  if(fV0cuts) delete fV0cuts;

  if(TESTBIT(fDestBits, 1)){
    if(fQA) delete fQA;
    if(fOutput) delete fOutput;
  }
}

//____________________________________________________________
void AliHFEV0pid::InitQA(){
  //
  // Initialize QA histograms
  //
  
  memset(&fDestBits, 0, sizeof(fDestBits));
  SETBIT(fDestBits, 1);

  fOutput = new TList();
  fOutput->SetOwner();

  fV0cuts = new AliHFEV0cuts();
  fV0cuts->Init("V0cuts");

  if(!fQA){
    fQA = new AliHFEcollection("v0pidQA", "QA histograms for V0 selection");

    fQA->CreateTH1F("h_nV0s", "No. of found and accepted V0s", 5, -0.5, 4.5);

    // QA histograms for invariant mass
    fQA->CreateTH1F("h_InvMassGamma", "Gamma invariant mass; inv mass [GeV/c^{2}]; counts", 100, 0, 0.25);
    fQA->CreateTH1F("h_InvMassK0s", "K0s invariant mass; inv mass [GeV/c^{2}]; counts", 200, 0.4, 0.65);
    fQA->CreateTH1F("h_InvMassL", "Lambda invariant mass; inv mass [GeV/c^{2}]; counts", 100, 1.05, 1.15);
    fQA->CreateTH1F("h_InvMassAL", "Lambda invariant mass; inv mass [GeV/c^{2}]; counts", 100, 1.05, 1.15);
    
    // QA histograms for p distribution (of the daughters)
    fQA->CreateTH1F("h_P_electron", "P distribution of the gamma electrons; p (GeV/c); counts", 100, 0.1, 10);
    fQA->CreateTH1F("h_P_K0pion", "P distribution of the K0 pions; p (GeV/c); counts", 100, 0.1, 10);
    fQA->CreateTH1F("h_P_Lpion", "P distribution of the Lambda pions; p (GeV/c); counts", 100, 0.1, 10);
    fQA->CreateTH1F("h_P_Lproton", "P distribution of the Lambda protons; p (GeV/c); counts", 100, 0.1, 10);

    // QA pt of the V0
    fQA->CreateTH1F("h_Pt_Gamma", "Pt of the gamma conversion; p_{T} (GeV/c); counts", 100, 0, 10);
    fQA->CreateTH1F("h_Pt_K0", "Pt of the K0; p_{T} (GeV/c); counts", 100, 0, 10);
    fQA->CreateTH1F("h_Pt_L", "Pt of the Lambda; p_{T} (GeV/c); counts", 100, 0, 10);    
    fQA->CreateTH1F("h_Pt_AL", "Pt of the Lambda; p_{T} (GeV/c); counts", 100, 0, 10);    
    
    // Armenteros plot V0 preselection
    fQA->CreateTH2F("h_AP_all_V0s", "armenteros plot for all V0 candidates; #alpha; Q_{T}", 200, -1, 1, 200, 0, 0.25);
    fQA->CreateTH2F("h_AP_selected_V0s", "armenteros plot for all V0 candidates; #alpha; Q_{T}", 200, -1, 1, 200, 0, 0.25);

    //
    // !!! MC plots !!!
    // 
    fQA->CreateTH2F("h_AP_MC_all_V0s", "armenteros plot for all MC tagged V0s; #alpha; Q_{T}", 200, -1, 1, 200, 0, 0.25);
    fQA->CreateTH2F("h_AP_MC_Gamma", "armenteros plot for MC tagged Gammas; #alpha; Q_{T}",  200, -1, 1, 200, 0, 0.25);
    fQA->CreateTH2F("h_AP_MC_K0", "armenteros plot for MC tagged K0s; #alpha; Q_{T}",  200, -1, 1, 200, 0, 0.25);
    fQA->CreateTH2F("h_AP_MC_Lambda", "armenteros plot for MC tagged Lambdas; #alpha; Q_{T}",  200, -1, 1, 200, 0, 0.25);
    fQA->CreateTH2F("h_AP_MC_ALambda", "armenteros plot for MC tagged A-Lambdass; #alpha; Q_{T}",  200, -1, 1, 200, 0, 0.25);
 

   
    // armenteros plots for different V0 momenta - MC signal only
    fQA->CreateTH2Fvector1(12, "h_AP_MC_Gamma_p", "armenteros plot for MC tagged Gammas; #alpha; Q_{T}",  200, -1., 1., 200, 0., 0.25);
    fQA->CreateTH2Fvector1(12, "h_AP_MC_K0_p", "armenteros plot for MC tagged K0s; #alpha; Q_{T}",  200, -1., 1., 200, 0., 0.25);
    fQA->CreateTH2Fvector1(12, "h_AP_MC_Lambda_p", "armenteros plot for MC tagged Lambdas; #alpha; Q_{T}",  200, -1., 1., 200, 0., 0.25);
   //

    
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
  fNtracks = fInputEvent->GetNumberOfTracks();
  fIndices->Init(fInputEvent->GetNumberOfV0s() * 2);
  fPrimaryVertex = new AliKFVertex(*(fInputEvent->GetPrimaryVertex()));
  if(!fPrimaryVertex) return;
  fV0cuts->SetInputEvent(fInputEvent);
  fV0cuts->SetPrimaryVertex(fPrimaryVertex);
  if(fMCEvent) fV0cuts->SetMCEvent(fMCEvent);
  Int_t check[fNtracks];
  memset(check, 0, sizeof(Int_t)*fNtracks);
  Int_t v0status = 0;

  //BenchmarkV0finder();

  for(Int_t iv0 = 0; iv0 < fInputEvent->GetNumberOfV0s(); iv0++){
    if(!TString(fInputEvent->IsA()->GetName()).CompareTo("AliESDEvent")){
      // case ESD
      SetESDanalysis();
      AliESDv0 *esdV0 = (static_cast<AliESDEvent *>(fInputEvent))->GetV0(iv0);
      if(!esdV0) continue;
      if(!esdV0->GetOnFlyStatus()) continue; // Take only V0s from the On-the-fly v0 finder
      v0status = ProcessV0(esdV0);
    } else {
      // case AOD
      SetAODanalysis();
      AliAODv0 *aodV0 = (static_cast<AliAODEvent *>(fInputEvent))->GetV0(iv0);
      if(aodV0->GetOnFlyStatus()) continue; // Take only V0s from the On-the-fly v0 finder
      v0status = ProcessV0(aodV0);
      if(AliHFEV0cuts::kUndef != v0status){
      }
    }
    switch(v0status){
    case AliHFEV0cuts::kRecoGamma: nGamma++; break;
    case AliHFEV0cuts::kRecoK0: nK0s++; break;
    case AliHFEV0cuts::kRecoPhi: nPhi++; break;  
    case AliHFEV0cuts::kRecoLambda: nLambda++; break;
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
  if(!v0)  return AliHFEV0cuts::kUndef;
  // CHECK
  AliESDv0* esdV0 =  dynamic_cast<AliESDv0 *>(v0);
  if( ! esdV0 ) {
    AliError("Unexpected v0 type.");
    return AliHFEV0cuts::kUndef;
  }  
  AliVTrack* daughter[2];
  daughter[0] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(esdV0->GetPindex()));
  daughter[1] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(esdV0->GetNindex()));
  if(!daughter[0] || !daughter[1]) return AliHFEV0cuts::kUndef;

  if(fMCEvent != NULL) fMCon = kTRUE;
  //printf("-D: fMCEvent %x, fMCon: %i\n", fMCEvent, fMCon);


  Int_t dMC[2] = {-1, -1};
  Int_t idMC = AliHFEV0cuts::kUndef; 

  if(IsESDanalysis()){
    if(fMCon){
      idMC = IdentifyV0(v0, dMC);
      //printf("--D: V0 pid: %i, P: %i, N: %i\n", id, d[0], d[1]);
      fV0cuts->SetCurrentV0id(idMC);
      fV0cuts->SetDaughtersID(dMC);
    }
    // check common single track cuts
    for(Int_t i=0; i<2; ++i){
      if(!fV0cuts->TrackCutsCommon(static_cast<AliESDtrack*>(daughter[i]))) return AliHFEV0cuts::kUndef;
    }
    // check commom V0 cuts
    // CHECK
    if(!fV0cuts->V0CutsCommon(esdV0)) return AliHFEV0cuts::kUndef;
  }

  // preselect the V0 candidates based on the Armenteros plot
  Int_t id = PreselectV0(esdV0, idMC);
  // store the resutls
  if(AliHFEV0cuts::kRecoGamma == id && IsGammaConv(v0)){
    fQA->Fill("h_nV0s", AliHFEV0cuts::kRecoGamma);
    return AliHFEV0cuts::kRecoGamma;
  }
  else if(AliHFEV0cuts::kRecoK0 == id && IsK0s(v0)){
    fQA->Fill("h_nV0s", AliHFEV0cuts::kRecoK0);
    return AliHFEV0cuts::kRecoK0;
  }
  else if(AliHFEV0cuts::kRecoLambda == TMath::Abs(id) && IsLambda(v0)){
    fQA->Fill("h_nV0s", AliHFEV0cuts::kRecoLambda);    
    return AliHFEV0cuts::kRecoLambda;
  }
  else return AliHFEV0cuts::kUndef; 


    
}
//____________________________________________________________
void AliHFEV0pid::Flush(){
  //
  // Clear the Lists
  //
  AliDebug(1, "Flushing containers");
  fProtons->Delete();
  fPionsK0->Delete();
  fPionsL->Delete();
  fElectrons->Delete();
  fIndices->Flush();
}
//____________________________________________________________
Int_t AliHFEV0pid::PreselectV0(AliESDv0 * const v0, Int_t idMC){
  //
  // Based on Armenteros plot preselet the possible identity of the V0 candidate
  //

  if(!v0) return -1;

  // momentum dependent armenteros plots
  ArmenterosPlotMC(v0, idMC);

  // comupte the armenteros variables
  Float_t ar[2];
  fV0cuts->Armenteros(v0, ar);
  // for clarity
  const Float_t alpha = ar[0];
  const Float_t qt = ar[1];
  //printf(" -D: Alpha: %f, QT: %f \n", alpha, qt);

  if(TMath::Abs(alpha) > 1) return AliHFEV0cuts::kUndef;

  fQA->Fill("h_AP_all_V0s", alpha, qt);

  // fill a few MC tagged histograms - AP plots
  if(fMCEvent){
    switch(idMC){
    case AliHFEV0cuts::kRecoGamma :
      fQA->Fill("h_AP_MC_all_V0s", alpha, qt);
      fQA->Fill("h_AP_MC_Gamma", alpha, qt);
      break;
    case AliHFEV0cuts::kRecoK0 :
      fQA->Fill("h_AP_MC_all_V0s", alpha, qt);
    fQA->Fill("h_AP_MC_K0", alpha, qt);
      break;
    case AliHFEV0cuts::kRecoLambda :
      fQA->Fill("h_AP_MC_all_V0s", alpha, qt);
      fQA->Fill("h_AP_MC_Lambda", alpha, qt);
      break;
    case AliHFEV0cuts::kRecoALambda :
      fQA->Fill("h_AP_MC_all_V0s", alpha, qt);
      fQA->Fill("h_AP_MC_ALambda", alpha, qt);
      break;
    }
  }

  // Gamma cuts
  const Double_t cutAlphaG = 0.35; 
  const Double_t cutQTG = 0.05;
  const Double_t cutAlphaG2[2] = {0.6, 0.8};
  const Double_t cutQTG2 = 0.04;

  // K0 cuts
  const Float_t cutQTK0[2] = {0.1075, 0.215};
  const Float_t cutAPK0[2] = {0.199, 0.8};   // parameters for curved QT cut
  
  // Lambda & A-Lambda cuts
  const Float_t cutQTL = 0.03;
  const Float_t cutAlphaL[2] = {0.35, 0.7};
  const Float_t cutAlphaAL[2] = {-0.7,  -0.35};
  const Float_t cutAPL[3] = {0.107, -0.69, 0.5};  // parameters fir curved QT cut


  // Check for Gamma candidates
  if(qt < cutQTG){
    if( (TMath::Abs(alpha) < cutAlphaG) ){
      fQA->Fill("h_AP_selected_V0s", alpha, qt);
      return  AliHFEV0cuts::kRecoGamma;
    }
  }
  // additional region - should help high pT gammas
  if(qt < cutQTG2){
    if( (TMath::Abs(alpha) > cutAlphaG2[0]) &&  (TMath::Abs(alpha) < cutAlphaG2[1]) ){
      fQA->Fill("h_AP_selected_V0s", alpha, qt);
      return  AliHFEV0cuts::kRecoGamma;
    }
  }


  // Check for K0 candidates
  Float_t q = cutAPK0[0] * TMath::Sqrt(TMath::Abs(1 - alpha*alpha/(cutAPK0[1]*cutAPK0[1])));
  if( (qt > cutQTK0[0]) && (qt < cutQTK0[1]) && (qt > q) ){
    fQA->Fill("h_AP_selected_V0s", alpha, qt);
    return AliHFEV0cuts::kRecoK0;
  }
  
  if( (alpha > 0) && (alpha > cutAlphaL[0])  && (alpha < cutAlphaL[1]) && (qt > cutQTL)){
    q = cutAPL[0] * TMath::Sqrt(1 - ( (alpha + cutAPL[1]) * (alpha + cutAPL[1]))  / (cutAPL[2]*cutAPL[2]) );
    if( qt < q  ){
      fQA->Fill("h_AP_selected_V0s", alpha, qt);
      return AliHFEV0cuts::kRecoLambda;
    }
  }

  // Check for A-Lambda candidates
  if( (alpha < 0) && (alpha > cutAlphaAL[0]) && (alpha < cutAlphaAL[1]) && (qt > cutQTL)){
    q = cutAPL[0] * TMath::Sqrt(1 - ( (alpha - cutAPL[1]) * (alpha - cutAPL[1]) ) / (cutAPL[2]*cutAPL[2]) );
    if( qt < q ){
      fQA->Fill("h_AP_selected_V0s", alpha, qt);
      return AliHFEV0cuts::kRecoLambda;
    }
  }
  
  return AliHFEV0cuts::kUndef;

}
//____________________________________________________________
Bool_t AliHFEV0pid::IsGammaConv(TObject *v0){
  //
  // Identify Gamma
  //

  if(!v0) return kFALSE;

  AliVTrack* daughter[2];
  Int_t pIndex = 0, nIndex = 0;
  Double_t invMass = 0.;
  Double_t mPt = 0.;
  Int_t v0id = -1;
  if(IsESDanalysis()){
    // ESD - cut V0
    AliESDv0 *esdV0 = static_cast<AliESDv0 *>(v0);
    v0id = esdV0->GetLabel();
    // apply FULL gamma cuts
    if(!fV0cuts->GammaCuts(esdV0)) return kFALSE;
    invMass = esdV0->GetEffMass(AliPID::kElectron, AliPID::kElectron);
    pIndex = esdV0->GetPindex();
    nIndex = esdV0->GetNindex();
    mPt = esdV0->Pt();
  } else {
    // AOD Analysis - not possible to cut
    AliAODv0 *aodV0 = static_cast<AliAODv0 *>(v0);
    v0id = aodV0->GetID();
    pIndex = aodV0->GetPosID();
    nIndex = aodV0->GetNegID();
    invMass = aodV0->InvMass2Prongs(0, 1, kElectron, kElectron);
  }
  daughter[0] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(pIndex));
  daughter[1] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(nIndex));
  if(!daughter[0] || !daughter[1]) return kFALSE;
 
  // DEBUG
  AliDebug(1, Form("Gamma identified, daughter IDs: %d,%d", daughter[0]->GetID(), daughter[1]->GetID()));

  // AFTER all gamma cuts
  fQA->Fill("h_Pt_Gamma", mPt);
  fQA->Fill("h_InvMassGamma", invMass);

  // Identified gamma - store tracks in the electron containers
  if(!fIndices->Find(daughter[0]->GetID())){
    AliDebug(1, Form("Gamma identified, daughter IDs: %d,%d", daughter[0]->GetID(), daughter[1]->GetID()));    
    fElectrons->Add(new AliHFEV0info(daughter[0], daughter[1]->GetID(), v0id));
    fIndices->Add(daughter[0]->GetID(), AliHFEV0cuts::kRecoElectron);
  }
  if(!fIndices->Find(daughter[1]->GetID())){
    AliDebug(1, Form("Gamma identified, daughter IDs: %d,%d", daughter[1]->GetID(), daughter[1]->GetID()));
    fElectrons->Add(new AliHFEV0info(daughter[1], daughter[0]->GetID(), v0id));
    fIndices->Add(daughter[1]->GetID(), AliHFEV0cuts::kRecoElectron);
  }
  fGammas->Add(v0);
  
  return kTRUE;
}
//____________________________________________________________
Bool_t AliHFEV0pid::IsK0s(TObject *v0){
  //
  // Identify K0s
  //

  if(!v0) return kFALSE;

  AliVTrack* daughter[2];
  Int_t pIndex = 0, nIndex = 0;
  Int_t v0id = -1;
  Double_t invMass = 0.;
  Double_t mPt = 0.;
  if(IsESDanalysis()){
    // ESD - cut V0
    AliESDv0 *esdV0 = static_cast<AliESDv0 *>(v0);
    if(!fV0cuts->K0Cuts(esdV0)) return kFALSE;
    v0id = esdV0->GetLabel();
    pIndex = esdV0->GetPindex();
    nIndex = esdV0->GetNindex();
    invMass = esdV0->GetEffMass(AliPID::kPion, AliPID::kPion);
    mPt = esdV0->Pt();
  } else {
    // AOD Analysis - not possible to cut
    AliAODv0 *aodV0 = static_cast<AliAODv0 *>(v0);
    aodV0->GetID();
    pIndex = aodV0->GetPosID();
    nIndex = aodV0->GetNegID();
    invMass = aodV0->MassK0Short();
  }
  daughter[0] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(pIndex));
  daughter[1] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(nIndex));
  if(!daughter[0] || !daughter[1]) return kFALSE;

  fQA->Fill("h_Pt_K0", mPt);
  fQA->Fill("h_InvMassK0s", invMass);

  AliDebug(1, Form("K0 identified, daughter IDs: %d,%d", daughter[0]->GetID(), daughter[1]->GetID()));

  // AFTER all K0 cuts
  // Identified gamma - store tracks in the electron containers
  if(!fIndices->Find(daughter[0]->GetID())){
    AliDebug(1, Form("Adding K0 Pion track with ID %d", daughter[0]->GetID()));
    fPionsK0->Add(new AliHFEV0info(daughter[0], daughter[1]->GetID(), v0id));
    fIndices->Add(daughter[0]->GetID(), AliHFEV0cuts::kRecoPionK0);
  }
  if(!fIndices->Find(daughter[1]->GetID())){
    AliDebug(1, Form("Adding K0 Pion track with ID %d", daughter[1]->GetID()));
    fPionsK0->Add(new AliHFEV0info(daughter[1], daughter[0]->GetID(), v0id));
    fIndices->Add(daughter[1]->GetID(), AliHFEV0cuts::kRecoPionK0);
  }
  fK0s->Add(v0);
  return kTRUE; 
}

//____________________________________________________________
Bool_t AliHFEV0pid::IsPhi(const TObject *v0) const {
  //
  // Identify Phi - very preliminary - requires diffrent approach (V0 fnder is not effective)
  //

  //const Double_t kPhiMass=TDatabasePDG::Instance()->GetParticle(333)->Mass();  // PDG phi mass
  //AliVTrack* daughter[2];
  //Double_t invMass = 0.;

  if(!v0) return kFALSE;
 
  //Int_t pIndex = 0, nIndex = 0;
  if(IsESDanalysis()){
    // ESD - cut V0
    //AliESDv0 *esdV0 = static_cast<AliESDv0 *>(v0);
    //pIndex = esdV0->GetPindex();
    //nIndex = esdV0->GetNindex();
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

  if(!v0) return kFALSE;

  AliVTrack* daughter[2];
  Int_t pIndex = 0, nIndex = 0;
  Double_t invMass = 0.;
  Bool_t isLambda = kTRUE; // Lambda - kTRUE, Anti Lambda - kFALSE
  Double_t mPt = 0.;
  Int_t v0id = -1;
  if(IsESDanalysis()){
    // ESD - cut V0
    AliESDv0 *esdV0 = static_cast<AliESDv0 *>(v0);
    v0id = esdV0->GetLabel();
    if(!fV0cuts->LambdaCuts(esdV0,isLambda)) return kFALSE; 
    mPt = esdV0->Pt();
    if(fV0cuts->CheckSigns(esdV0)){
      pIndex = esdV0->GetPindex();
      nIndex = esdV0->GetNindex();
    }
    else{
      pIndex = esdV0->GetNindex();
      nIndex = esdV0->GetPindex();      
    }
  } else {
    // PRELIMINARY - !!!
    // AOD Analysis - not possible to cut
    
    // again - two cases as above
    AliAODv0 *aodV0 = static_cast<AliAODv0 *>(v0);
    v0id = aodV0->GetID();
    pIndex = aodV0->GetPosID();
    nIndex = aodV0->GetNegID();
    invMass = aodV0->MassLambda();
  } 
  
  daughter[0] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(pIndex));
  daughter[1] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(nIndex));
  if(!daughter[0] || !daughter[1]) return kFALSE;
  
  // lambda
  if(isLambda){
    fQA->Fill("h_Pt_L", mPt);
    fQA->Fill("h_InvMassL", invMass);

    if(!fIndices->Find(daughter[0]->GetID())){
      fProtons->Add(new AliHFEV0info(daughter[0], daughter[1]->GetID(), v0id));
      fIndices->Add(daughter[0]->GetID(), AliHFEV0cuts::kRecoProton);
    }
    if(!fIndices->Find(daughter[1]->GetID())){
      fPionsL->Add(new AliHFEV0info(daughter[1], daughter[0]->GetID(), v0id));
      fIndices->Add(daughter[1]->GetID(), AliHFEV0cuts::kRecoPionL);
    }
  }
  // antilambda
  else{
    fQA->Fill("h_Pt_AL", mPt);
    fQA->Fill("h_InvMassAL", invMass);
    if(!fIndices->Find(daughter [1]->GetID())){
      fProtons->Add(new AliHFEV0info(daughter[1], daughter[0]->GetID(), v0id));
      fIndices->Add(daughter[1]->GetID(), AliHFEV0cuts::kRecoProton);
    }
    if(!fIndices->Find(daughter [0]->GetID())){
      fPionsL->Add(new AliHFEV0info(daughter[0], daughter[1]->GetID(), v0id));
      fIndices->Add(daughter [0]->GetID(), AliHFEV0cuts::kRecoPionL);
    }
  }
  if(isLambda) fLambdas->Add(v0);
  else fAntiLambdas->Add(v0);

  return kTRUE;
}

//____________________________________________________________
Int_t AliHFEV0pid::IdentifyV0(TObject *esdV0, Int_t d[2]){
  //
  // for MC only, returns the V0 Id
  //

  //
  // be carefull about changing the return values - they are used later selectively
  // In particulra "-2" means that identity of either of the daughters could not be
  // estimated
  //

  AliESDv0 *v0 = dynamic_cast<AliESDv0 *>(esdV0);
  
  if(!v0) return -1;
  AliESDtrack* dN, *dP; 
  Int_t iN, iP;
  iN = iP = -1;
  iP = v0->GetPindex();
  iN = v0->GetNindex();
  if(iN < 0 || iP < 0) return -1;
  if(iN >= fNtracks || iP >= fNtracks) return -1;
  dP = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(iP));
  dN = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(iN));  
  if(!dN || !dP) return -1;

  // as of 26/10/2010
  // there is still a problem with wrong assignment of positive and negative
  // V0 daughter in V0 finder - a check is necessary
  // if the V0 daughters are miss-assigned - swap their labels
  Bool_t sign = fV0cuts->CheckSigns(v0);

  // get the MC labels
  Int_t lN, lP;
  if(sign){
    lN = dN->GetLabel();
    lP = dP->GetLabel();
  }
  else{
    lP = dN->GetLabel();
    lN = dP->GetLabel();
  }
  if(lN < 0 || lP < 0) return -2;
  // get the associated MC particles
  AliMCParticle *mcP, *mcN;
  mcP = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(lP));
  mcN = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(lN));
  if(!mcP || !mcN) return -2;
  
  // identify the daughter tracks and their mothers
  Int_t pdgP, pdgN;
  pdgP = TMath::Abs(mcP->PdgCode());
  pdgN = TMath::Abs(mcN->PdgCode());
  // store the daughter ID for later use
  d[0] = pdgP;
  d[1] = pdgN;
  //printf(" -D: pdgP: %i, pdgN: %i\n", pdgP, pdgN);
  // lablel of the mother particle
  // -1 may mean it was a primary particle
  Int_t lPm, lNm;
  lPm = mcP->GetMother();
  lNm = mcN->GetMother();
  if(-1==lPm || -1==lNm) return -3;

  // return if mothers are not the same particle
  if(lPm != lNm) return -3;
  // get MC mother particle - now we need only 1
  AliMCParticle *m = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(lPm));
  if(!m) return -2;
  // mother PDG
  Int_t pdgM = m->PdgCode();

  //   if(3122 == TMath::Abs(pdgM)){
  //     printf("-D: v0 signs: %i\n", fV0cuts->CheckSigns(v0));
  //     printf("-D: pdgM: %i, pdgN: %i, pdgP: %i \n", pdgM, pdgN, pdgP);
  //   }
  
 
  // now check the mother and daughters identity
  if(22 == TMath::Abs(pdgM) && 11 == pdgN && 11 == pdgP) return AliHFEV0cuts::kRecoGamma;
  if(310 == TMath::Abs(pdgM) && 211 == pdgN && 211 == pdgP) return AliHFEV0cuts::kRecoK0;
  if(-3122 == pdgM && 2212 == pdgN && 211 == pdgP) return AliHFEV0cuts::kRecoALambda;
  if(3122 == pdgM && 211 == pdgN && 2212 == pdgP) return AliHFEV0cuts::kRecoLambda;
    
  return AliHFEV0cuts::kUndef;

}
//____________________________________________________________
void AliHFEV0pid::BenchmarkV0finder(){
  //
  // produce histograms for all findable V0s that are
  // were selected byt the (oline) V0 finder and can
  // be used to estimate the efficiency of teh V0 cuts
  //

  for(Int_t iv0 = 0; iv0 < fInputEvent->GetNumberOfV0s(); iv0++){
    AliESDv0 *esdV0 = (static_cast<AliESDEvent *>(fInputEvent))->GetV0(iv0);
    if(!esdV0) continue;
    if(!esdV0->GetOnFlyStatus()) continue; // Take only V0s from the On-the-fly v0 finder
    // indetify the V0 candidate
    Int_t idV0 = AliHFEV0cuts::kUndef;
    Int_t idD[2] = {-1, -1};
    idV0 = IdentifyV0(esdV0, idD);
  }
}
//____________________________________________________________
void   AliHFEV0pid::ArmenterosPlotMC(AliESDv0 * const v0, Int_t idMC){
  //
  // Armenteros plots as a function of Mohter Momentum
  //
  //const Float_t minP = 0.1;
  //const Float_t maxP = 10.;
  // approx log bins - over the 0.1 - 10 GeV/c
  const Float_t bins[13] = {0.1, 0.1468, 0.2154, 0.3162, 0.4642, 0.6813, 1.0, 1.4678, 2.1544, 3.1623, 4.6416, 6.8129, 10.0};
  
  Float_t ar[2];
  fV0cuts->Armenteros(v0, ar);
  Float_t p = v0->P();
 
  if( (p <=  bins[0]) || (p >= bins[12])) return;

  Int_t pBin = 0;
  Float_t tmp = bins[0];
  while(tmp < p){
    ++pBin;
    tmp = bins[pBin];
  }
  pBin--;

  if(AliHFEV0cuts::kRecoGamma == idMC) fQA->Fill("h_AP_MC_Gamma_p", pBin, ar[0], ar[1]);
  if(AliHFEV0cuts::kRecoK0 == idMC) fQA->Fill("h_AP_MC_K0_p", pBin, ar[0], ar[1]);
  if(AliHFEV0cuts::kRecoLambda == TMath::Abs(idMC)) fQA->Fill("h_AP_MC_Lambda_p", pBin, ar[0], ar[1]);
  
  

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
    case AliHFEV0cuts::kRecoElectron:
      fIndexElectron[fNElectrons++] = index;
      break;
    case AliHFEV0cuts::kRecoPionK0:
      fIndexPionK0[fNPionsK0++] = index;
      break;
    case AliHFEV0cuts::kRecoPionL:
      fIndexPionL[fNPionsL++] = index;
      break;
    case AliHFEV0cuts::kRecoProton:
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
  case AliHFEV0cuts::kRecoElectron:
    container = fIndexElectron;
    n = fNElectrons;
    break;
  case AliHFEV0cuts::kRecoPionK0:
    container = fIndexPionK0;
    n = fNPionsK0;
    break;
  case AliHFEV0cuts::kRecoPionL:
    container = fIndexPionL;
    n = fNPionsL;
    break;
  case AliHFEV0cuts::kRecoProton:
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
  if(Find(index, AliHFEV0cuts::kRecoElectron)) return kTRUE;
  else if(Find(index, AliHFEV0cuts::kRecoPionK0)) return kTRUE;
  else if(Find(index, AliHFEV0cuts::kRecoPionL)) return kTRUE;
  else return Find(index, AliHFEV0cuts::kRecoProton);
}

//____________________________________________________________
TList *AliHFEV0pid::GetListOfQAhistograms(){
  //
  // Getter for V0 PID QA histograms
  //

  CLRBIT(fDestBits, 1);

  TList *tmp = fV0cuts->GetList();
  tmp->SetName("V0cuts");
  fOutput->Add(tmp);
  if(fQA){
    tmp = 0x0;
    tmp = fQA->GetList();
    tmp->SetName("V0pid");
    fOutput->Add(tmp);
  } 
  tmp = 0x0;
  tmp = fV0cuts->GetListMC();
  tmp->SetName("V0cutsMC");
  //printf(" -D: adding MC V0 cuts stuff\n");
  fOutput->Add(tmp);
  
  return fOutput;
}
