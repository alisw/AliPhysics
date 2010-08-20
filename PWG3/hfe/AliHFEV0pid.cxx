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

#include "AliHFEV0cuts.h"
#include "AliHFEV0info.h"
#include "AliHFEcollection.h"

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
  , fGammas(NULL)
  , fK0s(NULL)
  , fLambdas(NULL)
  , fAntiLambdas(NULL)
  , fIndices(NULL)
  , fQA(NULL)
  , fV0cuts(NULL)
  , fOutput(NULL)
{
  //
  // Default constructor
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
  if(fQA) delete fQA;
  if(fV0cuts) delete fV0cuts;
  if(fOutput) delete fOutput;
}

//____________________________________________________________
void AliHFEV0pid::InitQA(){
  //
  // Initialize QA histograms
  //
  
  fOutput = new TList();

  fV0cuts = new AliHFEV0cuts();
  fV0cuts->Init("V0cuts");

  if(!fQA){
    fQA = new AliHFEcollection("v0pidQA", "QA histograms for V0 PID");

    fQA->CreateTH1F("h_nV0s", "No. of found and accepted V0s", 5, -0.5, 4.5);

    // QA histograms for invariant mass
    fQA->CreateTH1F("h_InvMassGamma", "Gamma invariant mass; inv mass [GeV/c^{2}]; counts", 100, 0, 0.25);
    fQA->CreateTH1F("h_InvMassK0s", "K0s invariant mass; inv mass [GeV/c^{2}]; counts", 200, 0.4, 0.65);
    fQA->CreateTH1F("h_InvMassLambda", "Lambda invariant mass; inv mass [GeV/c^{2}]; counts", 100, 1.05, 1.15);
    
    // QA histograms for p distribution (of the daughters)
    fQA->CreateTH1F("h_P_electron", "P distribution of the gamma electrons; p (GeV/c); counts", 100, 0.1, 10);
    fQA->CreateTH1F("h_P_K0pion", "P distribution of the K0 pions; p (GeV/c); counts", 100, 0.1, 10);
    fQA->CreateTH1F("h_P_Lpion", "P distribution of the Lambda pions; p (GeV/c); counts", 100, 0.1, 10);
    fQA->CreateTH1F("h_P_Lproton", "P distribution of the Lambda protons; p (GeV/c); counts", 100, 0.1, 10);

    // QA pt of the V0
    fQA->CreateTH1F("h_Pt_Gamma", "Pt of the gamma conversion; p_{T} (GeV/c); counts", 100, 0, 10);
    fQA->CreateTH1F("h_Pt_K0", "Pt of the K0; p_{T} (GeV/c); counts", 100, 0, 10);
    fQA->CreateTH1F("h_Pt_Lambda", "Pt of the Lambda; p_{T} (GeV/c); counts", 100, 0, 10);    
    
    
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
  if(!fPrimaryVertex) return;
  fV0cuts->SetInputEvent(fInputEvent);
  fV0cuts->SetPrimaryVertex(fPrimaryVertex);
  Int_t v0status = 0;
  for(Int_t iv0 = 0; iv0 < fInputEvent->GetNumberOfV0s(); iv0++){
    if(!TString(fInputEvent->IsA()->GetName()).CompareTo("AliESDEvent")){
      // case ESD
      SetESDanalysis();
      AliESDv0 *esdV0 = (dynamic_cast<AliESDEvent *>(fInputEvent))->GetV0(iv0);
      if(!esdV0) continue;
      if(!esdV0->GetOnFlyStatus()) continue; // Take only V0s from the On-the-fly v0 finder
      v0status = ProcessV0(esdV0);
    } else {
      // case AOD
      SetAODanalysis();
      AliAODv0 *aodV0 = (dynamic_cast<AliAODEvent *>(fInputEvent))->GetV0(iv0);
      if(aodV0->GetOnFlyStatus()) continue; // Take only V0s from the On-the-fly v0 finder
      v0status = ProcessV0(aodV0);
      if(kUndef != v0status){
      }
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
  daughter[1] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack((dynamic_cast<AliESDv0 *>(v0))->GetNindex()));
  if(!daughter[0] || !daughter[1]) return kUndef;

  
  if(IsESDanalysis()){
    for(Int_t i=0; i<2; ++i){
      // check common single track cuts
      if(!fV0cuts->TrackCutsCommon(dynamic_cast<AliESDtrack*>(daughter[i]))) return kUndef;
    }    
    // check commom V0 cuts
    if(!fV0cuts->V0CutsCommon(dynamic_cast<AliESDv0 *>(v0))) return kUndef;
  }


  // store the resutls
  if(IsGammaConv(v0)){
    fQA->Fill("h_nV0s", kRecoGamma);
    return kRecoGamma;
  }
  else if(IsK0s(v0)){
    fQA->Fill("h_nV0s", kRecoK0s);
    return kRecoK0s;
  }
  else if(IsLambda(v0)){
    fQA->Fill("h_nV0s", kRecoLambda);    
    return kRecoLambda;
  }
  else return kUndef;
    
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
Bool_t AliHFEV0pid::IsGammaConv(TObject *v0){
  //
  // Identify Gamma
  //
  AliVTrack* daughter[2];
  Int_t pIndex = 0, nIndex = 0;
  Double_t invMass = 0.;
  Double_t mPt = 0.;
  Int_t v0id = -1;
  if(IsESDanalysis()){
    // ESD - cut V0
    AliESDv0 *esdV0 = dynamic_cast<AliESDv0 *>(v0);
    v0id = esdV0->GetLabel();
    // apply FULL gamma cuts
    if(!fV0cuts->GammaCuts(esdV0)) return kFALSE;
    invMass = esdV0->GetEffMass(AliPID::kElectron, AliPID::kElectron);
    pIndex = esdV0->GetPindex();
    nIndex = esdV0->GetNindex();
    mPt = esdV0->Pt();
  } else {
    // AOD Analysis - not possible to cut
    AliAODv0 *aodV0 = dynamic_cast<AliAODv0 *>(v0);
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
    fIndices->Add(daughter[0]->GetID(), AliHFEV0pid::kRecoElectron);
  }
  if(!fIndices->Find(daughter[1]->GetID())){
    AliDebug(1, Form("Gamma identified, daughter IDs: %d,%d", daughter[1]->GetID(), daughter[1]->GetID()));
    fElectrons->Add(new AliHFEV0info(daughter[1], daughter[0]->GetID(), v0id));
    fIndices->Add(daughter[1]->GetID(), AliHFEV0pid::kRecoElectron);
  }
  fGammas->Add(v0);
  
  return kTRUE;
}
//____________________________________________________________
Bool_t AliHFEV0pid::IsK0s(TObject *v0){
  //
  // Identify K0s
  //
  AliVTrack* daughter[2];
  Int_t pIndex = 0, nIndex = 0;
  Int_t v0id = -1;
  Double_t invMass = 0.;
  Double_t mPt = 0.;
  if(IsESDanalysis()){
    // ESD - cut V0
    AliESDv0 *esdV0 = dynamic_cast<AliESDv0 *>(v0);
    if(!fV0cuts->K0Cuts(esdV0)) return kFALSE;
    v0id = esdV0->GetLabel();
    pIndex = esdV0->GetPindex();
    nIndex = esdV0->GetNindex();
    invMass = esdV0->GetEffMass(AliPID::kPion, AliPID::kPion);
    mPt = esdV0->Pt();
  } else {
    // AOD Analysis - not possible to cut
    AliAODv0 *aodV0 = dynamic_cast<AliAODv0 *>(v0);
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
    fIndices->Add(daughter[0]->GetID(), AliHFEV0pid::kRecoPionK0);
  }
  if(!fIndices->Find(daughter[1]->GetID())){
    AliDebug(1, Form("Adding K0 Pion track with ID %d", daughter[1]->GetID()));
    fPionsK0->Add(new AliHFEV0info(daughter[1], daughter[0]->GetID(), v0id));
    fIndices->Add(daughter[1]->GetID(), AliHFEV0pid::kRecoPionK0);
  }
  fK0s->Add(v0);
  return kTRUE; 
}

//____________________________________________________________
Bool_t AliHFEV0pid::IsPhi(TObject *v0){
  //
  // Identify Phi - very preliminary - requires diffrent approach (V0 fnder is not effective)
  //

  //const Double_t kPhiMass=TDatabasePDG::Instance()->GetParticle(333)->Mass();  // PDG phi mass
  //AliVTrack* daughter[2];
  //Double_t invMass = 0.;
 
  Int_t pIndex = 0, nIndex = 0;
  if(IsESDanalysis()){
    // ESD - cut V0
    AliESDv0 *esdV0 = dynamic_cast<AliESDv0 *>(v0);
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
  AliVTrack* daughter[2];
  Int_t pIndex = 0, nIndex = 0;
  Double_t invMass = 0.;
  Bool_t isLambda = kTRUE; // Lambda - kTRUE, Anti Lambda - kFALSE
  Int_t v0id = -1;
  if(IsESDanalysis()){
    // ESD - cut V0
    AliESDv0 *esdV0 = dynamic_cast<AliESDv0 *>(v0);
    v0id = esdV0->GetLabel();
    if(!fV0cuts->LambdaCuts(esdV0,isLambda)) return kFALSE; 
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
    AliAODv0 *aodV0 = dynamic_cast<AliAODv0 *>(v0);
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
    if(!fIndices->Find(daughter[0]->GetID())){
      fProtons->Add(new AliHFEV0info(daughter[0], daughter[1]->GetID(), v0id));
      fIndices->Add(daughter[0]->GetID(), AliHFEV0pid::kRecoProton);
    }
    if(!fIndices->Find(daughter[1]->GetID())){
      fPionsL->Add(new AliHFEV0info(daughter[1], daughter[0]->GetID(), v0id));
      fIndices->Add(daughter[1]->GetID(), AliHFEV0pid::kRecoPionL);
    }
  }
  // antilambda
  else{
    if(!fIndices->Find(daughter [1]->GetID())){
      fProtons->Add(new AliHFEV0info(daughter[1], daughter[0]->GetID(), v0id));
      fIndices->Add(daughter[1]->GetID(), AliHFEV0pid::kRecoProton);
    }
    if(!fIndices->Find(daughter [0]->GetID())){
      fPionsL->Add(new AliHFEV0info(daughter[0], daughter[1]->GetID(), v0id));
      fIndices->Add(daughter [0]->GetID(), AliHFEV0pid::kRecoPionL);
    }
  }
  if(isLambda) fLambdas->Add(v0);
  else fAntiLambdas->Add(v0);

  return kTRUE;
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
TList *AliHFEV0pid::GetListOfQAhistograms(){
  //
  // Getter for V0 PID QA histograms
  //
  
  TList *tmp = fV0cuts->GetList();
  tmp->SetName("V0cuts");
  fOutput->Add(tmp);
  if(fQA){
    tmp = 0x0;
    tmp = fQA->GetList();
    tmp->SetName("V0pid");
    fOutput->Add(tmp);
  } 
  return fOutput;
}
