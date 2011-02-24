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
// Task fir checking the performance of the V0 tender
// 
// 
// Authors
//   Matus Kalisky <matus.kalisky@cern.ch>
//

#include <sstream>

#include <TH1F.h>
#include <TList.h>

#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliESDv0.h"
#include "AliESDtrack.h"
#include "AliMCParticle.h"
#include "AliMCEvent.h"

#include "AliHFEtools.h"
#include "AliHFEcollection.h"

#include "AliAnalysisTaskCheckV0tender.h"

ClassImp(AliAnalysisTaskCheckV0tender)

//__________________________________________________________
AliAnalysisTaskCheckV0tender::AliAnalysisTaskCheckV0tender():
  AliAnalysisTaskSE("CheckV0tenderTask")
  , fOutput(0x0)
  , fColl(0x0)
  , fCollMC(0x0)
  , fEvents(0x0)
{
  //
  // Default Constructor
  //
}
//__________________________________________________________
AliAnalysisTaskCheckV0tender::AliAnalysisTaskCheckV0tender(const Char_t *name):
  AliAnalysisTaskSE(name)
  , fOutput(0x0)
  , fColl(0x0)
  , fCollMC(0x0)
  , fEvents(0x0)
{
  //
  // Default Constructor
  //
  DefineOutput(1, TH1F::Class());
  DefineOutput(2, TList::Class());

}
//__________________________________________________________
AliAnalysisTaskCheckV0tender::~AliAnalysisTaskCheckV0tender(){
  //
  // Destructor
  //

  if (fOutput) delete fOutput;
  //if (fColl) delete fColl;
  //if (fCollMC) delete fCollMC;
  if (fEvents) delete fEvents;
}
//__________________________________________________________
void AliAnalysisTaskCheckV0tender::UserCreateOutputObjects(){
  //
  // prepare output objects
  //

  fOutput = new TList();
  fOutput->SetOwner();
  // Counter for number of events
  fEvents = new TH1I("nEvents", "NumberOfEvents", 1, 1, 2);

  fColl = new AliHFEcollection("V0_QA", "tender V0s for data");
  fCollMC = new AliHFEcollection("V0_MC_QA", "tender V0s for MC");

  // 
  // Data histos
  //
  
  // total number of tagged V0s
  fColl->CreateTH1F("h_NumberOf_V0s", "Number of tagged V0s; type; counts", 4, -0.5, 3.5);
  // pT spectra of the tagged V0s
  fColl->CreateTH1F("h_Gamma_pt", "p_{T} spectrum of tagged gammas; p_{T} (GeV/c); counts", 20, 0.1, 20, 0);
  fColl->CreateTH1F("h_K0_pt", "p_{T} spectrum of tagged K0s; p_{T} (GeV/c); counts", 20, 0.1, 20, 0);
  fColl->CreateTH1F("h_Lambda_pt", "p_{T} spectrum of tagged Lambdas; p_{T} (GeV/c); counts", 20, 0.1, 20, 0);
  fColl->CreateTH1F("h_ALambda_pt", "p_{T} spectrum of tagged A-Lambdas; p_{T} (GeV/c); counts", 20, 0.1, 20, 0);
  // invariant mass of the V0s
  fColl->CreateTH1F("h_Gamma_mass", "Inv. Mass of the ", 100, 0., 0.1);
  fColl->CreateTH1F("h_K0_mass", "Inv. Mass of the ", 100, 0.45, 0.55);
  fColl->CreateTH1F("h_Lambda_mass", "Inv. Mass of the ", 100, 1.08, 1.14);
  fColl->CreateTH1F("h_ALambda_mass", "Inv. Mass of the ", 100, 1.08, 1.14);

  // total number of tagged daughter particles (should correlate with number of V0s !
  fColl->CreateTH1F("h_NumberOfDaughters", "Number of tagged daughters; type; counts", 3, -0.5, 2.5);
  // pT spectra of tagged daughter particles
  fColl->CreateTH1F("h_Electron_pt", "p_{T} spectrum of tagged; p_{T} (GeV/c); counts", 20, 0.1, 20, 0);
  fColl->CreateTH1F("h_Pion_pt", "p_{T} spectrum of tagged; p_{T} (GeV/c); counts", 20, 0.1, 20, 0);
  fColl->CreateTH1F("h_Proton_pt", "p_{T} spectrum of tagged; p_{T} (GeV/c); counts", 20, 0.1, 20, 0);

  //
  // MC histos
  //

  // pT spectra of the tagged V0s
  fCollMC->CreateTH1F("h_Gamma_pt_S", "MC-S: p_{T} spectrum of tagged gammas; p_{T} (GeV/c); counts", 20, 0.1, 20, 0);
  fCollMC->CreateTH1F("h_K0_pt_S", "MC-S: p_{T} spectrum of tagged K0s; p_{T} (GeV/c); counts", 20, 0.1, 20, 0);
  fCollMC->CreateTH1F("h_Lambda_pt_S", "MC-S: p_{T} spectrum of tagged Lambdas; p_{T} (GeV/c); counts", 20, 0.1, 20, 0);
  fCollMC->CreateTH1F("h_ALambda_pt_S", "MC-S: p_{T} spectrum of tagged A-Lambdas; p_{T} (GeV/c); counts", 20, 0.1, 20, 0); 
  fCollMC->CreateTH1F("h_Gamma_pt_B", "MC-B: p_{T} spectrum of tagged gammas; p_{T} (GeV/c); counts", 20, 0.1, 20, 0);
  fCollMC->CreateTH1F("h_K0_pt_B", "MC-B: p_{T} spectrum of tagged K0s; p_{T} (GeV/c); counts", 20, 0.1, 20, 0);
  fCollMC->CreateTH1F("h_Lambda_pt_B", "MC-B: p_{T} spectrum of tagged Lambdas; p_{T} (GeV/c); counts", 20, 0.1, 20, 0);
  fCollMC->CreateTH1F("h_ALambda_pt_B", "MC-B: p_{T} spectrum of tagged A-Lambdas; p_{T} (GeV/c); counts", 20, 0.1, 20, 0); 
  // invariant mass of the V0s
  fCollMC->CreateTH1F("h_Gamma_mass_S", "MC-S: Inv. Mass of the gamma; m (GeV/c^{2}); counts", 100, 0., 0.1);
  fCollMC->CreateTH1F("h_K0_mass_S", "MC-S: Inv. Mass of the K0; m (GeV/c^{2}); counts", 100, 0.45, 0.55);
  fCollMC->CreateTH1F("h_Lambda_mass_S", "MC-S: Inv. Mass of the Lambda; m (GeV/c^{2}); counts", 100, 1.08, 1.14);
  fCollMC->CreateTH1F("h_ALambda_mass_S", "MC-S: Inv. Mass of the A-Lambda; m (GeV/c^{2}); counts", 100, 1.08, 1.14);
  fCollMC->CreateTH1F("h_Gamma_mass_B", "MC-B: Inv. Mass of the gamma; m (GeV/c^{2}); counts", 100, 0., 0.1);
  fCollMC->CreateTH1F("h_K0_mass_B", "MC-B: Inv. Mass of the K0; m (GeV/c^{2}); counts", 100, 0.45, 0.55);
  fCollMC->CreateTH1F("h_Lambda_mass_B", "MC-B: Inv. Mass of the Lambda; m (GeV/c^{2}); counts", 100, 1.08, 1.14);
  fCollMC->CreateTH1F("h_ALambda_mass_B", "MC-B: Inv. Mass of the A-Lambda; m (GeV/c^{2}); counts", 100, 1.08, 1.14);
  // pT spectra of tagged daughter particles
  fCollMC->CreateTH1F("h_Electron_pt_S", "MC-S: p_{T} spectrum of tagged electrons; p_{T} (GeV/c); counts", 20, 0.1, 20, 0);
  fCollMC->CreateTH1F("h_Pion_pt_S", "MC-S: p_{T} spectrum of tagged pions; p_{T} (GeV/c); counts", 20, 0.1, 20, 0);
  fCollMC->CreateTH1F("h_Proton_pt_S", "MC-S: p_{T} spectrum of tagged protons; p_{T} (GeV/c); counts", 20, 0.1, 20, 0);
  fCollMC->CreateTH1F("h_Electron_pt_B", "MC-B: p_{T} spectrum of tagged electrons; p_{T} (GeV/c); counts", 20, 0.1, 20, 0);
  fCollMC->CreateTH1F("h_Pion_pt_B", "MC-B: p_{T} spectrum of tagged pions; p_{T} (GeV/c); counts", 20, 0.1, 20, 0);
  fCollMC->CreateTH1F("h_Proton_pt_B", "MC-B: p_{T} spectrum of tagged protons; p_{T} (GeV/c); counts", 20, 0.1, 20, 0);


  TList *tmp = fColl->GetList();
  tmp->SetName(fColl->GetName());
  fOutput->Add(tmp);
  tmp = 0x0;
  tmp = fCollMC->GetList();
  tmp->SetName(fCollMC->GetName());
  fOutput->Add(tmp);
  
  
}
//__________________________________________________________
void AliAnalysisTaskCheckV0tender::UserExec(Option_t *){
  //
  // Event Loop
  // 

  AliMCEventHandler* mcHandler = (dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()));
  //AliESDInputHandler *inh = dynamic_cast<AliESDInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  //AliESDpid *workingPID = NULL;
  //if(inh && (workingPID = inh->GetESDpid())) workingPID = inh->GetESDpid();
  //else workingPID = AliHFEtools::GetDefaultPID(mcHandler ? kTRUE : kFALSE);
     
  // check the MC data
  if(fMCEvent && !mcHandler ) return;
  if(fMCEvent &&  !mcHandler->InitOk() ) return;
  if(fMCEvent &&  !mcHandler->TreeK() ) return;
  if(fMCEvent &&  !mcHandler->TreeTR() ) return;

  ProcessV0s();
  ProcessDaughters();

  if(fMCEvent){
    ProcessV0sMC();
    ProcessDaughtersMC();
  }

  fEvents->Fill(1.1);
  PostData(1, fEvents);
  PostData(2, fOutput);
}
//__________________________________________________________
void AliAnalysisTaskCheckV0tender::Terminate(Option_t *){
  //
  // Do Post Processing
  //
}
//__________________________________________________________
void AliAnalysisTaskCheckV0tender::ProcessV0s(){
  //
  // loop over the V0s and extract the information about
  // the V0s tagged by the V0 tender
  //

  const TString type[4] = {"Gamma", "K0", "Lambda", "ALambda"};
  TString name;

  Int_t nV0s = fInputEvent->GetNumberOfV0s();
  for(Int_t i=0; i<nV0s; ++i){
    AliESDv0 *esdV0 = (static_cast<AliESDEvent *>(fInputEvent))->GetV0(i);
    if(!esdV0) continue;
    if(!esdV0->GetOnFlyStatus()) continue; // Take only V0s from the On-the-fly v0 finder
    Int_t pid = GetTenderPidV0(esdV0);
    if(pid < 0) continue;
    fColl->Fill("h_NumberOf_V0s", pid);
    Float_t pT = esdV0->Pt();
    name = "h_" + type[pid] + "_pt";
    fColl->Fill(name, pT);
    Float_t mass = MassV0(esdV0, pid);    
    name = "h_" + type[pid] + "_mass";
    //printf(" -D: name: %s \n", name.Data());
    fColl->Fill(name, mass);
  }

}
//__________________________________________________________
void AliAnalysisTaskCheckV0tender::ProcessDaughters(){
  //
  // produce some check plots for V0 tender tagged single tracks
  //

  const TString type[3] = {"Electron", "Pion", "Proton"};
  TString name;

  Int_t nTracks = fInputEvent->GetNumberOfTracks();
  for(Int_t i=0; i<nTracks; ++i){
    AliESDtrack *track = dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(i));
    if(!track) continue;
    Int_t  pid = GetTenderPidDaughter(track);
    if(pid < 0) continue;
    fColl->Fill("h_NumberOfDaughters", pid*1.0);
    Float_t pT = track->Pt();
    name = "h_" + type[pid] + "_pt";
    fColl->Fill(name, pT);
    
  }
}
//__________________________________________________________
void AliAnalysisTaskCheckV0tender::ProcessV0sMC(){
  //
  // check all V0tender selected V0 on their true identity
  // 

  const Int_t pid2pdg[4] = {22, 310, 3122, -3122};
  const TString type[4] = {"Gamma", "K0", "Lambda", "ALambda"};
  TString name;
  Int_t nV0s = fInputEvent->GetNumberOfV0s();

  Int_t nTracks = fInputEvent->GetNumberOfTracks();

  // V0 loop
  for(Int_t i=0; i<nV0s; ++i){
    Bool_t id = kFALSE;
    AliESDv0 *esdV0 = (static_cast<AliESDEvent *>(fInputEvent))->GetV0(i);
    if(!esdV0) continue;
    if(!esdV0->GetOnFlyStatus()) continue; // Take only V0s from the On-the-fly v0 finder
    Int_t pid = GetTenderPidV0(esdV0);
    if(pid < 0) continue;

    // both ESD daughtr tracks
    Int_t iN, iP;
    iN = iP = -1;
    iP = esdV0->GetPindex();
    iN = esdV0->GetNindex();    
    if(iN < 0 || iP < 0) continue;
    if(iN >= nTracks || iP >= nTracks) continue;    
    AliESDtrack *dP, *dN;
    dP = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(iP));
    dN = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(iN));  
    if(!dN || !dP) continue;

    // MC labels of the daughter tracks
    Int_t lN, lP;
    lN = dN->GetLabel();
    lP = dP->GetLabel();
    if(lN < 0 || lP < 0) continue;

    // MC daughter particles
    AliMCParticle *mcP, *mcN;
    mcP = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(lP));
    mcN = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(lN));
    if(!mcP || !mcN) continue;

    // labels of the mother particles
    Int_t lPm, lNm;
    lPm = mcP->GetMother();
    lNm = mcN->GetMother();
    if(lPm < 0) continue;
    AliMCParticle *m = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(lPm)); 
    if(!m) continue;
    Int_t pdg = m->PdgCode();
    if((lPm == lNm) && (pdg == pid2pdg[pid])) id = kTRUE;

    if(id) name = "h_" + type[pid] + "_pt_S"; 
    else name = "h_" + type[pid] + "_pt_B"; 
    Float_t pT = esdV0->Pt();
    fCollMC->Fill(name, pT);

    if(id) name = "h_" + type[pid] + "_mass_S";
    else name = "h_" + type[pid] + "_mass_B";
    Float_t mass = MassV0(esdV0, pid);
    fCollMC->Fill(name, mass);

   }
}
//__________________________________________________________
void AliAnalysisTaskCheckV0tender::ProcessDaughtersMC(){
  //
  // check the identity of the V0tender selected V0 daughters
  // !!! for positive check only the true identity plays a role here, 
  // not the true V0 mother identity (e.g. selected electron could come
  // from primary vertex or pion dalitz deca or true gamma conversion) !!!
  //

  const Int_t pid2pdg [3] = {11, 211, 2212};
  const TString type[3] = {"Electron", "Pion", "Proton"};
  TString name;


  Int_t nTracks = fInputEvent->GetNumberOfTracks();
  for(Int_t i=0; i<nTracks; ++i){
    Bool_t id = kFALSE;
    AliESDtrack *track = dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(i));
    if(!track) continue;
    Int_t  pid = GetTenderPidDaughter(track);
    if(pid < 0) continue;
    Float_t pT = track->Pt();
    Int_t label = track->GetLabel();
    if(label < 0) continue;
    AliMCParticle *mcp = 0x0;
    mcp = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(label));
    if(!mcp) continue;
    Int_t pdg = TMath::Abs(mcp->PdgCode());
    if(pdg == pid2pdg[pid]) id = kTRUE;
    if(id) name = "h_" + type[pid] + "_pt_S";
    else name = "h_" + type[pid] + "_pt_B";
    fCollMC->Fill(name, pT);
  }
}
//__________________________________________________________
Int_t AliAnalysisTaskCheckV0tender::GetTenderPidV0(AliESDv0 * const v0){
  //
  // retrieve the PID nformation stored in the status flags by the train tender
  // 
  Int_t pid = -1;
  if(!v0){
    return pid;
  }
  Int_t nTimes = 0;
  if(v0->TestBit(BIT(14))){
    pid = 0;
    nTimes++;
  }
  if(v0->TestBit(BIT(15))){
    pid = 1;
    nTimes++;
  }
  if(v0->TestBit(BIT(16))){
    pid = 2;
    nTimes++;
  }
  if(v0->TestBit(BIT(17))){
    pid = 3;
    nTimes++;
  }
  if(nTimes > 1){
    AliWarning("V0 track labeled multiple times by the V0 tender");
    pid = -1;
  }    

  //printf(" -D: pid: %i \n", pid);

  return pid;
}
//__________________________________________________________
Int_t AliAnalysisTaskCheckV0tender::GetTenderPidDaughter(AliESDtrack * const track){
  //
  // retrieve the PID nformation stored in the status flags by the train tender
  // 

  Int_t pid = -1;
  if(!track){
    return pid;
  }
  Int_t nTimes = 0;
  if(track->TestBit(BIT(14))){
    pid = 0;
    nTimes++;
  }
  if(track->TestBit(BIT(15))){
    pid = 1;
    nTimes++;
  }
  if(track->TestBit(BIT(16))){
    pid = 2;
    nTimes++;
  }
  if(nTimes > 1){
    AliWarning("V0 track labeled multiple times by the V0 tender");
    pid = -1;
  }    
  return pid;
}
//__________________________________________________________
Float_t AliAnalysisTaskCheckV0tender::MassV0(AliESDv0 * const v0, Int_t id){
  //
  // Get the V0 effective mass
  //

  Float_t mass = -0.1;
  Bool_t sign = CheckSigns(v0);
  if(0 == id){
    mass = v0->GetEffMass(0, 0);
  }
  else if(1 == id){
    mass = v0->GetEffMass(2, 2);
  }
  else if(2 == id){
    mass = (sign) ? v0->GetEffMass(4, 2) : v0->GetEffMass(2, 4);
  }
  else if(3 == id){
    mass = (sign) ? v0->GetEffMass(2, 4) : v0->GetEffMass(4, 2);
  }
  else{
    AliWarning(Form("Unrecognized V0 id: %i", id));
  }

  return mass;

}
//__________________________________________________________
Bool_t AliAnalysisTaskCheckV0tender::CheckSigns(AliESDv0 * const v0){
  //
  // check wheter the sign was correctly applied to 
  // V0 daughter tracks
  // This function should become obsolete once the V0 finder will be updated (fixed)
  //
  
  Bool_t correct = kFALSE;

  Int_t pIndex = 0, nIndex = 0;
  pIndex = v0->GetPindex();
  nIndex = v0->GetNindex();
  
  AliESDtrack* d[2];
  d[0] = dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(pIndex));
  d[1] = dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(nIndex));

  Int_t sign[2];
  sign[0] = (int)d[0]->GetSign();
  sign[1] = (int)d[1]->GetSign();
  
  if(-1 == sign[0] && 1 == sign[1]){
    correct = kFALSE;
    //v0->SetIndex(0, pIndex);  // set the index of the negative v0 track
    //v0->SetIndex(1, nIndex);  // set the index of the positive v0 track
  }
  else{
    correct = kTRUE;
  }
  
  //pIndex = v0->GetPindex();
  //nIndex = v0->GetNindex();
  //printf("-D2: P: %i, N: %i\n", pIndex, nIndex);

  return correct;
}
