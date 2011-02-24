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
// Task for checking the performance of the V0 tender
// 
// 
// Authors:
//   Matus Kalisky <matus.kalisky@cern.ch>
//

#include <TH1F.h>
#include <TList.h>

#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliESDv0.h"
#include "AliESDtrack.h"
#include "AliMCParticle.h"
#include "AliMCEvent.h"
#include "AliESDv0KineCuts.h"
#include "AliKFVertex.h"

#include "AliHFEtools.h"
#include "AliHFEcollection.h"

#include "AliAnalysisTaskCheckV0tenderII.h"

ClassImp(AliAnalysisTaskCheckV0tenderII)

//__________________________________________________________
AliAnalysisTaskCheckV0tenderII::AliAnalysisTaskCheckV0tenderII():
  AliAnalysisTaskSE("CheckV0tenderTask")
  , fOutput(0x0)
  , fColl(0x0)
  , fCollMC(0x0)
  , fV0cuts(0x0)
  , fPrimaryVertex(0x0)
  , fpdgV0(-1)
  , fpdgP(-1)
  , fpdgN(-1)
  , fEvents(0x0)
{
  //
  // Default Constructor
  //


  DefineOutput(1, TH1F::Class());
  DefineOutput(2, TList::Class());


}
//__________________________________________________________
AliAnalysisTaskCheckV0tenderII::AliAnalysisTaskCheckV0tenderII(const Char_t *name):
  AliAnalysisTaskSE(name)
  , fOutput(0x0)
  , fColl(0x0)
  , fCollMC(0x0)
  , fV0cuts(0x0)
  , fPrimaryVertex(0x0)
  , fpdgV0(0)
  , fpdgP(0)
  , fpdgN(0)
  , fEvents(0x0)
{
  //
  // Constructor
  //

  DefineOutput(1, TH1F::Class());
  DefineOutput(2, TList::Class());

}
//__________________________________________________________
AliAnalysisTaskCheckV0tenderII::~AliAnalysisTaskCheckV0tenderII(){
  //
  // Destructor
  //

  if (fOutput) delete fOutput;
  //if (fColl) delete fColl;
  //if (fCollMC) delete fCollMC;
  if (fV0cuts) delete fV0cuts;
  if (fPrimaryVertex) delete fPrimaryVertex;
  if (fEvents) delete fEvents;
}
//__________________________________________________________
void AliAnalysisTaskCheckV0tenderII::UserCreateOutputObjects(){
  //
  // prepare output objects
  //

  fOutput = new TList();
  //fOutput->SetOwner();
  // Counter for number of events
  fEvents = new TH1I("nEvents", "NumberOfEvents", 1, 1, 2);

  fColl = new AliHFEcollection("V0_QA", "tender V0s for data");
  fCollMC = new AliHFEcollection("V0_MC_QA", "tender V0s for MC");

  fV0cuts = new AliESDv0KineCuts();
  if(!fV0cuts){
    AliError("Failed to initialize AliESDv0KineCuts instance");
    return;
  }

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
  // background study - MC only
  fCollMC->CreateTH1Fvector1(6, "h_Gamma_Bg", "MC gamma bg. - different sources; p_{T} (GeV/c), counts", 20, 0.1, 20, 0);
  fCollMC->CreateTH1Fvector1(6, "h_K0_Bg", "MC K0 bg. - different sources; p_{T} (GeV/c), counts", 20, 0.1, 20, 0);
  fCollMC->CreateTH1Fvector1(6, "h_Lambda_Bg", "MC gamma bg. - different sources; p_{T} (GeV/c), counts", 20, 0.1, 20, 0);
  fCollMC->CreateTH1Fvector1(6, "h_ALambda_Bg", "MC gamma bg. - different sources; p_{T} (GeV/c), counts", 20, 0.1, 20, 0);
  

  TList *tmp = fColl->GetList();
  tmp->SetName(fColl->GetName());
  fOutput->Add(tmp);
  tmp = 0x0;
  tmp = fCollMC->GetList();
  tmp->SetName(fCollMC->GetName());
  fOutput->Add(tmp);
  
  
}
//__________________________________________________________
void AliAnalysisTaskCheckV0tenderII::UserExec(Option_t *){
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

  AliKFVertex *primaryVertex = new AliKFVertex(*(fInputEvent->GetPrimaryVertex()));
  if(!primaryVertex) return;
  
  fV0cuts->SetPrimaryVertex(primaryVertex);
  fV0cuts->SetEvent(fInputEvent);
  
  ProcessV0s();
  

  if (primaryVertex) delete primaryVertex;
  fEvents->Fill(1.1);
  PostData(1, fEvents);
  PostData(2, fOutput);
}
//__________________________________________________________
void AliAnalysisTaskCheckV0tenderII::Terminate(Option_t *){
  //
  // Do Post Processing
  //
}
//__________________________________________________________
void AliAnalysisTaskCheckV0tenderII::ProcessV0s(){
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

    ResetPDGcodes();

    // New standalone V0 selection software
    if( ! (fV0cuts->ProcessV0(esdV0, fpdgV0, fpdgP, fpdgN)) ) continue;
    
    Int_t pid = PDGtoPIDv0(fpdgV0);
    //printf(" -D: pdg: %i, pid: %i\n", fpdgV0, pid);
    if(pid <= -1) continue;
    
    fColl->Fill("h_NumberOf_V0s", pid);
    name = "h_" + type[pid] + "_pt";
    Float_t pT = esdV0->Pt();
    fColl->Fill(name, pT);
    Float_t mass = MassV0(esdV0, pid);
    name = "h_" + type[pid] + "_mass";
    fColl->Fill(name, mass);

    ProcessDaughters(esdV0);
    if(fMCEvent){
      ProcessV0sMC(esdV0);
      ProcessDaughtersMC(esdV0);
      ProcessBackground(esdV0);
    }

  }

}
//__________________________________________________________
void AliAnalysisTaskCheckV0tenderII::ProcessDaughters(AliESDv0 * const v0){
  //
  // produce some check plots for V0 tender tagged single tracks
  //

  const TString type[3] = {"Electron", "Pion", "Proton"};
  TString name;

  if(!v0) return;

  // daughter tracks
  const Int_t nTracks = fInputEvent->GetNumberOfTracks();
  AliESDtrack *d[2];
  Int_t iN, iP;
  iN = iP = -1;
  iP = v0->GetPindex();
  iN = v0->GetNindex();
  if(iN < 0 || iP < 0) return;
  if(iN >= nTracks || iP >= nTracks) return;
  d[0] = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(iP));
  d[1] = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(iN));  
  if(!d[0] || !d[1]) return;

  for(Int_t i=0; i<2; ++i){
    Int_t  pid = (i == 0) ? PDGtoPID(fpdgP) : PDGtoPID(fpdgN);
    if(pid < 0) continue;
    fColl->Fill("h_NumberOfDaughters", pid*1.0);
    Float_t pT = d[i]->Pt();
    name = "h_" + type[pid] + "_pt";
    fColl->Fill(name, pT);
  }
}
//__________________________________________________________
void AliAnalysisTaskCheckV0tenderII::ProcessV0sMC(AliESDv0 * const v0){
  //
  // check all V0tender selected V0 on their true identity
  // 

  const TString type[4] = {"Gamma", "K0", "Lambda", "ALambda"};
  TString name;
  
  Int_t pid = PDGtoPIDv0(fpdgV0);
  if(pid < 0) return;

  // true if fpdgV0 agrees with MC pdg of the mother
  Bool_t id = kFALSE;
  const Int_t nTracks = fInputEvent->GetNumberOfTracks();
  // both ESD daughtr tracks
  Int_t iN, iP;
  iN = iP = -1;
  iP = v0->GetPindex();
  iN = v0->GetNindex();    
  if(iN < 0 || iP < 0) return;
  if(iN >= nTracks || iP >= nTracks) return;    
  AliESDtrack *dP, *dN;
  dP = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(iP));
  dN = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(iN));  
  if(!dN || !dP) return;

  // MC labels of the daughter tracks
  Int_t lN, lP;
  lN = dN->GetLabel();
  lP = dP->GetLabel();
  if(lN < 0 || lP < 0) return;

  // MC daughter particles
  AliMCParticle *mcP, *mcN;
  mcP = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(lP));
  mcN = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(lN));
  if(!mcP || !mcN) return;

  // labels of the mother particles
  Int_t lPm, lNm;
  lPm = mcP->GetMother();
  lNm = mcN->GetMother();
  AliMCParticle *m = 0x0;
  // particles without mother particle are most probably
  // primary particles
  if(lPm >= 0){
    m = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(lPm)); 
    if(!m) return;
  }
  if(m){
    Int_t pdg = m->PdgCode();
    //if(lPm == lNm) printf(" -D: pdg: %i, fpdgV0: %i \n", pdg, fpdgV0);
    if((lPm == lNm) && (pdg == fpdgV0)) id = kTRUE;
  }
  
  if(id) name = "h_" + type[pid] + "_pt_S";
  else  name = "h_" + type[pid] + "_pt_B";
  Float_t pT = v0->Pt();
  fCollMC->Fill(name, pT);

  if(id) name = "h_" + type[pid] + "_mass_S";
  else name = "h_" + type[pid] + "_mass_B";
  Float_t mass = MassV0(v0, pid);
  fCollMC->Fill(name, mass);

  
}
//__________________________________________________________
void AliAnalysisTaskCheckV0tenderII::ProcessDaughtersMC(AliESDv0 * const v0){
  //
  // check the identity of the V0tender selected V0 daughters
  // !!! for positive check only the true identity plays a role here, 
  // not the true V0 mother identity (e.g. selected electron could come
  // from primary vertex or pion dalitz deca or true gamma conversion) !!!!
  //

  const TString type[3] = {"Electron", "Pion", "Proton"};
  TString name;

  if(!v0) return;

  // daughter tracks
  const Int_t nTracks = fInputEvent->GetNumberOfTracks();
  AliESDtrack *d[2];
  Int_t iN, iP;
  iN = iP = -1;
  iP = v0->GetPindex();
  iN = v0->GetNindex();
  if(iN < 0 || iP < 0) return;
  if(iN >= nTracks || iP >= nTracks) return;
  d[0] = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(iP));
  d[1] = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(iN));  
  if(!d[0] || !d[1]) return;

  //printf(" *** fpdgV0: %i, fpdgP: %i, fpdgN: %i \n", fpdgV0, fpdgP, fpdgN);

  for(Int_t i=0; i<2; ++i){
    Bool_t id = kFALSE;
    Int_t  pid = (i == 0) ? PDGtoPID(fpdgP) : PDGtoPID(fpdgN);
    Int_t  pdg  = (i == 0) ? fpdgP : fpdgN;
    if(pid < 0) continue;
    Float_t pT = d[i]->Pt();
    Int_t label = d[i]->GetLabel();
    if(label < 0) continue;
    AliMCParticle *mcp = 0x0;
    mcp = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(label));
    if(!mcp) continue;
    Int_t pdgMC = mcp->PdgCode();
    //printf("   * pid: %i, pdg: %i, pdgMC: %i \n", pid, pdg, pdgMC);
    if(pdgMC == pdg) id = kTRUE;
    if(id) name =  "h_" + type[pid] + "_pt_S";
    else name =  "h_" + type[pid] + "_pt_B";
    fCollMC->Fill(name, pT);
  }
}
//__________________________________________________________
void  AliAnalysisTaskCheckV0tenderII::ProcessBackground(AliESDv0 * const v0){
  //
  // look at the compostition and the properties of the
  // miss-identified V0
  //

  if(!v0) return;
  Float_t pt = v0->Pt();
  const TString type[4] = {"Gamma", "K0", "Lambda", "ALambda"};
  TString name;
  // daughter tracks
  const Int_t nTracks = fInputEvent->GetNumberOfTracks();
  AliESDtrack *d[2];
  Int_t index[2] = {-1, -1};
  index[0] = v0->GetPindex();
  index[1] = v0->GetNindex();
  if(index[0] < 0 || index[1] < 0) return;
  if(index[0] >= nTracks || index[1] >= nTracks) return;
  d[0] = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(index[0]));
  d[1] = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(index[1]));  
  if(!d[0] || !d[1]) return;

  // daughter MC particles
  Int_t label[2] = {d[0]->GetLabel(), d[1]->GetLabel()};
  if(label[0] < 0 || label[1] < 0) return;
  AliMCParticle *dmc[2] = {0x0, 0x0};
  dmc[0] = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(label[0]));
  dmc[1] = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(label[1]));
  if(!dmc[0]  || !dmc[1]) return;

  Bool_t primary[2] = {dmc[0]->Particle()->IsPrimary(), dmc[1]->Particle()->IsPrimary()};

  // mother MC particles
  Int_t labelM[2] = {dmc[0]->GetMother(), dmc[1]->GetMother()};
  AliMCParticle *mmc[2] = {0x0, 0x0};
  if(!primary[0])  mmc[0] = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(labelM[0]));
  if(!primary[1])  mmc[1] = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(labelM[1]));

  // if particle is not primary it MUST have a mother particle
  if(!primary[0] && !mmc[0]) return;
  if(!primary[1] && !mmc[1]) return;

  Int_t pdgM[2] = {0, 0};
  if(mmc[0]) pdgM[0] = mmc[0]->PdgCode();
  if(mmc[1]) pdgM[1] = mmc[1]->PdgCode();

  // select only miss-identified V0s
  if(labelM[0] == labelM[1] && pdgM[0] == fpdgV0) return;

  // indices for V0 background histograms
  // 0 - random match
  // 1 - gamma
  // 2 - K0
  // 3 - lambda 
  // 4 - a-lambda
  // 5 - other dacay or interaction

  Int_t ixM = -1;  

  Int_t typeV0  = (labelM[0] != labelM[1]) ? 0 : 1;
  if(0 == typeV0) ixM = 0;
  else{
    ixM = PDGtoPIDv0(pdgM[0]) + 1;  
    if(ixM < 0) ixM = 5; 
  }
  Int_t ix = PDGtoPIDv0(fpdgV0);
  if(0 <= ix ){
    name = "h_" + type[PDGtoPIDv0(fpdgV0)] + "_Bg"; 
    fCollMC->Fill(name, ixM, pt);
  }
  
  // now look at the daughter tracks
  

}
//__________________________________________________________
Float_t AliAnalysisTaskCheckV0tenderII::MassV0(AliESDv0 * const v0, Int_t id){
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
Bool_t AliAnalysisTaskCheckV0tenderII::CheckSigns(AliESDv0 * const v0){
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
//__________________________________________________________
Int_t  AliAnalysisTaskCheckV0tenderII::PDGtoPIDv0(Int_t pdgV0) const {
  //
  // convert thereconstructed V0 pdg to local pid
  //

  switch(pdgV0){
  case 22: return 0; break;
  case 310: return 1; break;
  case 3122: return 2; break;
  case -3122: return 3; break;
  }
  
  return -1;

}
//__________________________________________________________
Int_t  AliAnalysisTaskCheckV0tenderII::PDGtoPID(Int_t pdg) const {
  //
  // convert daughter pdg code to local pid
  //
  switch(TMath::Abs(pdg)){
  case 11: return 0; break;
  case 211: return 1; break;
  case 2212: return 2; break;
  }
  return -1;

}
//__________________________________________________________
void  AliAnalysisTaskCheckV0tenderII::ResetPDGcodes(){
  //
  // reset the PDG codes of the V0 and daughter
  // particles to default values
  //

  fpdgV0 = -1;
  fpdgP = -1;
  fpdgN = -1;

}
