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

#include <iostream>
#include <TList.h>

#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliStack.h"


#include "AliFilteredTreeEventCuts.h"

using namespace std;

ClassImp(AliFilteredTreeEventCuts)

Int_t AliFilteredTreeEventCuts::fgLastProcessType = -1;

//_____________________________________________________________________________
AliFilteredTreeEventCuts::AliFilteredTreeEventCuts(const Char_t* name,const Char_t *title) : 
AliAnalysisCuts(name, title)
, fTriggerRequired(kTRUE)
, fRecVertexRequired(kTRUE)
, fEventProcessType(kInvalidProcess)
, fMinNContributors(0)
, fMaxNContributors(0)
, fMaxR(0)
, fMinZv(0)
, fMaxZv(0)
, fMeanXv(0)
, fMeanYv(0)
, fMeanZv(0)
, fSigmaMeanXv(0)
, fSigmaMeanYv(0)
, fSigmaMeanZv(0)
, fRedoTPCVertex(kTRUE)
, fUseBeamSpotConstraint(kTRUE)
, fEventSelectedRequired(kFALSE)
{
  // default constructor 
  
  // init data members with defaults
  Init();
}

//_____________________________________________________________________________
AliFilteredTreeEventCuts::~AliFilteredTreeEventCuts()  
{
  // destructor
}

//_____________________________________________________________________________
void AliFilteredTreeEventCuts::Init()  
{
  // set default values
  SetTriggerRequired();
  SetRecVertexRequired();
  SetEventProcessType();
  SetNContributorsRange();
  SetMaxR();
  SetZvRange();
  SetMeanXYZv();
  SetSigmaMeanXYZv();
  SetRedoTPCVertex();
  SetUseBeamSpotConstraint();
}

//_____________________________________________________________________________
Bool_t AliFilteredTreeEventCuts::AcceptEvent(AliESDEvent *esdEvent,AliMCEvent *mcEvent, const AliESDVertex *vtx)
{
  // Check event selection cuts
  Bool_t retValue=kTRUE;

  if(!esdEvent) return kFALSE;
  if(!IsRecVertexRequired()) return kTRUE;
  if(!vtx) return kFALSE;
  if(!vtx->GetStatus()) return kFALSE;

  if(mcEvent) {
   // check MC event conditions
   AliHeader* header = mcEvent->Header();
   if(!header) return kFALSE;
  
    // select event type (ND-non diffractive, SD-single diffractive, DD-double diffractive)
    if(fEventProcessType == kInvalidProcess) { 
      retValue=kTRUE;
    } 
    else if(fEventProcessType == kSD || fEventProcessType == kDD) {
      MCProcessType processType = GetEventProcessType(header);
      if(processType == kND) retValue=kFALSE;
      else retValue=kTRUE;
    }
    else if(fEventProcessType == GetEventProcessType(header)) { 
      retValue=kTRUE;
    }
    else 
      retValue=kFALSE;
  }

  if(vtx->GetZ() < fMinZv) return kFALSE; 
  if(vtx->GetZ() > fMaxZv) return kFALSE; 

return retValue;  
}

//_____________________________________________________________________________
Bool_t AliFilteredTreeEventCuts::AcceptMCEvent(AliMCEvent *mcEvent)
{
  // Check event selection cuts
  if(!mcEvent) return kFALSE;

  Bool_t retValue=kTRUE;

  // check MC event conditions
  AliHeader* header = mcEvent->Header();
  if(!header) return kFALSE;

  AliGenEventHeader* genHeader = header->GenEventHeader();
  if (!genHeader) {
    AliDebug(AliLog::kError, "Could not retrieve genHeader from Header");
    return kFALSE;
  }
  TArrayF vtxMC(3);
  genHeader->PrimaryVertex(vtxMC);
  
  // select event type (ND-non diffractive, SD-single diffractive, DD-double diffractive)
  if(fEventProcessType == kInvalidProcess) { 
     retValue=kTRUE;
  } else {
     if(fEventProcessType == GetEventProcessType(header)) retValue=kTRUE;
     else retValue=kFALSE;
  }

  /*
  Float_t R = TMath::Sqrt(vtxMC[0]*vtxMC[0]+vtxMC[1]*vtxMC[1]);
  if(R > fMaxR) return kFALSE; 
  */

  if(vtxMC[2] < fMinZv) return kFALSE; 
  if(vtxMC[2] > fMaxZv) return kFALSE; 

return retValue;  
}

//_____________________________________________________________________________
Long64_t AliFilteredTreeEventCuts::Merge(TCollection* list) 
{
  // Merge list of objects (needed by PROOF)
  if (!list)
  return 0;

  if (list->IsEmpty())
  return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  Int_t count=0;
  while((obj = iter->Next()) != 0) 
  {
    AliFilteredTreeEventCuts* entry = dynamic_cast<AliFilteredTreeEventCuts*>(obj);
    if (entry == 0)  
      continue; 

  count++;
  }

return count;
}

//_____________________________________________________________________________
AliFilteredTreeEventCuts::MCProcessType AliFilteredTreeEventCuts::GetEventProcessType(AliHeader* aHeader, Bool_t adebug) {
  //
  // get the process type of the event.
  //


  // Check for simple headers first

  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(aHeader->GenEventHeader());
  if (pythiaGenHeader) {
    return GetPythiaEventProcessType(pythiaGenHeader,adebug);
  }

  AliGenDPMjetEventHeader* dpmJetGenHeader = dynamic_cast<AliGenDPMjetEventHeader*>(aHeader->GenEventHeader());
  if (dpmJetGenHeader) {
    return GetDPMjetEventProcessType(dpmJetGenHeader,adebug);
  }
  

  // check for cocktail

  AliGenCocktailEventHeader* genCocktailHeader = dynamic_cast<AliGenCocktailEventHeader*>(aHeader->GenEventHeader());
  if (!genCocktailHeader) {
    printf("AliFilteredTreeEventCuts::GetProcessType : Unknown header type (not Pythia or Cocktail). \n");
    return kInvalidProcess;
  }

  TList* headerList = genCocktailHeader->GetHeaders();
  if (!headerList) {
    return kInvalidProcess;
  }

  for (Int_t i=0; i<headerList->GetEntries(); i++) {

    pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(headerList->At(i));
    if (pythiaGenHeader) {
      return GetPythiaEventProcessType(pythiaGenHeader,adebug);
    }

    dpmJetGenHeader = dynamic_cast<AliGenDPMjetEventHeader*>(headerList->At(i));
    if (dpmJetGenHeader) {
      return GetDPMjetEventProcessType(dpmJetGenHeader,adebug);
    }
  }
  return kInvalidProcess;
}

//____________________________________________________________________
AliFilteredTreeEventCuts::MCProcessType AliFilteredTreeEventCuts::GetEventProcessType(AliESDEvent* esd, AliHeader* header, AliStack* stack, DiffTreatment diffTreatment)
{
  // 
  // get process type
  //
  // diffTreatment:
  //   kMCFlags: use MC flags
  //   kUA5Cuts: SD events are those that have the MC SD flag and fulfill M^2/s < 0.05; DD from MC flag; Remainder is ND
  //   kE710Cuts: SD events are those that have the MC SD flag and fulfill 2 < M^2 < 0.05s; DD from MC flag; Remainder is ND
  //   kALICEHadronLevel: SD events are those that fulfill M^2/s < 0.05; DD from MC flag; Remainder is ND
  //

  MCProcessType mcProcessType = GetEventProcessType(header);
  
  if (diffTreatment == kMCFlags)
    return mcProcessType;
    
  if (!esd)
  {
    Printf("ERROR: AliFilteredTreeEventCuts::GetEventProcessType: diffTreatment != kMCFlags and esd == 0");
    return kInvalidProcess;
  }
    
  Float_t cms = esd->GetESDRun()->GetBeamEnergy();
  if (esd->GetESDRun()->IsBeamEnergyIsSqrtSHalfGeV())
    cms *= 2;
  //Printf("cms = %f", cms);

  if (diffTreatment == kUA5Cuts && mcProcessType == kSD)
  {
    if (IsHadronLevelSingleDiffractive(stack, cms, 0, 0.05))
      return kSD;
  }
  else if (diffTreatment == kE710Cuts && mcProcessType == kSD)
  {
    if (IsHadronLevelSingleDiffractive(stack, cms, 2. / cms / cms, 0.05))
      return kSD;
  }
  else if (diffTreatment == kALICEHadronLevel)
  {
    if (IsHadronLevelSingleDiffractive(stack, cms, 0, 0.05))
      return kSD;
  }
  
  if (mcProcessType == kSD)
    return kND;
    
  return mcProcessType;
}


AliFilteredTreeEventCuts::MCProcessType AliFilteredTreeEventCuts::GetPythiaEventProcessType(AliGenEventHeader* aHeader, Bool_t adebug) {

  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(aHeader);

  if (!pythiaGenHeader) {
    printf("AliFilteredTreeEventCuts::GetProcessType : Unknown gen Header type). \n");
    return kInvalidProcess;
  }


  Int_t pythiaType = pythiaGenHeader->ProcessType();
  fgLastProcessType = pythiaType;
  MCProcessType globalType = kInvalidProcess;  


  if (adebug) {
    printf("AliFilteredTreeEventCuts::GetProcessType : Pythia process type found: %d \n",pythiaType);
  }


  if(pythiaType==92||pythiaType==93){
    globalType = kSD;
  }
  else if(pythiaType==94){
    globalType = kDD;
  }
  //else if(pythiaType != 91){ // also exclude elastic to be sure... CKB??}
  else {
    globalType = kND;
  }
  return globalType;
}


AliFilteredTreeEventCuts::MCProcessType AliFilteredTreeEventCuts::GetDPMjetEventProcessType(AliGenEventHeader* aHeader, Bool_t adebug) {
  //
  // get the process type of the event.
  //

  // can only read pythia headers, either directly or from cocktalil header
  AliGenDPMjetEventHeader* dpmJetGenHeader = dynamic_cast<AliGenDPMjetEventHeader*>(aHeader);

  if (!dpmJetGenHeader) {
    printf("AliFilteredTreeEventCuts::GetDPMjetProcessType : Unknown header type (not DPMjet or). \n");
    return kInvalidProcess;
  }

  Int_t dpmJetType = dpmJetGenHeader->ProcessType();
  fgLastProcessType = dpmJetType;
  MCProcessType globalType = kInvalidProcess;  


  if (adebug) {
    printf("AliFilteredTreeEventCuts::GetDPMJetProcessType : DPMJet process type found: %d \n",dpmJetType);
  }


  if (dpmJetType == 1 || dpmJetType == 4) { // explicitly inelastic plus central diffraction
    globalType = kND;
  }  
  else if (dpmJetType==5 || dpmJetType==6) {
    globalType = kSD;
  }
  else if (dpmJetType==7) {
    globalType = kDD;
  }
  return globalType;
}

//____________________________________________________________________
Bool_t AliFilteredTreeEventCuts::IsHadronLevelSingleDiffractive(AliStack* stack, Float_t cms, Float_t xiMin, Float_t xiMax)
{
  //
  // return if process is single diffractive on hadron level
  // 
  // xiMax and xiMin cut on M^2/s
  //
  // Based on code from Martin Poghoysan
  //
  
  TParticle* part1 = 0;
  TParticle* part2 = 0;
  
  Double_t smallestY = 1e10;
  Double_t largestY  = -1e10;
  
  for (Int_t iParticle = 0; iParticle < stack->GetNprimary(); iParticle++)
  {
    TParticle* part = stack->Particle(iParticle);
    if (!part)
      continue;

    Int_t pdg = TMath::Abs(part->GetPdgCode());

    Int_t child1 = part->GetFirstDaughter();
    Int_t ist    = part->GetStatusCode();

    Int_t mfl  = Int_t (pdg / TMath::Power(10, Int_t(TMath::Log10(pdg))));
    if (child1 > -1 || ist != 1)
      mfl = 0; // select final state charm and beauty
    if (!(stack->IsPhysicalPrimary(iParticle) || pdg == 111 || pdg == 3212 || pdg==3124 || mfl >= 4)) 
      continue;
    Int_t imother = part->GetFirstMother();
    if (imother>0)
    {
      TParticle *partM = stack->Particle(imother);
      Int_t pdgM=TMath::Abs(partM->GetPdgCode());
      if (pdgM==111 || pdgM==3124 || pdgM==3212) 
        continue;
    }
    
    Double_t y = 0;

    // fix for problem with getting mass of particle 3124
    if (pdg != 3124)
      y = Rapidity(part->Pt(), part->Pz(), part->GetMass());
    else 
      y = Rapidity(part->Pt(), part->Pz(), 1.5195);
      
    if (y < smallestY)
    {
      smallestY = y;
      part1 = part;
    }
    
    if (y > largestY)
    {
      largestY = y;
      part2 = part;
    }
  }
  
  if (part1 == 0 || part2 == 0)
    return kFALSE;

  Int_t pdg1 = part1->GetPdgCode();
  Int_t pdg2 = part2->GetPdgCode();

  Double_t pt1 = part1->Pt();
  Double_t pt2 = part2->Pt();
  Double_t pz1 = part1->Pz();
  Double_t pz2 = part2->Pz();
  
  Double_t y1 = TMath::Abs(Rapidity(pt1, pz1, 0.938));
  Double_t y2 = TMath::Abs(Rapidity(pt2, pz2, 0.938));
  
  Int_t arm = -99999;
  if (pdg1 == 2212 && pdg2 == 2212)
  {
    if (y1 > y2) 
      arm = 0;
    else
      arm = 1;
  }
  else if (pdg1 == 2212) 
    arm = 0;
  else if (pdg2 == 2212) 
    arm = 1;

  Double_t M02s = 1. - 2 * part1->Energy() / cms;
  Double_t M12s = 1. - 2 * part2->Energy() / cms;

  if (arm == 0 && M02s > xiMin && M02s < xiMax)
    return kTRUE;
  else if (arm == 1 && M12s > xiMin && M12s < xiMax)
    return kTRUE;

  return kFALSE;
}

//____________________________________________________________________
Double_t AliFilteredTreeEventCuts::Rapidity(Double_t pt, Double_t pz, Double_t m)
{
  //
  // calculates rapidity keeping the sign in case E == pz
  //

  Double_t energy = TMath::Sqrt(pt*pt+pz*pz+m*m);
  if (energy != TMath::Abs(pz))
    return 0.5*TMath::Log((energy+pz)/(energy-pz));

  Printf("W- mt=0");
  return TMath::Sign(1.e30,pz);
}

