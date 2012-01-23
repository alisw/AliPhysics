/**************************************************************************
 * Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//
// class used to extract and store info of MC particle
//
// Author: X-M. Zhang, zhang@clermont.in2p3.fr
//                     zhangxm@iopp.ccnu.edu.cn
/////////////////////////////////////////////////////////////

#include <TParticle.h>

#include "AliMCEvent.h"
#include "AliAODMCParticle.h"
#include "AliESDMuonTrack.h"
#include "AliAODTrack.h"
#include "AliMuonInfoStoreRD.h"
#include "AliMuonInfoStoreMC.h"
#include "AliMCEvent.h"

ClassImp(AliMuonInfoStoreMC)

const TString AliMuonInfoStoreMC::fgkStdBranchName("MuonMC");
const Int_t   AliMuonInfoStoreMC::fgkSourcesN = 6;

//-----------------------------------------------------------------------------
AliMuonInfoStoreMC::AliMuonInfoStoreMC() :
AliMuonInfoStoreRD(),
fIsFull(kFALSE),
fLorentzP(),
fTrackIndex(-1),
fTrackPDGCode(0),
fSource(-1),
fNParents(0),
fOscillation(kFALSE),
fWeight(0.)
{
  //
  // default constructor
  //
  for (Int_t i=5; i--;) { fParentIndex[i] = -1; fParentPDGCode[i] = 0; }
  for (Int_t i=4; i--;) { fQuarkIndex[i]  = -1; fQuarkPDGCode[i]  = 0; }
}

//-----------------------------------------------------------------------------
AliMuonInfoStoreMC::AliMuonInfoStoreMC(AliAODTrack *trkAOD, AliMCEvent *mcEvent, Bool_t full) :
AliMuonInfoStoreRD(trkAOD),
fIsFull(full),
fLorentzP(),
fTrackIndex(-1),
fTrackPDGCode(0),
fSource(-1),
fNParents(0),
fOscillation(kFALSE),
fWeight(0.)
{
  //
  // default constructor
  //
  for (Int_t i=5; i--;) { fParentIndex[i] = -1; fParentPDGCode[i] = 0; }
  for (Int_t i=4; i--;) { fQuarkIndex[i]  = -1; fQuarkPDGCode[i]  = 0; }

  this->SetMCInfoAOD(mcEvent, trkAOD->GetLabel());
}

//-----------------------------------------------------------------------------
AliMuonInfoStoreMC::AliMuonInfoStoreMC(AliESDMuonTrack *trkESD, AliMCEvent *mcEvent, Bool_t full) :
AliMuonInfoStoreRD(trkESD),
fIsFull(full),
fLorentzP(),
fTrackIndex(-1),
fTrackPDGCode(0),
fSource(-1),
fNParents(0),
fOscillation(kFALSE),
fWeight(0.)
{
  //
  // default constructor
  //
  for (Int_t i=5; i--;) { fParentIndex[i] = -1; fParentPDGCode[i] = 0; }
  for (Int_t i=4; i--;) { fQuarkIndex[i]  = -1; fQuarkPDGCode[i]  = 0; }

  this->SetMCInfoESD(mcEvent, trkESD->GetLabel());
}

//-----------------------------------------------------------------------------
AliMuonInfoStoreMC::AliMuonInfoStoreMC(const AliMuonInfoStoreMC &src) :
AliMuonInfoStoreRD(src),
fIsFull(src.fIsFull),
fLorentzP(src.fLorentzP),
fTrackIndex(src.fTrackIndex),
fTrackPDGCode(src.fTrackPDGCode),
fSource(src.fSource),
fNParents(src.fNParents),
fOscillation(src.fOscillation),
fWeight(src.fWeight)
{
  //
  // copy constructor
  //
  for (Int_t i=5; i--;) {
    fParentIndex[i]   = src.fParentIndex[i];
    fParentPDGCode[i] = src.fParentPDGCode[i];
  }
  for (Int_t i=4; i--;) {
    fQuarkIndex[i]    = src.fQuarkIndex[i];
    fQuarkPDGCode[i]  = src.fQuarkPDGCode[i];
  }
}

//-----------------------------------------------------------------------------
AliMuonInfoStoreMC& AliMuonInfoStoreMC::operator=(const AliMuonInfoStoreMC &src)
{
  //
  // assignment constructor
  //
  if(&src==this) return *this;
  AliMuonInfoStoreRD::operator=(src);

  fIsFull       = src.fIsFull;
  fLorentzP     = src.fLorentzP;
  fTrackIndex   = src.fTrackIndex;
  fTrackPDGCode = src.fTrackPDGCode;
  fSource       = src.fSource;
  fNParents     = src.fNParents;
  fOscillation  = src.fOscillation;
  fWeight       = src.fWeight;
  for (Int_t i=5; i--;) {
    fParentIndex[i]   = src.fParentIndex[i];
    fParentPDGCode[i] = src.fParentPDGCode[i];
  }
  for (Int_t i=4; i--;) {
    fQuarkIndex[i]    = src.fQuarkIndex[i];
    fQuarkPDGCode[i]  = src.fQuarkPDGCode[i];
  }

  return *this;
}

//-----------------------------------------------------------------------------
AliMuonInfoStoreMC::~AliMuonInfoStoreMC()
{
  //
  // destructor
  //
}

//-----------------------------------------------------------------------------
void AliMuonInfoStoreMC::SetMCInfoAOD(AliMCEvent *mcEvent, Int_t label)
{
  // fill track MC info with AOD base
  fTrackIndex = label;
  if (fTrackIndex<0) { fSource=5; return; }

  AliAODMCParticle *pMC = (AliAODMCParticle*)mcEvent->GetTrack(fTrackIndex);
  fLorentzP.SetPxPyPzE(pMC->Px(), pMC->Py(), pMC->Pz(), pMC->E());

  fTrackPDGCode = pMC->PdgCode();
  if (TMath::Abs(fTrackPDGCode)!=13) { fSource=4; return; } 

  Int_t lineM = pMC->GetMother();
  if (lineM<0) { fSource=2; return; }

  Bool_t isPrimary = ((AliAODMCParticle*)mcEvent->GetTrack(lineM))->IsPrimary();
  if (!isPrimary) { fSource=3; return; }

  Int_t countP=-1, pdg=0;
  Int_t parents[10], parLine[10];
  AliAODMCParticle *mother = 0;
  while(lineM>=0){
    mother = (AliAODMCParticle*)mcEvent->GetTrack(lineM);
    pdg = mother->GetPdgCode();
    if(pdg==92 || pdg==21 || TMath::Abs(pdg)<10 || IsDiquark(pdg)) break;
    parents[++countP] = pdg;
    parLine[countP] = lineM;
    lineM = mother->GetMother();
  }
  for(Int_t i=0; i<=countP; i++){
    fParentIndex[i] = parLine[countP-i];
    fParentPDGCode[i] = parents[countP-i];
  }
  fNParents = countP + 1;

  if (fIsFull && lineM>=0) this->FillHistoryQuarksAOD(mcEvent, lineM);

  fSource = this->SelectHFMuon();
  return;
}

//-----------------------------------------------------------------------------
void AliMuonInfoStoreMC::SetMCInfoESD(AliMCEvent *mcEvent, Int_t label)
{
  // fill track MC info with ESD base
  fTrackIndex = label;
  if (fTrackIndex<0) { fSource=5; return; }

  TParticle *pMC = ((AliMCParticle*)mcEvent->GetTrack(fTrackIndex))->Particle();
  fLorentzP.SetPxPyPzE(pMC->Px(), pMC->Py(), pMC->Pz(), pMC->Energy());

  fTrackPDGCode = pMC->GetPdgCode();
  if (TMath::Abs(fTrackPDGCode)!=13) { fSource=4; return; }

  Int_t lineM = pMC->GetFirstMother();
  if (lineM<0) { fSource=2; return; }

  if (lineM>=mcEvent->Stack()->GetNprimary()) { fSource=3; return; }

  Int_t countP=-1, pdg=0;
  Int_t parents[10], parLine[10];
  TParticle *mother = 0;
  while(lineM>=0){
    mother = ((AliMCParticle*)mcEvent->GetTrack(lineM))->Particle();
    pdg = mother->GetPdgCode();
    if(pdg==92 || pdg==21 || TMath::Abs(pdg)<10 || IsDiquark(pdg)) break;
    parents[++countP] = pdg;
    parLine[countP] = lineM;
    lineM = mother->GetFirstMother();
  }
  for(Int_t i=0; i<=countP; i++){
    fParentIndex[i] = parLine[countP-i];
    fParentPDGCode[i] = parents[countP-i];
  }
  fNParents = countP + 1;

  if (fIsFull && lineM>=0) this->FillHistoryQuarksESD(mcEvent, lineM);

  fSource = this->SelectHFMuon();
  return;
}

//-----------------------------------------------------------------------------
void AliMuonInfoStoreMC::FillHistoryQuarksAOD(AliMCEvent* const mcEvent, Int_t lineM)
{
  // method in $ALICE_ROOT/MUON/AliMUONTrackLight.cxx 

  if (lineM<0) return;
  Int_t countP=-1, pdg=0;
  AliAODMCParticle *mother = 0;
  while(lineM>=0){
    mother = (AliAODMCParticle*)mcEvent->GetTrack(lineM);
    pdg = mother->GetPdgCode();
    fQuarkIndex[++countP] = lineM;
    fQuarkPDGCode[countP] = pdg;
    lineM = mother->GetMother();
  }

  // for PYTHIA checking
  countP = 1;
  for(Int_t par=0; par<4; par++) {
    if(TMath::Abs(this->QuarkPDGCode(par))<6) { countP=par; break; }
  }
  if(this->QuarkIndex(countP)>-1 && (this->ParentFlavour(0)==4 || this->ParentFlavour(0)==5)) {
    if(this->ParentFlavour(0)!=TMath::Abs(this->QuarkPDGCode(countP))) {
      AliWarning(Form("quark flavour of parent and that of quark do not correspond: %d %d --> correcting\n",
                 this->ParentFlavour(0), TMath::Abs(this->QuarkPDGCode(countP))));

      pdg = this->QuarkPDGCode(countP);
      Int_t line = this->QuarkIndex(countP);
      this->ResetQuarkInfo();
      while(TMath::Abs(pdg)!=this->ParentFlavour(0)) {
        pdg = ((AliAODMCParticle*)mcEvent->GetTrack(++line))->GetPdgCode();
      }
      while(line>=0){
        mother = (AliAODMCParticle*)mcEvent->GetTrack(line);
        pdg = mother->GetPdgCode();
        fQuarkIndex[countP] = line;
        fQuarkPDGCode[countP++] = pdg;
        line = mother->GetMother();
      }
    }
  }
  return;
}

//-----------------------------------------------------------------------------
void AliMuonInfoStoreMC::FillHistoryQuarksESD(AliMCEvent* const mcEvent, Int_t lineM)
{
  // method in $ALICE_ROOT/MUON/AliMUONTrackLight.cxx 

  if (lineM<0) return;
  Int_t countP=-1, pdg=0;
  TParticle *mother = 0;
  while(lineM>=0){
    mother = ((AliMCParticle*)mcEvent->GetTrack(lineM))->Particle();
    pdg = mother->GetPdgCode();
    fQuarkIndex[++countP] = lineM;
    fQuarkPDGCode[countP] = pdg;
    lineM = mother->GetFirstMother();
  }

  // for PYTHIA checking
  countP = 1;
  for(Int_t par=0; par<4; par++) {
    if(TMath::Abs(this->QuarkPDGCode(par))<6) { countP=par; break; }
  }
  if(this->QuarkIndex(countP)>-1 && (this->ParentFlavour(0)==4 || this->ParentFlavour(0)==5)) {
    if(this->ParentFlavour(0)!=TMath::Abs(this->QuarkPDGCode(countP))) {
      AliWarning(Form("quark flavour of parent and that of quark do not correspond: %d %d --> correcting\n",
                 this->ParentFlavour(0), TMath::Abs(this->QuarkPDGCode(countP))));

      pdg = this->QuarkPDGCode(countP);
      Int_t line = this->QuarkIndex(countP);
      this->ResetQuarkInfo();
      while(TMath::Abs(pdg)!=this->ParentFlavour(0)) {
        pdg = ((AliMCParticle*)mcEvent->GetTrack(++lineM))->Particle()->GetPdgCode();
      }
      while(line>=0){
        mother = ((AliMCParticle*)mcEvent->GetTrack(++lineM))->Particle();
        pdg = mother->GetPdgCode();
        fQuarkIndex[countP] = line;
        fQuarkPDGCode[countP++] = pdg;
        line = mother->GetFirstMother();
      }
    }
  }
  return;
}

//-----------------------------------------------------------------------------
Int_t AliMuonInfoStoreMC::SelectHFMuon()
{
  // set info of muon from HF

  Int_t flv = ParentFlavour(0);
  if (flv!=4 && flv!=5) return 2;

  Bool_t isRes = kFALSE;
  Int_t i=0, nparents=this->ParentsN();
  while (i<nparents && !isRes) {
    isRes = IsMotherAResonance(i++);
  }

  if (isRes) return 2;
  if (flv==5) return 0;
  else return 1;
}

//-----------------------------------------------------------------------------
Bool_t AliMuonInfoStoreMC::IsDiquark(Int_t pdg)
{
  // copy from $ALICE_ROOT/MUON/AliMUONTrackLight.cxx
  pdg = TMath::Abs(pdg);
  if(pdg>1000 && (pdg%100)<10) return kTRUE;
  else return kFALSE;
}

//-----------------------------------------------------------------------------
void AliMuonInfoStoreMC::ResetQuarkInfo()
{
  // copy from $ALICE_ROOT/MUON/AliMUONTrackLight.cxx
  for(Int_t pos=1; pos<4; pos++) {
    fQuarkIndex[pos] = -1;
    fQuarkPDGCode[pos] = 0;
  }
  return;
}

//-----------------------------------------------------------------------------
Int_t AliMuonInfoStoreMC::ParentFlavour(Int_t i) const
{
  // copy from $ALICE_ROOT/MUON/AliMUONTrackLight.cxx
  Int_t pdg = ParentPDGCode(i);
  pdg = TMath::Abs(pdg/100);
  if(pdg>9) pdg /= 10;
  return pdg;
}

//-----------------------------------------------------------------------------
Bool_t AliMuonInfoStoreMC::IsMotherAResonance(Int_t i) const
{
  // copy from $ALICE_ROOT/MUON/AliMUONTrackLight.cxx
  Int_t pdg = ParentPDGCode(i);
  Int_t id=pdg%100000;
  return (!((id-id%10)%110));
}
