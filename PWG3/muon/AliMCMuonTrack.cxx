#include <TParticle.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAODMCParticle.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDMuonCluster.h"
#include "AliAODTrack.h"
#include "AliAODMuonTrack.h"
#include "AliMCMuonTrack.h"

ClassImp(AliMCMuonTrack)

const Double_t AliMCMuonTrack::fgkSigmaCut = 10.;

//-----------------------------------------------------------------------------
AliMCMuonTrack::AliMCMuonTrack() :
AliAODMuonTrack(),
fIsFull(kFALSE),
fPGen(),
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
  for (Int_t i=0; i<fgkNParentsMax; i++) {
    fParentIndex[i] = -1;
    fParentPDGCode[i] = 0;
  }
  for (Int_t i=0; i<4; i++) {
    fQuarkIndex[i] = -1;
    fQuarkPDGCode[i] = 0;
  }
}

//-----------------------------------------------------------------------------
AliMCMuonTrack::AliMCMuonTrack(AliAODTrack *trkAOD, TClonesArray *mcClArr, Bool_t full) :
AliAODMuonTrack(trkAOD),
fIsFull(full),
fPGen(),
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
  AliAODMCParticle *pMC = this->FindTrackRef(trkAOD, mcClArr);
  if (pMC) this->SetMCInfo(pMC, mcClArr);
}

//-----------------------------------------------------------------------------
AliMCMuonTrack::AliMCMuonTrack(AliESDMuonTrack *trkESD, AliESDEvent *esd, AliMCEventHandler *mcH, Bool_t full) :
AliAODMuonTrack(trkESD),
fIsFull(full),
fPGen(),
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
  TParticle *pMC = this->FindTrackRef(trkESD, esd, mcH);  // do nothing when running on official train
  if (pMC) this->SetMCInfo(pMC, mcH);
  // MC infomation is not implement when running on the official train for they depend on the MUON module.
  // In the case of private running, one can get the information from ESD by uncommenting the definitions
  // of FindTrackRef() and CovESDtoMuonTrack()
  // and including the following head files in both AliMCMuonTrack.h and this file,
  // #include "AliMUONRecoCheck.h"
  // #include "AliMUONVTrackStore.h"
  // #include "AliMUONVClusterStore.h"
  // #include "AliMUONTrackParam.h"
  // #include "AliMUONTrack.h"
  // #include "AliMUONVCluster.h"
  // #include "AliMUONESDInterface.h"
}



//-----------------------------------------------------------------------------
AliMCMuonTrack::~AliMCMuonTrack()
{
  //
  // destructor
  //
}

//-----------------------------------------------------------------------------
AliAODMCParticle* AliMCMuonTrack::FindTrackRef(AliAODTrack *trkAOD, TClonesArray *mcClArr)
{
  AliAODMCParticle *pMC = 0;
  fTrackIndex = trkAOD->GetLabel();
  if (fTrackIndex>=0)
    pMC = (AliAODMCParticle*)mcClArr->At(fTrackIndex);
  return pMC;
}

//-----------------------------------------------------------------------------
TParticle* AliMCMuonTrack::FindTrackRef(AliESDMuonTrack *trkESD, AliESDEvent *esd, AliMCEventHandler *mcH)
{
  TParticle *pMCRef = 0;

  /*AliMUONTrack *trkMuon = CovESDtoMuonTrack(*trkESD);
  AliMUONRecoCheck rc(esd,mcH);
  AliMUONVTrackStore* trkRefArr = rc.TrackRefs(-1);
  TIter next(trkRefArr->CreateIterator());
  AliMUONTrack* trkRef = 0;
  while ((trkRef=static_cast<AliMUONTrack*>(next()))) {
    Bool_t trkCompArr[10];
    Int_t nMatchClusters = trkMuon->CompatibleTrack(trkRef, fgkSigmaCut, trkCompArr);
    Double_t matchClusterFrac = ((Double_t)nMatchClusters) / ((Double_t)trkMuon->GetNClusters());
    if ((trkCompArr[0] || trkCompArr[1] || trkCompArr[2] || trkCompArr[3]) &&
        (trkCompArr[6] || trkCompArr[7] || trkCompArr[8] || trkCompArr[9]) && matchClusterFrac>0.5) {
      fTrackIndex = trkRef->GetUniqueID();
      trkRefArr->Remove(*trkRef);
      break;
    }
  }
  if (fTrackIndex>=0) pMCRef = mcH->MCEvent()->Stack()->Particle(fTrackIndex);*/

  return pMCRef;
}

//-----------------------------------------------------------------------------
void AliMCMuonTrack::SetMCInfo(AliAODMCParticle *pMC, TClonesArray *mcClArr)
{
  fPGen.SetPxPyPzE(pMC->Px(), pMC->Py(), pMC->Pz(), pMC->E());
  fTrackPDGCode = pMC->GetPdgCode();
  if (TMath::Abs(fTrackPDGCode)!=13) {
    fSource = 4;
    return;
  } 

  Int_t lineM = pMC->GetMother();
  if (lineM<0) {
    fSource = 2;
    return;
  }

  Bool_t isPrimary =
    ((AliAODMCParticle*)mcClArr->At(lineM))->IsPrimary();
  if (!isPrimary) {
    fSource = 3;
    return;
  }

  this->FillHistoryParents(lineM, mcClArr);
  fSource = this->SelectHFMuon();
  return;
}

//-----------------------------------------------------------------------------
void AliMCMuonTrack::SetMCInfo(TParticle *pMC, AliMCEventHandler *mcH)
{
  fPGen.SetPxPyPzE(pMC->Px(), pMC->Py(), pMC->Pz(), pMC->Energy());
  fTrackPDGCode = pMC->GetPdgCode();
  if (TMath::Abs(fTrackPDGCode)!=13) {
    fSource = 4;
    return;
  }

  Int_t lineM = pMC->GetFirstMother();
  if (lineM<0) {
    fSource = 2;
    return;
  }

  AliStack *stack = mcH->MCEvent()->Stack();
  if (lineM>=stack->GetNprimary()) {
    fSource = 3;
    return;
  }

  this->FillHistoryParents(lineM, stack);
  fSource = this->SelectHFMuon();
  return;
}

//-----------------------------------------------------------------------------
void AliMCMuonTrack::FillHistoryParents(Int_t lineM, TClonesArray *mcClArr)
{
  Int_t countP=-1, pdg=0;
  Int_t parents[10], parLine[10];
  AliAODMCParticle *mother = 0;
  while(lineM>=0){
    mother = (AliAODMCParticle*)mcClArr->At(lineM);
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

  if (fIsFull && lineM>=0) this->FillHistoryQuarks(lineM, mcClArr);
  return;
}

//-----------------------------------------------------------------------------
void AliMCMuonTrack::FillHistoryParents(Int_t lineM, AliStack *stack)
{
  Int_t countP=-1, pdg=0;
  Int_t parents[10], parLine[10];
  TParticle *mother = 0;
  while(lineM>=0){
    mother = stack->Particle(lineM);
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

  if (fIsFull && lineM>=0) this->FillHistoryQuarks(lineM, stack);
  return;
}

//-----------------------------------------------------------------------------
void AliMCMuonTrack::FillHistoryQuarks(Int_t lineM, TClonesArray *mcClArr)
{
  // method in $ALICE_ROOT/MUON/AliMUONTrackLight.cxx 
  if (lineM<0) return;
  Int_t countP=-1, pdg=0;
  AliAODMCParticle *mother = 0;
  while(lineM>=0){
    mother = (AliAODMCParticle*)mcClArr->At(lineM);
    pdg = mother->GetPdgCode();
    fQuarkIndex[++countP] = lineM;
    fQuarkPDGCode[countP] = pdg;
    lineM = mother->GetMother();
  }

  // for PYTHIA checking
  countP = 1;
  for(Int_t par=0; par<4; par++) {
    if(TMath::Abs(this->GetQuarkPDGCode(par))<6) { countP=par; break; }
  }
  if(this->GetQuarkIndex(countP)>-1 && (this->GetParentFlavour(0)==4 || this->GetParentFlavour(0)==5)) {
    if(this->GetParentFlavour(0)!=TMath::Abs(this->GetQuarkPDGCode(countP))) {
      AliWarning(Form("quark flavour of parent and that of quark do not correspond: %d %d --> correcting\n",
                 this->GetParentFlavour(0), TMath::Abs(this->GetQuarkPDGCode(countP))));

      pdg = this->GetQuarkPDGCode(countP);
      Int_t line = this->GetQuarkIndex(countP);
      this->ResetQuarkInfo();
      while(TMath::Abs(pdg)!=this->GetParentFlavour(0)) {
        pdg = ((AliAODMCParticle*)mcClArr->At(++line))->GetPdgCode();
      }
      while(line>=0){
        mother = (AliAODMCParticle*)mcClArr->At(line);
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
void AliMCMuonTrack::FillHistoryQuarks(Int_t lineM, AliStack *stack)
{
  // method in $ALICE_ROOT/MUON/AliMUONTrackLight.cxx 
  if (lineM<0) return;
  Int_t countP=-1, pdg=0;
  TParticle *mother = 0;
  while(lineM>=0){
    mother = stack->Particle(lineM);
    pdg = mother->GetPdgCode();
    fQuarkIndex[++countP] = lineM;
    fQuarkPDGCode[countP] = pdg;
    lineM = mother->GetFirstMother();
  }

  // for PYTHIA checking
  countP = 1;
  for(Int_t par=0; par<4; par++) {
    if(TMath::Abs(this->GetQuarkPDGCode(par))<6) { countP=par; break; }
  }
  if(this->GetQuarkIndex(countP)>-1 && (this->GetParentFlavour(0)==4 || this->GetParentFlavour(0)==5)) {
    if(this->GetParentFlavour(0)!=TMath::Abs(this->GetQuarkPDGCode(countP))) {
      AliWarning(Form("quark flavour of parent and that of quark do not correspond: %d %d --> correcting\n",
                 this->GetParentFlavour(0), TMath::Abs(this->GetQuarkPDGCode(countP))));

      pdg = this->GetQuarkPDGCode(countP);
      Int_t line = this->GetQuarkIndex(countP);
      this->ResetQuarkInfo();
      while(TMath::Abs(pdg)!=this->GetParentFlavour(0)) {
        pdg = stack->Particle(++line)->GetPdgCode();
      }
      while(line>=0){
        mother = stack->Particle(line);
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
/*AliMUONTrack* AliMCMuonTrack::CovESDtoMuonTrack(AliESDMuonTrack &trkESD)
{
  // method in $ALICE_ROOT/PWG3/muondep/AliAnalysisTaskESDMCLabelAddition.cxx
  AliMUONTrack *trkMuon = new AliMUONTrack();
  if(!trkESD.ClustersStored()) return trkMuon;

  AliMUONTrackParam param;
  AliMUONESDInterface::GetParamAtFirstCluster(trkESD, param);
  AliMUONESDInterface::GetParamCov(trkESD, param);
  AliMUONVClusterStore *cStore = AliMUONESDInterface::NewClusterStore();
  AliMUONVCluster *cluster = cStore->CreateCluster(0,0,0);

  AliESDMuonCluster *esdCluster = (AliESDMuonCluster*)trkESD.GetClusters().First();
  while (esdCluster) {
    AliMUONESDInterface::ESDToMUON(*esdCluster, *cluster);
    param.SetZ(cluster->GetZ());
    trkMuon->AddTrackParamAtCluster(param, *cluster, kTRUE);
    esdCluster = (AliESDMuonCluster*)trkESD.GetClusters().After(esdCluster);
  }

  delete cluster; cluster = 0;
  delete cStore;  cStore = 0;
  return trkMuon;
}*/

//-----------------------------------------------------------------------------
Int_t AliMCMuonTrack::SelectHFMuon()
{
  Int_t flv = GetParentFlavour(0);
  if (flv!=4 && flv!=5) return 2;

  Bool_t isRes = kFALSE;
  Int_t i=0, nparents=this->GetNParents();
  while (i<nparents && !isRes) {
    isRes = IsMotherAResonance(i++);
  }

  if (isRes) return 2;
  if (flv==5) return 0;
  else return 1;
}

//-----------------------------------------------------------------------------
Bool_t AliMCMuonTrack::IsDiquark(Int_t pdg)
{
  // copy from $ALICE_ROOT/MUON/AliMUONTrackLight.cxx
  pdg = TMath::Abs(pdg);
  if(pdg>1000 && (pdg%100)<10) return kTRUE;
  else return kFALSE;
}

//-----------------------------------------------------------------------------
void AliMCMuonTrack::ResetQuarkInfo()
{
  for(Int_t pos=1; pos<4; pos++) {
    fQuarkIndex[pos] = -1;
    fQuarkPDGCode[pos] = 0;
  }
  return;
}

//-----------------------------------------------------------------------------
Int_t AliMCMuonTrack::GetParentFlavour(Int_t i) const
{
  Int_t pdg = GetParentPDGCode(i);
  pdg = TMath::Abs(pdg/100);
  if(pdg>9) pdg /= 10;
  return pdg;
}

//-----------------------------------------------------------------------------
Bool_t AliMCMuonTrack::IsMotherAResonance(Int_t i) const
{
  // copy from $ALICE_ROOT/MUON/AliMUONTrackLight.cxx
  Int_t pdg = GetParentPDGCode(i);
  Int_t id=pdg%100000;
  return (!((id-id%10)%110));
}
