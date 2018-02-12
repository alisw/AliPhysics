/*
 * AliFemtoDreamv0.cxx
 *
 *  Created on: Dec 12, 2017
 *      Author: gu74req
 */

#include <vector>
#include "AliAODMCParticle.h"
#include "AliFemtoDreamv0.h"
#include "TClonesArray.h"

ClassImp(AliFemtoDreamv0)
AliFemtoDreamv0::AliFemtoDreamv0()
:AliFemtoDreamBasePart()
,fOnlinev0(false)
,fHasDaughter(false)
,fpDaug(new AliFemtoDreamTrack())
,fnDaug(new AliFemtoDreamTrack())
,fv0Mass(0)
,fdcav0Daug(0)
,fdcaPrim(0)
,fdcaPrimPos(0)
,fdcaPrimNeg(0)
,flenDecay(0)
,fTransRadius(0)
{
  for (int i=0;i<3;++i) {
    fv0Vtx[i]=0;
  }
}

AliFemtoDreamv0::~AliFemtoDreamv0() {
  if (fpDaug) {
    delete fpDaug;
  }
  if (fnDaug) {
    delete fnDaug;
  }
}

void AliFemtoDreamv0::Setv0(AliAODEvent *evt,AliAODv0* v0) {
  if (!fGTI) {
    AliFatal("no GTI Array set");
  }
  if (!v0) {
    AliFatal("SetProng No v0 to work with");
  }
  Reset();
  if (v0->GetNProngs()==2&&v0->GetNDaughters()==2) {
    fIsReset=false;
    if (v0->GetOnFlyStatus()) {
      this->fOnlinev0=true;
    } else {
      this->fOnlinev0=false;
    }
    this->SetMotherInfo(evt,v0);
    if (fIsMC) {
      this->SetMCMotherInfo(evt,v0);
    }
    this->SetDaughter(v0);
  } else {
    this->SetUse(false);
  }
}

void AliFemtoDreamv0::SetDaughter(AliAODv0 *v0) {
  if (v0->GetPosID()>=fTrackBufferSize||v0->GetNegID()>=fTrackBufferSize) {
    AliFatal("fGTI too small");
  } else {
    fpDaug->SetGlobalTrackInfo(fGTI,fTrackBufferSize);
    fnDaug->SetGlobalTrackInfo(fGTI,fTrackBufferSize);
    if (fGTI[v0->GetPosID()]&&fGTI[v0->GetNegID()]) {
      if (fGTI[v0->GetPosID()]->Charge() > 0 && fGTI[v0->GetNegID()]->Charge() < 0) {
        fnDaug->SetTrack(fGTI[v0->GetNegID()]);
        fpDaug->SetTrack(fGTI[v0->GetPosID()]);
        this->SetDaughterInfo(v0);
        this->fHasDaughter=true;
      } else if(fGTI[v0->GetPosID()]->Charge() < 0 && fGTI[v0->GetNegID()]->Charge() > 0) {
        fnDaug->SetTrack(fGTI[v0->GetPosID()]);
        fpDaug->SetTrack(fGTI[v0->GetNegID()]);
        this->SetDaughterInfo(v0);
        this->fHasDaughter=true;
      } else {
        this->fHasDaughter=false;
      }
    } else {
      this->fHasDaughter=false;
    }
  }
}

void AliFemtoDreamv0::SetDaughterInfo(AliAODv0 *v0) {
  this->SetEta(fnDaug->GetEta().at(0));
  this->SetEta(fpDaug->GetEta().at(0));

  this->SetTheta(fnDaug->GetTheta().at(0));
  this->SetTheta(fpDaug->GetTheta().at(0));

  this->SetPhi(fnDaug->GetPhi().at(0));
  this->SetPhi(fpDaug->GetPhi().at(0));

  this->SetIDTracks(fnDaug->GetIDTracks().at(0));
  this->SetIDTracks(fpDaug->GetIDTracks().at(0));

  this->SetCharge(fnDaug->GetCharge().at(0));
  this->SetCharge(fpDaug->GetCharge().at(0));
  //here we have to set the momentum from the v0. The Momentum of the Global
  //track is different to the one of the daughter track and is also
  //used in the official v0 class for the invarian mass
  //calculation.
  fnDaug->SetMomentum(v0->PxProng(1), v0->PyProng(1), v0->PzProng(1));
  fpDaug->SetMomentum(v0->PxProng(0), v0->PyProng(0), v0->PzProng(0));

  if(fIsMC) {
    if (fnDaug->IsSet()) {
      this->SetMCTheta(fnDaug->GetMCTheta().at(0));
      this->SetMCPhi(fnDaug->GetMCPhi().at(0));
    }
    if (fpDaug->IsSet()) {
      this->SetMCTheta(fpDaug->GetMCTheta().at(0));
      this->SetMCPhi(fpDaug->GetMCPhi().at(0));
    }
  }
}

void AliFemtoDreamv0::SetMotherInfo(AliAODEvent *evt,AliAODv0 *v0) {
  this->SetCharge(v0->GetCharge());
  this->SetPt(v0->Pt());
  this->SetMomentum(v0->Px(),v0->Py(),v0->Pz());
  this->SetEta(v0->Eta());
  this->SetPhi(v0->Phi());
  this->SetTheta(v0->Theta());
  double xvP = evt->GetPrimaryVertex()->GetX();
  double yvP = evt->GetPrimaryVertex()->GetY();
  double zvP = evt->GetPrimaryVertex()->GetZ();
  double vecTarget[3]={xvP,yvP,zvP};
  v0->GetXYZ(fv0Vtx);
  this->fdcav0Daug=v0->DcaV0Daughters();
  this->fdcaPrim=v0->DcaV0ToPrimVertex();
  this->fdcaPrimPos=v0->DcaPosToPrimVertex();
  this->fdcaPrimNeg=v0->DcaNegToPrimVertex();
  this->flenDecay=v0->DecayLengthV0(vecTarget);
  this->fCPA=v0->CosPointingAngle(vecTarget);
  this->fTransRadius=v0->DecayLengthXY(vecTarget);
}

void AliFemtoDreamv0::SetMCMotherInfo(AliAODEvent *evt,AliAODv0 *v0) {
  TClonesArray *mcarray = dynamic_cast<TClonesArray*>(
      evt->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!mcarray) {
    AliFatal("No MC Array found");
  }
  int PDGDaug[2];
  PDGDaug[0]=TMath::Abs(fpDaug->GetPDGCode());
  PDGDaug[1]=TMath::Abs(fnDaug->GetPDGCode());
  int label=v0->MatchToMC(TMath::Abs(fPDGCode),mcarray,2,PDGDaug);
  if (label<0) {
    //label will be -1 if there was not 'real' candidate for matching,
    //therefore we have the case of a contamination/background v0
    this->SetParticleOrigin(AliFemtoDreamBasePart::kContamination);
  } else {
    AliAODMCParticle* mcPart=(AliAODMCParticle*)mcarray->At(label);
    if (!mcPart) {
      this->SetUse(false);
    } else {
      this->SetMCPDGCode(mcPart->GetPdgCode());
      double mcMom[3]={0.,0.,0.};
      mcPart->PxPyPz(mcMom);
      this->SetMCMomentum(mcMom[0],mcMom[1],mcMom[2]);
      this->SetMCPt(mcPart->Pt());
      this->SetMCPhi(mcPart->Phi());
      this->SetMCTheta(mcPart->Theta());
      if (mcPart->IsPhysicalPrimary()&&!(mcPart->IsSecondaryFromWeakDecay())) {
        this->SetParticleOrigin(AliFemtoDreamBasePart::kPhysPrimary);
      } else if (mcPart->IsSecondaryFromWeakDecay()&&
          !(mcPart->IsSecondaryFromMaterial())) {
        this->SetParticleOrigin(AliFemtoDreamBasePart::kWeak);
        this->SetPDGMotherWeak(((AliAODMCParticle*)mcarray->At(
            mcPart->GetMother()))->PdgCode());
      } else if (mcPart->IsSecondaryFromMaterial()) {
        this->SetParticleOrigin(AliFemtoDreamBasePart::kMaterial);
      } else {
        this->SetParticleOrigin(AliFemtoDreamBasePart::kUnknown);
      }
    }
  }
}
void AliFemtoDreamv0::Reset() {
  if (!fIsReset) {
    fOnlinev0=false;
    fHasDaughter=false;
    //daughters don't need to be reset, are reset while setting a new track
    fv0Mass=0;
    fv0Vtx[0]=99;
    fv0Vtx[1]=99;
    fv0Vtx[2]=99;
    fdcav0Daug=99;
    fdcaPrim=0;
    fdcaPrimPos=0;
    fdcaPrimNeg=0;
    flenDecay=0;
    fTransRadius=0;
    fP.SetXYZ(0,0,0);
    fMCP.SetXYZ(0,0,0);
    fPt=0;
    fMCPt=0;
    fMCPt=0;
    fP_TPC=0;
    fEta.clear();
    fTheta.clear();
    fMCTheta.clear();
    fPhi.clear();
    fMCPhi.clear();
    fIDTracks.clear();
    fCharge.clear();
    fCPA=0;
    fOrigin=AliFemtoDreamBasePart::kUnknown;
    //we don't want to reset the fPDGCode
    fMCPDGCode=0;
    fPDGMotherWeak=0;
    //we don't want to reset isMC
    fUse=false;
    fIsSet=true;
    fIsReset=true;
  }
}
