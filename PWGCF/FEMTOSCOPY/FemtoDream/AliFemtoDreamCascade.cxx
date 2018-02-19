/*
 * AliFemtoDreamCascade.cxx
 *
 *  Created on: Jan 11, 2018
 *      Author: gu74req
 */

#include "AliFemtoDreamCascade.h"
#include "TMath.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODMCParticle.h"
ClassImp(AliFemtoDreamCascade)
AliFemtoDreamCascade::AliFemtoDreamCascade()
:AliFemtoDreamBasePart()
,fPosDaug(new AliFemtoDreamTrack())
,fNegDaug(new AliFemtoDreamTrack())
,fBach(new AliFemtoDreamTrack())
,fXiMass(0)
,fOmegaMass(0)
,fDCAXiDaug(0)
,fTransRadius(0)
,fRapXi(0)
,fAlphaXi(0)
,fPtArmXi(0)
,fXiLength(0)
,fOmegaLength(0)
,fMassv0(0)
,fv0DCADaug(0)
,fv0DCAPrimVtx(0)
,fv0Momentum()
,fDCABachPrimVtx(0)
,fv0TransRadius(0)
,fv0CPA(0)
,fv0PosToPrimVtx(0)
,fv0NegToPrimVtx(0)
,fv0ToXiPointAngle(0)
,fv0Length(0)
,fDCAv0Xi(0)
{
}

AliFemtoDreamCascade::~AliFemtoDreamCascade() {
  if (fPosDaug) {
    delete fPosDaug;
  }
  if (fNegDaug) {
    delete fNegDaug;
  }
  if (fBach) {
    delete fBach;
  }
}

void AliFemtoDreamCascade::SetCascade(AliAODEvent *evt,AliAODcascade *casc) {
  //xi business
  Reset();
  fIsReset=false;
  double PrimVtx[3]={99.,99.,99};
  double decayPosXi[3]={casc->DecayVertexXiX(),casc->DecayVertexXiY(),
      casc->DecayVertexXiZ()};

  fXiMass=casc->MassXi();
  fOmegaMass=casc->MassOmega();
  fDCAXiDaug=casc->DcaXiDaughters();
  evt->GetPrimaryVertex()->GetXYZ(PrimVtx);
  fCPA=casc->CosPointingAngleXi(PrimVtx[0],PrimVtx[1],PrimVtx[2]);
  fTransRadius=TMath::Sqrt(
      decayPosXi[0]*decayPosXi[0]+decayPosXi[1]*decayPosXi[1]);
  this->SetCharge(casc->ChargeXi());
  this->SetMomentum(casc->MomXiX(),casc->MomXiY(),casc->MomXiZ());
  this->SetPt(casc->MomXiX()*casc->MomXiX()+casc->MomXiY()*casc->MomXiY());
  fRapXi=casc->RapXi();
  this->SetEta(casc->Eta());
  this->SetTheta(casc->Theta());
  this->SetPhi(casc->Phi());
  fAlphaXi=casc->AlphaXi();
  fPtArmXi=casc->PtArmXi();
  fXiLength=TMath::Sqrt(TMath::Power((decayPosXi[0]-PrimVtx[0]),2)+
                        TMath::Power((decayPosXi[1]-PrimVtx[1]),2)+
                        TMath::Power((decayPosXi[2]-PrimVtx[2]),2));
  fOmegaLength=fXiLength;
  double XiPDGMass=1.321;
  double OmegaPDGMass=1.672;
  if (this->GetMomentum().Mag()>0) {
    fXiLength=fXiLength*XiPDGMass/this->GetMomentum().Mag()>0;
    fOmegaLength=fOmegaLength*OmegaPDGMass/this->GetMomentum().Mag()>0;
  } else {
    fXiLength=-1.;
  }

  //From here daughter business
  AliAODTrack *nTrackXi=dynamic_cast<AliAODTrack*>(casc->GetDaughter(1));
  AliAODTrack *pTrackXi=dynamic_cast<AliAODTrack*>(casc->GetDaughter(0));
  AliAODTrack *bachTrackXi=dynamic_cast<AliAODTrack*>(
      casc->GetDecayVertexXi()->GetDaughter(0));
  fNegDaug->SetTrack(nTrackXi);
  fPosDaug->SetTrack(pTrackXi);
  fBach->SetTrack(bachTrackXi);

  fNegDaug->SetMomentum(casc->MomNegX(),casc->MomNegY(),casc->MomNegZ());
  fPosDaug->SetMomentum(casc->MomPosX(),casc->MomPosY(),casc->MomPosZ());
  fBach->SetMomentum(casc->MomBachX(),casc->MomBachY(),casc->MomBachZ());

  this->SetIDTracks(nTrackXi->GetID());
  this->SetIDTracks(pTrackXi->GetID());
  this->SetIDTracks(bachTrackXi->GetID());

  this->SetEta(nTrackXi->Eta());
  this->SetEta(pTrackXi->Eta());
  this->SetEta(bachTrackXi->Eta());

  this->SetTheta(nTrackXi->Theta());
  this->SetTheta(pTrackXi->Theta());
  this->SetTheta(bachTrackXi->Theta());

  this->SetPhi(nTrackXi->Phi());
  this->SetPhi(pTrackXi->Phi());
  this->SetPhi(bachTrackXi->Phi());

  this->SetCharge(nTrackXi->Charge());
  this->SetCharge(pTrackXi->Charge());
  this->SetCharge(bachTrackXi->Charge());

  //v0 business
  if (casc->ChargeXi()<0) {//Xi minus
    fMassv0=casc->MassLambda();
  } else {//Xi plus case
    fMassv0=casc->MassAntiLambda();
  }
  fv0Momentum.SetXYZ(casc->MomV0X(),casc->MomV0Y(),casc->MomV0Z());
  fv0DCADaug=casc->DcaV0Daughters();
  fv0DCAPrimVtx=casc->DcaV0ToPrimVertex();
  fDCABachPrimVtx=casc->DcaBachToPrimVertex();
  double decayPosV0[3]={casc->DecayVertexV0X(),
      casc->DecayVertexV0Y(),casc->DecayVertexV0Z()};
  fv0TransRadius=TMath::Sqrt(decayPosV0[0]*decayPosV0[0]+
                             decayPosV0[1]*decayPosV0[1]);
  fv0CPA=casc->CosPointingAngle(evt->GetPrimaryVertex());
  fv0PosToPrimVtx=casc->DcaPosToPrimVertex();
  fv0NegToPrimVtx=casc->DcaNegToPrimVertex();
  fv0ToXiPointAngle=casc->CosPointingAngle(casc->GetDecayVertexXi());
  double LamPDGMass=1.115683;
  fv0Length=TMath::Sqrt(
      TMath::Power((decayPosV0[0]-decayPosXi[0]),2)+
      TMath::Power((decayPosV0[1]-decayPosXi[1]),2)+
      TMath::Power((decayPosV0[2]-decayPosXi[2]),2));
  //  Float_t lctauV0 = -1.;
  if (fv0Momentum.Mag()>0) {
    fv0Length=fv0Length*LamPDGMass/fv0Momentum.Mag();
  } else {
    fv0Length=-1.;
  }
  fDCAv0Xi=TMath::Sqrt(
      TMath::Power((decayPosV0[0]-decayPosXi[0]),2)+
      TMath::Power((decayPosV0[1]-decayPosXi[1]),2));
  if (fIsMC) {
    SetMCMotherInfo(evt,casc);
  }
}

void AliFemtoDreamCascade::Reset() {
  if (!fIsReset) {
    fXiMass=0;
    fOmegaMass=0;
    fDCAXiDaug=0;
    fTransRadius=0;
    fRapXi=-99;
    fAlphaXi=99;
    fPtArmXi=99;
    fXiLength=0;
    fOmegaLength=0;
    fMassv0=0;
    fv0DCADaug=0;
    fv0DCAPrimVtx=0;
    fv0Momentum.SetXYZ(0,0,0);
    fDCABachPrimVtx=0;
    fv0TransRadius=0;
    fv0CPA=0;
    fv0PosToPrimVtx=0;
    fv0NegToPrimVtx=0;
    fv0ToXiPointAngle=0;
    fv0Length=0;
    fDCAv0Xi=0;
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
    fIsSet=false;
    fIsReset=true;
  }
}

void AliFemtoDreamCascade::SetMCMotherInfo(
    AliAODEvent *evt, AliAODcascade *casc) {
  TClonesArray *mcarray = dynamic_cast<TClonesArray*>(
      evt->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!mcarray) {
    AliFatal("No MC Array found");
  }
  int PDGDaug[3];
  PDGDaug[0]=TMath::Abs(fPosDaug->GetPDGCode());
  PDGDaug[1]=TMath::Abs(fNegDaug->GetPDGCode());
  PDGDaug[2]=TMath::Abs(fBach->GetPDGCode());
  int label=casc->MatchToMC(TMath::Abs(fPDGCode),mcarray,2,PDGDaug);
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
