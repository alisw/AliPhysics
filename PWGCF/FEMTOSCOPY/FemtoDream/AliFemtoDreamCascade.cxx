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
,fRapOmega(0)
,fAlphaXi(0)
,fPtArmXi(0)
,fDCAXiPrimVtx(0)
,fXiLength(0)
,fOmegaLength(0)
,fMassv0(0)
,fv0DCADaug(0)
,fv0DCAPrimVtx(0)
,fv0Momentum()
,fv0Pt(0)
,fDCABachPrimVtx(0)
,fv0TransRadius(0)
,fv0CPA(0)
,fv0PosToPrimVtx(0)
,fv0NegToPrimVtx(0)
,fv0ToXiPointAngle(0)
,fv0Length(0)
,fDCAv0Xi(0)
,fv0PDG(0)
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
  this->SetEvtNumber(evt->GetRunNumber());
  fCPA=casc->CosPointingAngleXi(PrimVtx[0],PrimVtx[1],PrimVtx[2]);
  fTransRadius=TMath::Sqrt(
      decayPosXi[0]*decayPosXi[0]+decayPosXi[1]*decayPosXi[1]);
  this->SetCharge(casc->ChargeXi());
  this->SetMomentum(casc->MomXiX(),casc->MomXiY(),casc->MomXiZ());
  this->SetPt(casc->MomXiX()*casc->MomXiX()+casc->MomXiY()*casc->MomXiY());
  fRapXi=casc->RapXi();
  fRapOmega=casc->RapOmega();
  this->SetEta(casc->Eta());
  this->SetTheta(casc->Theta());
  this->SetPhi(casc->Phi());
  fAlphaXi=casc->AlphaXi();
  fPtArmXi=casc->PtArmXi();
  fDCAXiPrimVtx=casc->DcaXiToPrimVertex(PrimVtx[0],PrimVtx[1],PrimVtx[2]);
  fXiLength=TMath::Sqrt(TMath::Power((decayPosXi[0]-PrimVtx[0]),2)+
                        TMath::Power((decayPosXi[1]-PrimVtx[1]),2)+
                        TMath::Power((decayPosXi[2]-PrimVtx[2]),2));
  fOmegaLength=fXiLength;
  float XiPDGMass=1.321;
  float OmegaPDGMass=1.672;
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
  if (fIsMC) {
    if (fNegDaug->IsSet())this->SetPhiAtRadius(fNegDaug->GetPhiAtRaidius().at(0));
    if (fPosDaug->IsSet())this->SetPhiAtRadius(fPosDaug->GetPhiAtRaidius().at(0));
    if (fBach->IsSet())this->SetPhiAtRadius(fBach->GetPhiAtRaidius().at(0));
  }
  //v0 business
  if (casc->ChargeXi()<0) {//Xi minus
    fMassv0=casc->MassLambda();
  } else {//Xi plus case
    fMassv0=casc->MassAntiLambda();
  }
  fv0Momentum.SetXYZ(casc->MomV0X(),casc->MomV0Y(),casc->MomV0Z());
  fv0Pt=casc->MomV0X()*casc->MomV0X()+casc->MomV0Y()*casc->MomV0Y();
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
  float LamPDGMass=1.115683;
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
    fDCAXiPrimVtx=0;
    fDCAXiDaug=0;
    fTransRadius=0;
    fRapXi=-99;
    fRapOmega=-99;
    fAlphaXi=99;
    fPtArmXi=99;
    fXiLength=0;
    fOmegaLength=0;
    fMassv0=0;
    fv0DCADaug=0;
    fv0DCAPrimVtx=0;
    fv0Momentum.SetXYZ(0,0,0);
    fv0Pt=0;
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
    fPhiAtRadius.clear();
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
  if (fBach->IsSet()&&fPosDaug->IsSet()&&fNegDaug->IsSet()) {
    //look if the bachelor is from a weak decay and find the label of the
    //mother
    int labelBachMother=-1;
    if (fBach->GetParticleOrigin()==AliFemtoDreamBasePart::kWeak) {
//      std::cout << "Bachelor ID" << fBach->GetIDTracks().at(0) << "\n";
      int labelBach=dynamic_cast<AliAODTrack*>(
            casc->GetDecayVertexXi()->GetDaughter(0))->GetLabel();
      labelBachMother=((AliAODMCParticle*)mcarray->At(labelBach))->GetMother();
//      std::cout << "PDG Bach: "<<((AliAODMCParticle*)mcarray->At(labelBach))->GetPdgCode() << "\n";
//      std::cout << "PDG Mother: "<<((AliAODMCParticle*)mcarray->At(labelBachMother))->GetPdgCode() << "\n";
//      std::cout << "Found a weakling and his mother is: " << labelBachMother <<'\n';
      if (labelBachMother < 0) {
        //This should not happen, just in case somethings very fishy
        this->fIsSet=false;
      } else {
        //look if a v0 exists, and if this v0 stems from a weak decay, and
        //get the label of the mother particle
        int labelv0Mother=-1;
        int PDGDaug[2];
        //order doesn't matter for the MatchToMC
        PDGDaug[0]=TMath::Abs(fPosDaug->GetPDGCode());
        PDGDaug[1]=TMath::Abs(fNegDaug->GetPDGCode());
        int labelv0=casc->MatchToMC(TMath::Abs(fv0PDG),mcarray,2,PDGDaug);
        if (labelv0<0) {
          //both daughters exist as an MC Particle, therfore this can only
          //be a fake v0, and therefore a fake candidate overall
          this->SetParticleOrigin(AliFemtoDreamBasePart::kFake);
        } else {
//          std::cout <<"Potential weakling with ID: " << labelv0 << " found \n";
          AliAODMCParticle* mcv0=(AliAODMCParticle*)mcarray->At(labelv0);
          if (!mcv0) {
            this->fIsSet=false;
          } else {
            if (mcv0->IsSecondaryFromWeakDecay()&&
                !(mcv0->IsSecondaryFromMaterial())) {
              //the v0 candidate has to be from a weak decay itself
              labelv0Mother=mcv0->GetMother();
//              std::cout << "Indeed a weakling his mothers number is " << labelv0Mother << "\n";
              if (labelv0Mother < 0) {
                //Again, this should not happen, just in case somethings
                //very fishy
                this->fIsSet=false;
              } else {
                //now, for a real candidate bachelor and v0 need to have the
                //same mom, NO ADOPTION OR PATCHWORK :(
//                std::cout << labelv0Mother << '\t' << labelBachMother << '\n';
                if (labelv0Mother==labelBachMother) {
                  //MOMMY?
                  AliAODMCParticle* mcPart=(AliAODMCParticle*)mcarray->At(labelv0Mother);
//                  std::cout << "MOOOOOM: " << mcPart->GetPdgCode() << "\n";
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
                } else {
                  //combinatorial background
                  this->SetParticleOrigin(AliFemtoDreamBasePart::kFake);
                }
              }
            } else {
              //if its not from a secondary decay, this cant be a candidate
              this->SetParticleOrigin(AliFemtoDreamBasePart::kFake);
            }
          }
        }
      }
    } else {
      //if the bachlor is not from a weak decay, then it is a ...
      this->SetParticleOrigin(AliFemtoDreamBasePart::kFake);
    }
  } else {
    //this means the bachelor and daughter tracks have no corresponding MC Particle
    this->fIsSet=false;
  }
}
