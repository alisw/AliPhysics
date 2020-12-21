#include "AliAnalysisPIDCascadeV0.h"

AliAnalysisPIDCascadeV0::AliAnalysisPIDCascadeV0():
  TObject(),
  fPosAnalysisPIDCascadeTrack(),
  fNegAnalysisPIDCascadeTrack(),
  fInvMK0s(0),
  fInvML(0),
  fInvMAL(0),
  fRadius(0),
  fDCAV0Daughters(0),
  fV0CosinePA(0),
  fPt(0),
  fEta(0),
  fDCAPV(-999)
{
};
AliAnalysisPIDCascadeV0::AliAnalysisPIDCascadeV0(const AliAnalysisPIDCascadeV0 &source):
  TObject(source),
  fPosAnalysisPIDCascadeTrack(source.fPosAnalysisPIDCascadeTrack),
  fNegAnalysisPIDCascadeTrack(source.fNegAnalysisPIDCascadeTrack),
  fInvMK0s(source.fInvMK0s),
  fInvML(source.fInvML),
  fInvMAL(source.fInvMAL),
  fRadius(source.fRadius) ,
  fDCAV0Daughters(source.fDCAV0Daughters),
  fV0CosinePA(source.fV0CosinePA),
  fPt(source.fPt),
  fEta(source.fEta),
  fDCAPV(source.fDCAPV)
{
};
AliAnalysisPIDCascadeV0 &AliAnalysisPIDCascadeV0::operator=(const AliAnalysisPIDCascadeV0 &source) {
  if(&source == this) return *this;
  TObject::operator=(source);
  fPosAnalysisPIDCascadeTrack = source.fPosAnalysisPIDCascadeTrack;
  fNegAnalysisPIDCascadeTrack = source.fNegAnalysisPIDCascadeTrack;
  fInvMK0s = source.fInvMK0s;
  fInvML = source.fInvML;
  fInvMAL = source.fInvMAL;
  fRadius = source.fRadius;
  fDCAV0Daughters = source.fDCAV0Daughters;
  fV0CosinePA = source.fV0CosinePA;
  fPt = source.fPt;
  fEta = source.fEta;
  fDCAPV = source.fDCAPV;
  return *this;
};
void AliAnalysisPIDCascadeV0::Update(AliAnalysisPIDCascadeTrack *PosTrack, AliAnalysisPIDCascadeTrack *NegTrack, Double_t *InvMasses, Double_t Radius, Double_t DaughterDCA, Double_t CosinePA, Double_t pT, Double_t Eta, Double_t DCAPV) {
  fPosAnalysisPIDCascadeTrack = PosTrack;
  fNegAnalysisPIDCascadeTrack = NegTrack;
  fInvMK0s=InvMasses[0];
  fInvML=InvMasses[1];
  fInvMAL=InvMasses[2];
  fRadius = Radius; 
  fDCAV0Daughters = DaughterDCA;
  fV0CosinePA = CosinePA;
  fPt = pT;
  fEta = Eta;
  fDCAPV = DCAPV;
};
AliAnalysisPIDCascadeV0::~AliAnalysisPIDCascadeV0()
{
  delete fPosAnalysisPIDCascadeTrack;
  delete fNegAnalysisPIDCascadeTrack;
};
Int_t AliAnalysisPIDCascadeV0::GetMCPdgCode() {
  if(!fPosAnalysisPIDCascadeTrack || !fNegAnalysisPIDCascadeTrack) return 0;
  if(fPosAnalysisPIDCascadeTrack->GetMCMotherPdgCode()==0) return 0;
  if(fPosAnalysisPIDCascadeTrack->GetMCMotherPdgCode()==fNegAnalysisPIDCascadeTrack->GetMCMotherPdgCode())
    if(fPosAnalysisPIDCascadeTrack->GetMCMotherLabel()==fNegAnalysisPIDCascadeTrack->GetMCMotherLabel())
      return fPosAnalysisPIDCascadeTrack->GetMCMotherPdgCode();
  return 0;
};
Double_t AliAnalysisPIDCascadeV0::CalculateDCAPV() {
  return fRadius * ( TMath::Sqrt(1. - fV0CosinePA * fV0CosinePA));
};
