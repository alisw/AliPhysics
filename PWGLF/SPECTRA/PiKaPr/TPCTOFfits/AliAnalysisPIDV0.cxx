#include "AliAnalysisPIDV0.h"

AliAnalysisPIDV0::AliAnalysisPIDV0():
  TObject(),
  fPosAnalysisPIDTrack(),
  fNegAnalysisPIDTrack(),
  fInvMK0s(0),
  fInvML(0),
  fInvMAL(0),
  fRadius(0),
  fDCAV0Daughters(0),
  fV0CosinePA(0),
  fPt(0),
  fEta(0)
{
};
AliAnalysisPIDV0::AliAnalysisPIDV0(const AliAnalysisPIDV0 &source):
  TObject(source),
  fPosAnalysisPIDTrack(source.fPosAnalysisPIDTrack),
  fNegAnalysisPIDTrack(source.fNegAnalysisPIDTrack),
  fInvMK0s(source.fInvMK0s),
  fInvML(source.fInvML),
  fInvMAL(source.fInvMAL),
  fRadius(source.fRadius) ,
  fDCAV0Daughters(source.fDCAV0Daughters),
  fV0CosinePA(source.fV0CosinePA),
  fPt(source.fPt),
  fEta(source.fEta)
{
};
AliAnalysisPIDV0 &AliAnalysisPIDV0::operator=(const AliAnalysisPIDV0 &source) {
  if(&source == this) return *this;
  TObject::operator=(source);
  fPosAnalysisPIDTrack = source.fPosAnalysisPIDTrack;
  fNegAnalysisPIDTrack = source.fNegAnalysisPIDTrack;
  fInvMK0s = source.fInvMK0s;
  fInvML = source.fInvML;
  fInvMAL = source.fInvMAL;
  fRadius = source.fRadius;
  fDCAV0Daughters = source.fDCAV0Daughters;
  fV0CosinePA = source.fV0CosinePA;
  fPt = source.fPt;
  fEta = source.fEta;
  return *this;
};
void AliAnalysisPIDV0::Update(AliAnalysisPIDTrack *PosTrack, AliAnalysisPIDTrack *NegTrack, Double_t *InvMasses, Double_t Radius, Double_t DaughterDCA, Double_t CosinePA, Double_t pT, Double_t Eta) {
  fPosAnalysisPIDTrack = PosTrack;
  fNegAnalysisPIDTrack = NegTrack;
  fInvMK0s=InvMasses[0];
  fInvML=InvMasses[1];
  fInvMAL=InvMasses[2];
  fRadius = Radius; 
  fDCAV0Daughters = DaughterDCA;
  fV0CosinePA = CosinePA;
  fPt = pT;
  fEta = Eta;
};
AliAnalysisPIDV0::~AliAnalysisPIDV0()
{
  delete fPosAnalysisPIDTrack;
  delete fNegAnalysisPIDTrack;
};
Int_t AliAnalysisPIDV0::GetMCPdgCode() {
  if(!fPosAnalysisPIDTrack || !fNegAnalysisPIDTrack) return 0;
  if(fPosAnalysisPIDTrack->GetMCMotherPdgCode()==0) return 0;
  if(fPosAnalysisPIDTrack->GetMCMotherPdgCode()==fNegAnalysisPIDTrack->GetMCMotherPdgCode())
    if(fPosAnalysisPIDTrack->GetMCMotherLabel()==fNegAnalysisPIDTrack->GetMCMotherLabel())
      return fPosAnalysisPIDTrack->GetMCMotherPdgCode();
  return 0;
};
