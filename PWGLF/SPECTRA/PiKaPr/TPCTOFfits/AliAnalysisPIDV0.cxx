#include "AliAnalysisPIDV0.h"

AliAnalysisPIDV0::AliAnalysisPIDV0():
  TObject(),
  fPosAnalysisPIDTrack(),
  fNegAnalysisPIDTrack(),
  fInvMG(0),
  fInvMK0s(0),
  fInvML(0),
  fInvMAL(0),
  fRadius(0),
  fDCAV0Daughters(0),
  fV0CosinePA(0)
{
};
AliAnalysisPIDV0::AliAnalysisPIDV0(const AliAnalysisPIDV0 &source):
  TObject(source),
  fPosAnalysisPIDTrack(source.fPosAnalysisPIDTrack),
  fNegAnalysisPIDTrack(source.fNegAnalysisPIDTrack),
  fInvMG(source.fInvMG),
  fInvMK0s(source.fInvMK0s),
  fInvML(source.fInvML),
  fInvMAL(source.fInvMAL),
  fRadius(source.fRadius) ,
  fDCAV0Daughters(source.fDCAV0Daughters),
  fV0CosinePA(source.fV0CosinePA)
{
};
AliAnalysisPIDV0 &AliAnalysisPIDV0::operator=(const AliAnalysisPIDV0 &source) {
  if(&source == this) return *this;
  TObject::operator=(source);
  fPosAnalysisPIDTrack = source.fPosAnalysisPIDTrack;
  fNegAnalysisPIDTrack = source.fNegAnalysisPIDTrack;
  fInvMG = source.fInvMG;
  fInvMK0s = source.fInvMK0s;
  fInvML = source.fInvML;
  fInvMAL = source.fInvMAL;
  fRadius = source.fRadius;
  fDCAV0Daughters = source.fDCAV0Daughters;
  fV0CosinePA = source.fV0CosinePA;
  return *this;
};
void AliAnalysisPIDV0::Update(AliAnalysisPIDTrack *PosTrack, AliAnalysisPIDTrack *NegTrack, Double_t *InvMasses, Double_t Radius, Double_t DaughterDCA, Double_t CosinePA) {
  fPosAnalysisPIDTrack = PosTrack;
  fNegAnalysisPIDTrack = NegTrack;
  fInvMG=InvMasses[0];
  fInvMK0s=InvMasses[1];
  fInvML=InvMasses[2];
  fInvMAL=InvMasses[3];
  fRadius = Radius; 
  fDCAV0Daughters = DaughterDCA;
  fV0CosinePA = CosinePA;
};
AliAnalysisPIDV0::~AliAnalysisPIDV0()
{
  delete fPosAnalysisPIDTrack;
  delete fNegAnalysisPIDTrack;
};
