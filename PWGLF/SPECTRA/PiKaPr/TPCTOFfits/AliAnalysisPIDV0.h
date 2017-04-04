#ifndef ALIANALYSISPIDV0__H
#define ALIANALYSISPIDV0__H
#include "TObject.h"
#include "AliAnalysisPIDTrack.h"

class AliAnalysisPIDTrack;
class AliAnalysisPIDV0:
public TObject
{
 public:
  AliAnalysisPIDV0();
  AliAnalysisPIDV0(const AliAnalysisPIDV0 &source); //copy const
  AliAnalysisPIDV0 &operator=(const AliAnalysisPIDV0 &source); // operator =
  virtual ~AliAnalysisPIDV0();
  void Update(AliAnalysisPIDTrack *PosTrack, AliAnalysisPIDTrack *NegTrack, Double_t *InvMasses, Double_t Radius, Double_t DaughterDCA, Double_t CosinePA);
  AliAnalysisPIDTrack *GetPosAnalysisTrack() { return fPosAnalysisPIDTrack; };
  AliAnalysisPIDTrack *GetNegAnalysisTrack() { return fNegAnalysisPIDTrack; };
  Double_t GetRadius() { return fRadius; };
  Double_t GetIMGamma() { return fInvMG; }; //M Gamma
  Double_t GetIMK0s() { return fInvMK0s; }; //M K0s
  Double_t GetIML() { return fInvML; }; //M Lambda
  Double_t GetIMAL() { return fInvMAL; };//M Antilambda
  Double_t GetDCAV0Daughters() { return fDCAV0Daughters; };
  Double_t GetV0CosinePA() { return fV0CosinePA; };
 protected:  
  Double_t fInvMG,fInvMK0s,fInvML,fInvMAL;
  Double_t fRadius;
  Double_t fDCAV0Daughters, fV0CosinePA;
  AliAnalysisPIDTrack *fPosAnalysisPIDTrack;
  AliAnalysisPIDTrack *fNegAnalysisPIDTrack;
  ClassDef(AliAnalysisPIDV0, 2);
};
#endif
