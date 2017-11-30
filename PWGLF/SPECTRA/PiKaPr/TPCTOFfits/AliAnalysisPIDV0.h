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
  void Update(AliAnalysisPIDTrack *PosTrack, AliAnalysisPIDTrack *NegTrack, Double_t *InvMasses, Double_t Radius, Double_t DaughterDCA, Double_t CosinePA, Double_t pT, Double_t Eta, Double_t DCAPV);
  AliAnalysisPIDTrack *GetPosAnalysisTrack() { return fPosAnalysisPIDTrack; };
  AliAnalysisPIDTrack *GetNegAnalysisTrack() { return fNegAnalysisPIDTrack; };
  Double_t GetRadius() { return fRadius; };
  Double_t GetIMK0s() { return fInvMK0s; }; //M K0s
  Double_t GetIML() { return fInvML; }; //M Lambda
  Double_t GetIMAL() { return fInvMAL; };//M Antilambda
  Double_t GetDCAV0Daughters() { return fDCAV0Daughters; };
  Double_t GetV0CosinePA() { return fV0CosinePA; };
  Double_t GetPt() { return fPt; };
  Double_t GetEta() { return fEta; };
  Int_t GetMCPdgCode();
  Double_t GetDCAPV() { return fDCAPV; };
  Double_t CalculateDCAPV();
 protected:  
  AliAnalysisPIDTrack *fPosAnalysisPIDTrack;
  AliAnalysisPIDTrack *fNegAnalysisPIDTrack;
  Double_t fInvMK0s,fInvML,fInvMAL;
  Double_t fRadius;
  Double_t fDCAV0Daughters, fV0CosinePA;
  Double_t fPt;
  //Int_t fMCPdgCode;
  Double_t fEta;
  Double_t fDCAPV;
  ClassDef(AliAnalysisPIDV0, 5);
};
#endif
