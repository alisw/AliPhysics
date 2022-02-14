#ifndef ALIANALYSISPIDCASCADEV0__H
#define ALIANALYSISPIDCASCADEV0__H
#include "TObject.h"
#include "AliAnalysisPIDCascadeTrack.h"

class AliAnalysisPIDCascadeTrack;
class AliAnalysisPIDCascadeV0:
public TObject
{
 public:
  AliAnalysisPIDCascadeV0();
  AliAnalysisPIDCascadeV0(const AliAnalysisPIDCascadeV0 &source); //copy const
  AliAnalysisPIDCascadeV0 &operator=(const AliAnalysisPIDCascadeV0 &source); // operator =
  virtual ~AliAnalysisPIDCascadeV0();
  void Update(AliAnalysisPIDCascadeTrack *PosTrack, AliAnalysisPIDCascadeTrack *NegTrack, Double_t *InvMasses, Double_t Radius, Double_t DaughterDCA, Double_t CosinePA, Double_t pT, Double_t Eta, Double_t DCAPV);
  AliAnalysisPIDCascadeTrack *GetPosAnalysisTrack() { return fPosAnalysisPIDCascadeTrack; };
  AliAnalysisPIDCascadeTrack *GetNegAnalysisTrack() { return fNegAnalysisPIDCascadeTrack; };
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
  AliAnalysisPIDCascadeTrack *fPosAnalysisPIDCascadeTrack;
  AliAnalysisPIDCascadeTrack *fNegAnalysisPIDCascadeTrack;
  Double_t fInvMK0s,fInvML,fInvMAL;
  Double_t fRadius;
  Double_t fDCAV0Daughters, fV0CosinePA;
  Double_t fPt;
  //Int_t fMCPdgCode;
  Double_t fEta;
  Double_t fDCAPV;
  ClassDef(AliAnalysisPIDCascadeV0, 2);
};
#endif
