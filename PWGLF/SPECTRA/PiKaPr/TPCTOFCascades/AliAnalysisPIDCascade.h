#ifndef ALIANALYSISPIDCascade__H
#define ALIANALYSISPIDCascade__H
#include "TObject.h"
#include "AliAnalysisPIDTrack.h"
#include "AliAnalysisPIDV0.h"

class AliAnalysisPIDTrack;
class AliAnalysisPIDV0;
class AliAnalysisPIDCascade:
public TObject
{
 public:
  AliAnalysisPIDCascade();
  AliAnalysisPIDCascade(const AliAnalysisPIDCascade &source); //copy const
  AliAnalysisPIDCascade &operator=(const AliAnalysisPIDCascade &source); // operator =
  virtual ~AliAnalysisPIDCascade();
  void Update(AliAnalysisPIDV0* V0, AliAnalysisPIDTrack* BachTrack, Double_t *InvMasses, Double_t CascRadius, Double_t CascDCA,  Double_t CascCosinePA, Double_t PtCasc, Double_t EtaCasc, Double_t CascDCAPV, Int_t Charge);

  AliAnalysisPIDV0    *GetV0()  { return fV0AnalysisPIDV0;};
  AliAnalysisPIDTrack *GetBachAnalysisTrack() { return fBachAnalysisPIDTrack; };
  Double_t GetIMXi()         { return fInvMXi; }; //M Xi
  Double_t GetIMO()          { return fInvMO; }; //M Omega
  Double_t GetCascRadius()   { return fCascRadius; };
  Double_t GetCascDCA()   { return fCascDCA; } ;
  Double_t GetCascCosinePA() { return fCascCosinePA; };
  Double_t GetPtCasc() { return fPtCasc; };
  Double_t GetEtaCasc() { return fEtaCasc; };
  Double_t GetCascDCAPV()    { return fCascDCAPV; };
  Double_t CalculateCascDCAPV();
  Int_t GetCharge() { return fCharge;};

 protected:  
  AliAnalysisPIDV0 * fV0AnalysisPIDV0;
  AliAnalysisPIDTrack *fBachAnalysisPIDTrack;
  Double_t fInvMXi,fInvMO;
  Double_t fCascRadius;
  Double_t fCascDCA, fCascCosinePA;
  Double_t fPtCasc;
  Double_t fEtaCasc;
  Double_t fCascDCAPV;
  Int_t fCharge;
  
  ClassDef(AliAnalysisPIDCascade, 2);
};
#endif
