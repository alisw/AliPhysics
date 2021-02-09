#ifndef ALIPICODQ_H
#define ALIPICODQ_H

#include <TObject.h>
#include <TLorentzVector.h>

class AliPicoDQ : public TObject {

 public:

  AliPicoDQ();
  AliPicoDQ(const TLorentzVector &v);
  AliPicoDQ(const AliPicoDQ &a);
  AliPicoDQ &operator=(const AliPicoDQ &a);
  virtual ~AliPicoDQ();
//=============================================================================

  void GetKine(TLorentzVector &v) { v.SetPtEtaPhiM(fPt,fEta,fPhi,fM); }
//=============================================================================

 private:

  Double_t fPt;
  Double_t fEta;
  Double_t fPhi;
  Double_t fM;
//=============================================================================

  ClassDef(AliPicoDQ, 2);
};

#endif

