#ifndef ALIPICOHEADERJET_H
#define ALIPICOHEADERJET_H

#include <AliPicoHeaderV0.h>

class AliPicoHeaderJet : public AliPicoHeaderV0 {

 public :

  AliPicoHeaderJet(const TString s="");
  AliPicoHeaderJet(const AliPicoHeaderJet &src);
  AliPicoHeaderJet& operator=(const AliPicoHeaderJet &src);
  virtual ~AliPicoHeaderJet();
//=============================================================================

  Double_t BackgroundRho() const { return fBkgRho; }
  void BackgroundRho(const Double_t d) { fBkgRho = d; }

  virtual void Reset();
//=============================================================================

 private :

  Double_t fBkgRho; //

  ClassDef(AliPicoHeaderJet, 3);
};

#endif
