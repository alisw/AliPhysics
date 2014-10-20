#ifndef ALIPICOJET_H
#define ALIPICOJET_H

#include <TObject.h>

class TLorentzVector;

class AliEmcalJet;

class AliPicoJet : public TObject {

 public :

  AliPicoJet();
  AliPicoJet(AliEmcalJet const *pJet, Double_t dLeadingPt);
  AliPicoJet(const AliPicoJet &src);
  AliPicoJet& operator=(const AliPicoJet &src);
  virtual ~AliPicoJet();

  Double_t       Area() const { return fArea;      }
  Double_t  LeadingPt() const { return fLeadingPt; }
  TLorentzVector Kine() const { return fKine;      }

 private :

  TLorentzVector fKine; //

  Double_t fArea;       //
  Double_t fLeadingPt;  //

  ClassDef(AliPicoJet, 3)
};

#endif
