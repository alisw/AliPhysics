#ifndef AliTRDsim_H
#define AliTRDsim_H

#include "TObject.h"

class AliTRDsim : public TObject {
private:
  Int_t fNj;            // Number of channel of the histogram
  Float_t fU, fL, fBin;
  Double_t fRo1, fRo2, fOmega1, fOmega2;
  Int_t fIrst;

public:
  AliTRDsim();
  virtual ~AliTRDsim() {}
  virtual void trd_sim();
  virtual void xtr(Double_t gamma, Double_t omega1, Double_t omega2,
		    Float_t *sigmaRad, Int_t &np, Float_t *trEn);
  virtual Float_t fsigmaRad(Int_t ifl, Int_t ig, Float_t o);
  virtual Int_t locate(Double_t *xv, Int_t n, Double_t xval, 
		       Int_t &kl, Double_t &dx);
  virtual Float_t hisran(Float_t *y, Int_t n, Float_t xlo, Float_t xwid);

  
  ClassDef(AliTRDsim,1)  //Base class for all Alice hits
};
#endif
