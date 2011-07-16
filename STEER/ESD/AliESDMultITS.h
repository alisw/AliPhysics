#ifndef AliESDMultITS_H
#define AliESDMultITS_H

#include <TObject.h>
#include <TMath.h>

class AliESDMultITS : public TObject 
{
public:
  AliESDMultITS():
    TObject(),
    fPhi(999),
    fEta(999){
    // Default constructor
  }
  AliESDMultITS(Float_t phi, Float_t eta):
    TObject(),
    fPhi(phi),
    fEta(eta) {
    // Constructor
  }
    
  Float_t GetPhi() const {return fPhi;}
  Float_t GetEta() const {return fEta;}
  Float_t GetTheta() const {return 2*TMath::ATan(TMath::Exp(-fEta));}
  
private:

  Float_t  fPhi; // Phi angle
  Float_t  fEta; // Pseudo-rapidity

  ClassDef(AliESDMultITS,1)
};

#endif
