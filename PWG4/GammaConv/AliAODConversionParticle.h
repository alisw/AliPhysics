#ifndef ALIAODCONVERSIONPARTICLE_H
#define ALIAODCONVERSIONPARTICLE_H

#include "AliKFParticle.h"
#include "TLorentzVector.h"


class AliAODConversionParticle : public TLorentzVector {

 public: 

  //Constructors
  AliAODConversionParticle();    
  AliAODConversionParticle(AliKFParticle *kfparticle);

  //Copy Constructor
  AliAODConversionParticle(const AliAODConversionParticle & g);           
  //assignment operator
  AliAODConversionParticle & operator = (const AliAODConversionParticle & g);

  //Destructor
  virtual ~AliAODConversionParticle();

  //Overwrite Phi
  Double_t Phi() const;

  virtual Int_t GetLabel(Int_t i) const = 0;
  virtual Int_t GetLabel1() const { return GetLabel(0); }
  virtual Int_t GetLabel2() const { return GetLabel(1); }

 private:

  ClassDef(AliAODConversionParticle,1)

};

#endif
