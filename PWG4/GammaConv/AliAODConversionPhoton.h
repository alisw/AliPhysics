#ifndef ALIAODCONVERSIONPHOTON_H
#define ALIAODCONVERSIONPHOTON_H

#include "AliConversionPhotonBase.h"
#include "AliKFConversionPhoton.h"
#include "AliAODConversionParticle.h"

class AliAODConversionPhoton : public AliAODConversionParticle, public AliConversionPhotonBase {

 public: 

  //Constructors
  AliAODConversionPhoton();    
  AliAODConversionPhoton(AliKFConversionPhoton *kfphoton);

  //Copy Constructor
  AliAODConversionPhoton(const AliAODConversionPhoton & g);           
  //assignment operator
  AliAODConversionPhoton & operator = (const AliAODConversionPhoton & g);

  //Destructor
  virtual ~AliAODConversionPhoton();

  // Overwrite GetLabelFunctions to Make it accessible via AliAODConversionParticle
  virtual Int_t GetLabel(Int_t i) const { return AliConversionPhotonBase::GetTrackLabel(i); };
  virtual Int_t GetLabel1() const { return AliConversionPhotonBase::GetTrackLabelPositive(); };
  virtual Int_t GetLabel2() const { return AliConversionPhotonBase::GetTrackLabelNegative(); };

  virtual Double_t GetPhotonMass() const {return AliAODConversionParticle::M();}
  virtual Double_t GetPhotonPt() const {return AliAODConversionParticle::Pt();}
  virtual Double_t GetPhotonP() const {return AliAODConversionParticle::P();}
  virtual Double_t GetPhotonEta() const {return AliAODConversionParticle::Eta();}

  ClassDef(AliAODConversionPhoton,1)
};


#endif



