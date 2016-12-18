#ifndef ALIV0PARTICLESTRANGE_H
#define ALIV0PARTICLESTRANGE_H

#include "AliConversionPhotonBase.h"
#include "AliKFConversionPhoton.h"
#include "AliAODConversionParticle.h"

class AliV0ParticleStrange : public AliAODConversionParticle, public AliConversionPhotonBase {

  public: 
    
    //Constructors
    AliV0ParticleStrange();    
    AliV0ParticleStrange(AliKFParticle *kfparticle);
    AliV0ParticleStrange(TLorentzVector *vec);

    //Copy Constructor
    AliV0ParticleStrange(const AliV0ParticleStrange & g);           
    //assignment operator
    AliV0ParticleStrange & operator = (const AliV0ParticleStrange & g);

    //Destructor
    virtual ~AliV0ParticleStrange();

    // Overwrite GetLabelFunctions to Make it accessible via AliAODConversionParticle
    virtual Int_t GetLabel(Int_t i) const { return AliConversionPhotonBase::GetTrackLabel(i); }
    virtual Int_t GetLabel1() const { return AliConversionPhotonBase::GetTrackLabelPositive(); }
    virtual Int_t GetLabel2() const { return AliConversionPhotonBase::GetTrackLabelNegative(); }

    virtual Double_t GetPhotonMass() const {return AliAODConversionParticle::M();}
    virtual Double_t GetPhotonPt() const {return AliAODConversionParticle::Pt();}
    virtual Double_t GetPhotonP() const {return AliAODConversionParticle::P();}
    virtual Double_t GetPhotonEta() const {return AliAODConversionParticle::Eta();}
    virtual Double_t GetPhotonTheta() const {return AliAODConversionParticle::Theta();}
    virtual Double_t GetPhotonPhi() const {return AliAODConversionParticle::Phi();}
    virtual Double_t GetPx() const { return AliAODConversionParticle::Px();}
    virtual Double_t GetPy() const { return AliAODConversionParticle::Py();}
    virtual Double_t GetPz() const { return AliAODConversionParticle::Pz();}
    void CalculateDistanceOfClossetApproachToPrimVtx(const AliVVertex* primVertex);
    void SetMassToZero() { SetE(P()); }
    
    void SetInvMassPair(Float_t mass) {fInvMassPair=mass;}
    Float_t GetInvMassPair(){return fInvMassPair;}
    
    
    Float_t GetDCAzToPrimVtx()const {return fDCAzPrimVtx;}
    Float_t GetDCArToPrimVtx()const {return fDCArPrimVtx;}
    
    
  Float_t fDCArPrimVtx;
    Float_t fDCAzPrimVtx;
    Float_t fInvMassPair;
      
	
    ClassDef(AliV0ParticleStrange,5)
};


#endif



