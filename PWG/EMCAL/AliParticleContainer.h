#ifndef AliParticleContainer_H
#define AliParticleContainer_H

//
// container with name, TClonesArray and cuts for particles
//

class AliVEvent;
class AliVParticle;

#include "AliEmcalContainer.h"

class AliParticleContainer : public AliEmcalContainer {
 public:
  AliParticleContainer();
  AliParticleContainer(const char *name); 
  virtual ~AliParticleContainer(){;}

  void  SetParticleArray(AliVEvent *event);

  void  SetParticlePtCut(Double_t cut)                   { fParticlePtCut = cut ; }
  void  SetParticleEtaLimits(Double_t min, Double_t max) { fParticleMaxEta = max ; fParticleMinEta = min ; }
  void  SetParticlePhiLimits(Double_t min, Double_t max) { fParticleMaxPhi = max ; fParticleMinPhi = min ; }
  void  SetTrackBitMap(UInt_t m)                          { fTrackBitMap     = m ; }
  void  SetMCTrackBitMap(UInt_t m)                        { fMCTrackBitMap   = m ; }
  void  SetMinMCLabel(Int_t s)                            { fMinMCLabel      = s ; }

  AliVParticle               *GetLeadingParticle(const char* opt="")   const;
  AliVParticle               *GetParticle(Int_t i)                     const;
  AliVParticle               *GetAcceptParticle(Int_t i)               const;
  AliVParticle               *GetNextAcceptParticle(Int_t i=-1)        const;
  void                        GetMomentum(TLorentzVector &mom, Int_t i) const;
  Bool_t                      AcceptParticle(AliVParticle         *vp) const;
  Int_t                       GetNParticles()                          const   {return GetNEntries();}

 protected:
  Double_t                    fParticlePtCut;                 // cut on particle pt
  Double_t                    fParticleMinEta;                // cut on particle eta
  Double_t                    fParticleMaxEta;                // cut on particle eta
  Double_t                    fParticleMinPhi;                // cut on particle phi
  Double_t                    fParticleMaxPhi;                // cut on particle phi
  UInt_t                      fTrackBitMap;                   // bit map of accepted tracks (non MC)
  UInt_t                      fMCTrackBitMap;                 // bit map of accepted MC tracks
  Int_t                       fMinMCLabel;                    // minimum MC label value for the tracks/clusters being considered MC particles

 private:
  AliParticleContainer(const AliParticleContainer& obj); // copy constructor
  AliParticleContainer& operator=(const AliParticleContainer& other); // assignment

  ClassDef(AliParticleContainer,1);

};

#endif

