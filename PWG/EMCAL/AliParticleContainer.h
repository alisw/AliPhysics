#ifndef AliParticleContainer_H
#define AliParticleContainer_H

// $Id$

class AliVEvent;
class AliVParticle;

#include "AliEmcalContainer.h"

class AliParticleContainer : public AliEmcalContainer {
 public:
  AliParticleContainer();
  AliParticleContainer(const char *name); 
  virtual ~AliParticleContainer(){;}

  Bool_t                      AcceptParticle(AliVParticle         *vp)  const;
  Double_t                    GetParticlePtCut()                        const   { return fParticlePtCut; }
  Double_t                    GetParticleEtaMin()                       const   { return fParticleMinEta; }
  Double_t                    GetParticleEtaMax()                       const   { return fParticleMaxEta; }
  Double_t                    GetParticlePhiMin()                       const   { return fParticleMinPhi; }
  Double_t                    GetParticlePhiMax()                       const   { return fParticleMaxPhi; }
  AliVParticle               *GetLeadingParticle(const char* opt="")         ;
  AliVParticle               *GetParticle(Int_t i)                      const;
  AliVParticle               *GetAcceptParticle(Int_t i)                const;
  AliVParticle               *GetParticleWithLabel(Int_t lab)           const;
  AliVParticle               *GetAcceptParticleWithLabel(Int_t lab)     const;
  AliVParticle               *GetNextAcceptParticle(Int_t i=-1)              ;
  AliVParticle               *GetNextParticle(Int_t i=-1)                    ;
  void                        GetMomentum(TLorentzVector &mom, Int_t i) const;
  Int_t                       GetNParticles()                           const   {return GetNEntries();}
  Int_t                       GetNAcceptedParticles()                   ;
  void                        SetClassName(const char *clname);
  void                        SetMCTrackBitMap(UInt_t m)                        { fMCTrackBitMap   = m ; }
  void                        SetMinMCLabel(Int_t s)                            { fMinMCLabel      = s ; }
  void                        SetParticlePtCut(Double_t cut)                    { fParticlePtCut = cut ; }
  void                        SetParticleEtaLimits(Double_t min, Double_t max)  { fParticleMaxEta = max ; fParticleMinEta = min ; }
  void                        SetParticlePhiLimits(Double_t min, Double_t max)  { fParticleMaxPhi = max ; fParticleMinPhi = min ; }
  void                        SetTrackBitMap(UInt_t m)                          { fTrackBitMap     = m ; }

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

