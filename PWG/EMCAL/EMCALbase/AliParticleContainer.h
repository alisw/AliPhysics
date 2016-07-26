#ifndef ALIPARTICLECONTAINER_H
#define ALIPARTICLECONTAINER_H

#include <AliVParticle.h>

class AliVEvent;
class AliTLorentzVector;

#include "AliEmcalContainer.h"

class AliParticleContainer : public AliEmcalContainer {
 public:

  AliParticleContainer();
  AliParticleContainer(const char *name);
  virtual ~AliParticleContainer(){;}

  virtual Bool_t              ApplyParticleCuts(const AliVParticle* vp);
  virtual Bool_t              ApplyKinematicCuts(const AliTLorentzVector& mom);
  virtual Bool_t              AcceptObject(Int_t i)              { return AcceptParticle(i);}
  virtual Bool_t              AcceptObject(const TObject* obj)   { return AcceptParticle(dynamic_cast<const AliVParticle*>(obj));}
  virtual Bool_t              AcceptParticle(const AliVParticle* vp)         ;
  virtual Bool_t              AcceptParticle(Int_t i)                        ;
  Double_t                    GetParticlePtCut()                        const   { return GetMinPt()     ; }
  Double_t                    GetParticleEtaMin()                       const   { return GetMinEta()    ; }
  Double_t                    GetParticleEtaMax()                       const   { return GetMaxEta()    ; }
  Double_t                    GetParticlePhiMin()                       const   { return GetMinPhi()    ; }
  Double_t                    GetParticlePhiMax()                       const   { return GetMaxPhi()    ; }
  void                        SetParticlePtCut(Double_t cut)                    { SetMinPt(cut)         ; }
  void                        SetParticleEtaLimits(Double_t min, Double_t max)  { SetEtaLimits(min, max); }
  void                        SetParticlePhiLimits(Double_t min, Double_t max)  { SetPhiLimits(min, max); }
  virtual AliVParticle       *GetLeadingParticle(const char* opt="")         ;
  virtual AliVParticle       *GetParticle(Int_t i=-1)                   const;
  virtual AliVParticle       *GetAcceptParticle(Int_t i=-1)                  ;
  virtual AliVParticle       *GetNextAcceptParticle()                        ;
  virtual AliVParticle       *GetNextParticle()                              ;
  virtual Bool_t              GetMomentum(TLorentzVector &mom, const AliVParticle* part, Double_t mass);
  virtual Bool_t              GetMomentum(TLorentzVector &mom, const AliVParticle* part);
  virtual Bool_t              GetMomentum(TLorentzVector &mom, Int_t i);
  virtual Bool_t              GetAcceptMomentum(TLorentzVector &mom, Int_t i);
  virtual Bool_t              GetNextMomentum(TLorentzVector &mom);
  virtual Bool_t              GetNextAcceptMomentum(TLorentzVector &mom);
  Int_t                       GetNParticles()                           const   {return GetNEntries();}
  Int_t                       GetNAcceptedParticles()                   ;

  void                        SetClassName(const char *clname);
  void                        SetMinDistanceTPCSectorEdge(Double_t min)         { fMinDistanceTPCSectorEdge = min; }
  void                        SetCharge(Short_t c)                              { fCharge = c         ; }
  void                        SelectHIJING(Bool_t s)                            { if (s) fGeneratorIndex = 0; else fGeneratorIndex = -1; }
  void                        SetGeneratorIndex(Short_t i)                      { fGeneratorIndex = i  ; }

  const char*                 GetTitle() const;

 protected:

  Double_t                    fMinDistanceTPCSectorEdge;      // require minimum distance to edge of TPC sector edge
  Short_t                     fCharge;                        // select particles with charge=fCharge
  Short_t                     fGeneratorIndex;                // select MC particles with generator index (default = -1 = switch off selection)

 private:
  AliParticleContainer(const AliParticleContainer& obj); // copy constructor
  AliParticleContainer& operator=(const AliParticleContainer& other); // assignment

  ClassDef(AliParticleContainer,9);

};

#endif

