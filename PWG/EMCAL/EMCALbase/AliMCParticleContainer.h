#ifndef ALIMCPARTICLECONTAINER_H
#define ALIMCPARTICLECONTAINER_H

class AliVEvent;
class AliVParticle;
class AliTLorentzVector;

#include <TArrayC.h>

#include "AliAODMCParticle.h"
#include "AliParticleContainer.h"

class AliMCParticleContainer : public AliParticleContainer {
 public:

  AliMCParticleContainer();
  AliMCParticleContainer(const char *name);
  virtual ~AliMCParticleContainer(){;}

  virtual Bool_t              ApplyMCParticleCuts(const AliAODMCParticle* vp);
  virtual Bool_t              AcceptObject(Int_t i)                           { return AcceptMCParticle(i);}
  virtual Bool_t              AcceptObject(const TObject* obj)                { return AcceptMCParticle(dynamic_cast<const AliAODMCParticle*>(obj));}
  virtual Bool_t              AcceptParticle(Int_t i)                         { return AcceptMCParticle(i);}
  virtual Bool_t              AcceptParticle(const AliVParticle* vp)          { return AcceptMCParticle(dynamic_cast<const AliAODMCParticle*>(vp));}
  virtual Bool_t              AcceptMCParticle(const AliAODMCParticle* vp)   ;
  virtual Bool_t              AcceptMCParticle(Int_t i)                      ;
  virtual AliAODMCParticle   *GetMCParticleWithLabel(Int_t lab)         const;
  virtual AliAODMCParticle   *GetAcceptMCParticleWithLabel(Int_t lab)        ;
  virtual AliAODMCParticle   *GetLeadingMCParticle(const char* opt="")        { return static_cast<AliAODMCParticle*>(GetLeadingParticle(opt)); }
  virtual AliAODMCParticle   *GetMCParticle(Int_t i=-1)                 const;
  virtual AliAODMCParticle   *GetAcceptMCParticle(Int_t i=-1)                ;
  virtual AliAODMCParticle   *GetNextAcceptMCParticle()                      ;
  virtual AliAODMCParticle   *GetNextMCParticle()                            ;
  virtual AliVParticle       *GetParticle(Int_t i=-1)                   const { return GetMCParticle(i)           ; }
  virtual AliVParticle       *GetAcceptParticle(Int_t i=-1)                   { return GetAcceptMCParticle(i)     ; }
  virtual AliVParticle       *GetNextAcceptParticle()                         { return GetNextAcceptMCParticle()  ; }
  virtual AliVParticle       *GetNextParticle()                               { return GetNextMCParticle()        ; }
  virtual Bool_t              GetMomentum(TLorentzVector &mom, const AliAODMCParticle* part, Double_t mass);
  virtual Bool_t              GetMomentum(TLorentzVector &mom, const AliAODMCParticle* part);
  virtual Bool_t              GetMomentum(TLorentzVector &mom, Int_t i);
  virtual Bool_t              GetAcceptMomentum(TLorentzVector &mom, Int_t i);
  virtual Bool_t              GetNextMomentum(TLorentzVector &mom);
  virtual Bool_t              GetNextAcceptMomentum(TLorentzVector &mom);

  void                        SetClassName(const char *clname);
  void                        SetMCFlag(UInt_t m)                               { fMCFlag          = m ; }
  void                        SelectPhysicalPrimaries(Bool_t s)                 { if (s) fMCFlag |=  AliAODMCParticle::kPhysicalPrim ;   }

  const char*                 GetTitle() const;

 protected:
  UInt_t                      fMCFlag;                        // select MC particles with flags

 private:
  AliMCParticleContainer(const AliMCParticleContainer& obj); // copy constructor
  AliMCParticleContainer& operator=(const AliMCParticleContainer& other); // assignment

  ClassDef(AliMCParticleContainer,1);

};

#endif

