#ifndef ALIMCPARTICLECONTAINER_H
#define ALIMCPARTICLECONTAINER_H

class AliVEvent;
class AliVParticle;
class AliTLorentzVector;

#include <TArrayC.h>

#include "AliAODMCParticle.h"
#include "AliParticleContainer.h"

typedef AliEmcalIterableContainerT<AliAODMCParticle> AliMCParticleIterableContainer;

/**
 * @class AliMCParticleContainer
 * @brief Container for MC-true particles within the EMCAL framework
 * @ingroup EMCALCOREFW
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 */
class AliMCParticleContainer : public AliParticleContainer {
 public:

  AliMCParticleContainer();
  AliMCParticleContainer(const char *name);
  virtual ~AliMCParticleContainer(){;}

  virtual Bool_t              ApplyMCParticleCuts(const AliAODMCParticle* vp, UInt_t &rejectionReason) const;
  virtual Bool_t              AcceptObject(Int_t i, UInt_t &rejectionReason) const { return AcceptMCParticle(i, rejectionReason);}
  virtual Bool_t              AcceptObject(const TObject* obj, UInt_t &rejectionReason) const { return AcceptMCParticle(dynamic_cast<const AliAODMCParticle*>(obj), rejectionReason);}
  virtual Bool_t              AcceptParticle(Int_t i, UInt_t &rejectionReason) const { return AcceptMCParticle(i, rejectionReason);}
  virtual Bool_t              AcceptParticle(const AliVParticle* vp, UInt_t &rejectionReason) const { return AcceptMCParticle(dynamic_cast<const AliAODMCParticle*>(vp), rejectionReason);}
  virtual Bool_t              AcceptMCParticle(const AliAODMCParticle* vp, UInt_t &rejectionReason) const;
  virtual Bool_t              AcceptMCParticle(Int_t i, UInt_t &rejectionReason) const;
  virtual AliAODMCParticle   *GetMCParticleWithLabel(Int_t lab)         const;
  virtual AliAODMCParticle   *GetAcceptMCParticleWithLabel(Int_t lab)        ;
  virtual AliAODMCParticle   *GetLeadingMCParticle(const char* opt="")        { return static_cast<AliAODMCParticle*>(GetLeadingParticle(opt)); }
  virtual AliAODMCParticle   *GetMCParticle(Int_t i=-1)                 const;
  virtual AliAODMCParticle   *GetAcceptMCParticle(Int_t i=-1)           const;
  virtual AliAODMCParticle   *GetNextAcceptMCParticle()                      ;
  virtual AliAODMCParticle   *GetNextMCParticle()                            ;
  virtual AliVParticle       *GetParticle(Int_t i=-1)                   const { return GetMCParticle(i)           ; }
  virtual AliVParticle       *GetAcceptParticle(Int_t i=-1)             const { return GetAcceptMCParticle(i)     ; }
  virtual AliVParticle       *GetNextAcceptParticle()                         { return GetNextAcceptMCParticle()  ; }
  virtual AliVParticle       *GetNextParticle()                               { return GetNextMCParticle()        ; }

  void                        SetMCFlag(UInt_t m)                             { fMCFlag          = m ; }
  void                        SelectPhysicalPrimaries(Bool_t s)               { if (s) fMCFlag |=  AliAODMCParticle::kPhysicalPrim ;   }

  const char*                 GetTitle() const;

  const AliMCParticleIterableContainer      all() const;
  const AliMCParticleIterableContainer      accepted() const;

  AliMCParticleIterableContainer::iterator  accept_begin()  const { return accepted().begin()   ; }
  AliMCParticleIterableContainer::iterator  accept_end()    const { return accepted().end()     ; }
  AliMCParticleIterableContainer::iterator  accept_rbegin() const { return accepted().rbegin()  ; }
  AliMCParticleIterableContainer::iterator  accept_rend()   const { return accepted().rend()    ; }

  AliMCParticleIterableContainer::iterator  begin()         const { return all().begin()        ; }
  AliMCParticleIterableContainer::iterator  end()           const { return all().end()          ; }
  AliMCParticleIterableContainer::iterator  rbegin()        const { return all().rbegin()       ; }
  AliMCParticleIterableContainer::iterator  rend()          const { return all().rend()         ; }

 protected:
  UInt_t                      fMCFlag;                        ///< select MC particles with flags

 private:
  AliMCParticleContainer(const AliMCParticleContainer& obj); // copy constructor
  AliMCParticleContainer& operator=(const AliMCParticleContainer& other); // assignment

  /// \cond CLASSIMP
  ClassDef(AliMCParticleContainer,1);
  /// \endcond
};

#endif

