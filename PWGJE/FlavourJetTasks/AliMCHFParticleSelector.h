#ifndef ALIMCHFPARTICLESELECTOR_H
#define ALIMCHFPARTICLESELECTOR_H

#include "AliEmcalMCTrackSelector.h"

class AliMCHFParticleSelector : public AliEmcalMCTrackSelector {
 public:
  AliMCHFParticleSelector();
  AliMCHFParticleSelector(const char *name);
  virtual ~AliMCHFParticleSelector();

  void   SetSpecialPDG(Int_t pdg)          { fSpecialPDG = pdg            ; }

  void   SetRejectQuarkNotFound(Bool_t c)  { fRejectQuarkNotFound = c     ; }
  Bool_t GetRejectQuarkNotFound() const    { return fRejectQuarkNotFound  ; }

  void   SetRejectDfromB(Bool_t c)         { fRejectDfromB = c            ; }
  Bool_t GetRejectDfromB() const           { return fRejectDfromB         ; }

  void   SetKeepOnlyDfromB(Bool_t c)       { fKeepOnlyDfromB = c          ; }
  Bool_t GetKeepOnlyDfromB() const         { return fKeepOnlyDfromB       ; }
  
 protected:
  Bool_t          AcceptParticle(AliAODMCParticle* part) const;

  Bool_t          IsSpecialPDGDaughter(AliAODMCParticle* part) const;
  Bool_t          IsSpecialPDGDaughter(Int_t iPart) const;
 
  Int_t           fSpecialPDG;             //  include particles with this PDG code even if they are not primary particles (and exclude their daughters)
  Bool_t          fRejectQuarkNotFound;    //  reject D mesons for which the original charm or bottom quark could not be found (MC)
  Bool_t          fRejectDfromB;           //  reject D mesons coming from a B meson decay (MC)
  Bool_t          fKeepOnlyDfromB;         //  only accept D mesons coming from a B meson decay (MC)
  
 private:
  AliMCHFParticleSelector(const AliMCHFParticleSelector&);            // not implemented
  AliMCHFParticleSelector &operator=(const AliMCHFParticleSelector&); // not implemented

  ClassDef(AliMCHFParticleSelector, 1); // Task to select particle in MC events
};
#endif
