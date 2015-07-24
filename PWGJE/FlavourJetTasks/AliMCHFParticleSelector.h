#ifndef ALIMCHFPARTICLESELECTOR_H
#define ALIMCHFPARTICLESELECTOR_H

#include "AliEmcalMCTrackSelector.h"

class AliMCHFParticleSelector : public AliEmcalMCTrackSelector {
 public:
  AliMCHFParticleSelector();
  AliMCHFParticleSelector(const char *name);
  virtual ~AliMCHFParticleSelector();

  void SetSpecialPDG(Int_t pdg)                         { fSpecialPDG       = pdg  ; }
  
 protected:
  Bool_t                    AcceptParticle(AliAODMCParticle* part) const;

  Bool_t                    IsSpecialPDGDaughter(AliAODMCParticle* part) const;
  Bool_t                    IsSpecialPDGDaughter(Int_t iPart) const;
 
  Int_t                     fSpecialPDG;           // include particles with this PDG code even if they are not primary particles (and exclude their daughters)
  
 private:
  AliMCHFParticleSelector(const AliMCHFParticleSelector&);            // not implemented
  AliMCHFParticleSelector &operator=(const AliMCHFParticleSelector&); // not implemented

  ClassDef(AliMCHFParticleSelector, 1); // Task to select particle in MC events
};
#endif
