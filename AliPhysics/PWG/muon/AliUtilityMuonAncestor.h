#ifndef ALIUTILITYMUONANCESTOR_H
#define ALIUTILITYMUONANCESTOR_H

/* $Id: AliUtilityMuonAncestor.h 47782 2011-02-24 18:37:31Z martinez $ */ 

//
// MC utility to classify single muons
//
// Author: Diego Stocco
//

#include "TObject.h"

class AliMCEvent;
class AliVParticle;

class AliUtilityMuonAncestor : public TObject {
public:
  
  AliUtilityMuonAncestor();
  ~AliUtilityMuonAncestor();
  AliUtilityMuonAncestor(const AliUtilityMuonAncestor& obj);
  AliUtilityMuonAncestor& operator=(const AliUtilityMuonAncestor& obj);
  
  enum {
    kIsID,
    kIsMuon,
    kIsSecondary,
    kHasLightParent,
    kHasCharmParent,
    kHasBeautyParent,
    kHasQuarkoniumParent,
    kHasTauParent
  };
  
  Bool_t CheckAncestor ( const AliVParticle* track, const AliMCEvent* mcEvent, Int_t ancestorPdg, Bool_t matchAbsPdg = kTRUE );
  
  Int_t GetAncestor ( const AliVParticle* track, const AliMCEvent* mcEvent );
  Int_t GetAncestorPdg ( const AliVParticle* track, const AliMCEvent* mcEvent );
  Long64_t GetMask ( const AliVParticle* track, const AliMCEvent* mcEvent );

  Bool_t IsBeautyMu ( const AliVParticle* track, const AliMCEvent* mcEvent );
  Bool_t IsBeautyChainMu ( const AliVParticle* track, const AliMCEvent* mcEvent );
  Bool_t IsBJpsiMu ( const AliVParticle* track, const AliMCEvent* mcEvent );
  Bool_t IsCharmMu ( const AliVParticle* track, const AliMCEvent* mcEvent );
  Bool_t IsCharmChainMu ( const AliVParticle* track, const AliMCEvent* mcEvent );
  Bool_t IsDecayMu ( const AliVParticle* track, const AliMCEvent* mcEvent );
  Bool_t IsHadron ( const AliVParticle* track, const AliMCEvent* mcEvent );
  Bool_t IsMuon ( const AliVParticle* track, const AliMCEvent* mcEvent );
  Bool_t IsQuarkoniumMu ( const AliVParticle* track, const AliMCEvent* mcEvent );
  Bool_t IsSecondaryMu ( const AliVParticle* track, const AliMCEvent* mcEvent );
  Bool_t IsUnidentified ( const AliVParticle* track, const AliMCEvent* mcEvent );
  Bool_t IsWBosonMu ( const AliVParticle* track, const AliMCEvent* mcEvent );
  Bool_t IsZBosonMu ( const AliVParticle* track, const AliMCEvent* mcEvent );

  
private:
  Bool_t BuildAncestor ( const AliVParticle* track, const AliMCEvent* mcEvent );
  
  
  
  Double_t fPx; ///< Particle px
  Double_t fPy; ///< Particle py
  Double_t fPz; ///< Particle pz
  Long64_t fMask; ///< Mask
  Int_t fAncestor; ///< Ancestor position in stack
  
  ClassDef(AliUtilityMuonAncestor, 0);
};

#endif
