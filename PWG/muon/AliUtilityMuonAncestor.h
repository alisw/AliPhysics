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
    kHasQuarkoniumParent
  };
  
  Bool_t CheckAncestor ( AliVParticle* track, const AliMCEvent* mcEvent, Int_t ancestorPdg, Bool_t matchAbsPdg = kTRUE );
  
  Int_t GetAncestor ( AliVParticle* track, const AliMCEvent* mcEvent );
  Int_t GetAncestorPdg ( AliVParticle* track, const AliMCEvent* mcEvent );
  Long64_t GetMask ( AliVParticle* track, const AliMCEvent* mcEvent );

  Bool_t IsBeautyMu ( AliVParticle* track, const AliMCEvent* mcEvent );
  Bool_t IsBJpsiMu ( AliVParticle* track, const AliMCEvent* mcEvent );
  Bool_t IsCharmMu ( AliVParticle* track, const AliMCEvent* mcEvent );
  Bool_t IsDecayMu ( AliVParticle* track, const AliMCEvent* mcEvent );
  Bool_t IsHadron ( AliVParticle* track, const AliMCEvent* mcEvent );
  Bool_t IsMuon ( AliVParticle* track, const AliMCEvent* mcEvent );
  Bool_t IsQuarkoniumMu ( AliVParticle* track, const AliMCEvent* mcEvent );
  Bool_t IsSecondaryMu ( AliVParticle* track, const AliMCEvent* mcEvent );
  Bool_t IsUnidentified ( AliVParticle* track, const AliMCEvent* mcEvent );
  Bool_t IsWBosonMu ( AliVParticle* track, const AliMCEvent* mcEvent );
  Bool_t IsZBosonMu ( AliVParticle* track, const AliMCEvent* mcEvent );

  
private:
  Bool_t BuildAncestor ( AliVParticle* track, const AliMCEvent* mcEvent );
  
  
  
  Double_t fPx; ///< Particle px
  Double_t fPy; ///< Particle py
  Double_t fPz; ///< Particle pz
  Long64_t fMask; ///< Mask
  Int_t fAncestor; ///< Ancestor position in stack
  
  ClassDef(AliUtilityMuonAncestor, 0);
};

#endif
