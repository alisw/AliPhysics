#ifndef ALIUTILITYDIMUONSOURCE_H
#define ALIUTILITYDIMUONSOURCE_H

//
// MC utility to classify di-muons
//
// Author: Diego Stocco
//

#include "TObject.h"
#include "AliUtilityMuonAncestor.h"

class TString;
class AliMCEvent;
class AliVParticle;


class AliUtilityDimuonSource : public TObject {
public:
  
  AliUtilityDimuonSource();
  ~AliUtilityDimuonSource();

  enum {
    kStopAtLightQuark = 1<<0,
    kStopAtQuark = 1<<1
  };
  
  Int_t GetCommonAncestor ( const AliVParticle* track1, const AliVParticle* track2, const AliMCEvent* mcEvent, UInt_t mask = kStopAtLightQuark ) const;
  TString GetPairType ( Int_t partType1, Int_t partType2, Int_t commonAncestor, const AliMCEvent* mcEvent ) const;
  TString GetPairType ( const AliVParticle* part1, const AliVParticle* part2, const AliMCEvent* mcEvent );

  enum {
    kCharmChainMu,  ///< Mu from charm decay chain
    kBeautyChainMu, ///< Mu from beauty
    kSecondaryMu,   ///< Secondary mu
    kRecoHadron,    ///< Reconstructed hadron
    kOtherMu,       ///< other muon type
    kUnidentified,  ///< Particle that fails matching kine
    kNtrackSources  ///< Total number of track sources
  };

  Int_t GetParticleType ( const AliVParticle* track, const AliMCEvent* mcEvent);
  
private:

  AliUtilityMuonAncestor fUtilityMuonAncestor;
  
  ClassDef(AliUtilityDimuonSource, 0);
};

#endif
