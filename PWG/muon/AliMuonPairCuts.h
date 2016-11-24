#ifndef ALIMUONPAIRCUTS_H
#define ALIMUONPAIRCUTS_H

#include "AliAnalysisCuts.h"
#include "AliMuonTrackCuts.h"

class TList;
class TVector3;
class TArrayI;
class AliVEventHandler;

class AliMuonPairCuts : public AliAnalysisCuts
{
 public:
  
  enum {
    kBothMuEta = BIT(0),
    kBothMuThetaAbs = BIT(1),
    kBothMuPdca = BIT(2),
    kBothMuTrackChiSquare = BIT(3),
    kBothMuMatchApt = BIT(4),
    kBothMuMatchLpt = BIT(5),
    kBothMuMatchHpt = BIT(6),
    kOneMuMatchApt = BIT(7),
    kOneMuMatchLpt = BIT(8),
    kOneMuMatchHpt = BIT(9),
    kDimuUnlikeSign = BIT(30),
    kDimuRapidity = BIT(31)
  };

  AliMuonPairCuts();
  AliMuonPairCuts(const char* name, const char* title);
  AliMuonPairCuts(const char* name, const char* title, const AliMuonTrackCuts& trackCuts);
  AliMuonPairCuts(const AliMuonPairCuts& obj);
  AliMuonPairCuts& operator=(const AliMuonPairCuts& obj);

  virtual ~AliMuonPairCuts();

  virtual UInt_t GetSelectionMask ( const TObject* obj );
  virtual Bool_t IsSelected ( TObject* /*obj*/ );
  virtual Bool_t IsSelected ( TList* list );
  virtual Bool_t IsSelected ( TObject* track1, TObject* track2 );

  UInt_t GetSelectionMask(const TObject* track1, const TObject* track2);
  
  void SetDefaultFilterMask();

  Bool_t SetRun ( const AliVEventHandler* eventHandler );
  void SetIsMC ( Bool_t isMC = kTRUE );

  void Print ( Option_t* option = "" ) const;
  
  Double_t MuonMass2() const;
  
  /// Apply also sharp pt cut when matching with trigger
  void ApplySharpPtCutInMatching ( Bool_t sharpPtCut = kTRUE ) { fMuonTrackCuts.ApplySharpPtCutInMatching(sharpPtCut); }
  /// Get flag to apply the sharp pt cut when matching with trigger
  Bool_t IsApplySharpPtCutInMatching () const { return fMuonTrackCuts.IsApplySharpPtCutInMatching(); }
  
  Bool_t TrackPtCutMatchTrigClass ( const AliVParticle* track1, const AliVParticle* track2, const TArrayI ptCutFromClass ) const;
  
  /// Return the MuonTrackCuts (normally the standard user do not need to modify this data member)
  AliMuonTrackCuts& GetMuonTrackCuts () { return fMuonTrackCuts; }

 private:
  
  AliMuonTrackCuts fMuonTrackCuts; ///< Muon track cuts

  ClassDef(AliMuonPairCuts, 1); // Class for muon pair filters

};
 
#endif

