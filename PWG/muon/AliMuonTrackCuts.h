#ifndef ALIMUONTRACKCUTS_H
#define ALIMUONTRACKCUTS_H

#include "AliAnalysisCuts.h"
#include "AliOADBMuonTrackCutsParam.h"

class TList;
class TArrayI;
class AliVParticle;
class AliInputEventHandler;


class AliMuonTrackCuts : public AliAnalysisCuts
{
 public:
  
  enum {
    kMuEta = BIT(0),
    kMuThetaAbs = BIT(1),
    kMuPdca = BIT(2),
    kMuMatchApt = BIT(3),
    kMuMatchLpt = BIT(4),
    kMuMatchHpt = BIT(5),
    kMuTrackChiSquare = BIT(6)
  };
  
  AliMuonTrackCuts();
  AliMuonTrackCuts(const char* name, const char* title);
  AliMuonTrackCuts(const AliMuonTrackCuts& obj);
  AliMuonTrackCuts& operator=(const AliMuonTrackCuts& obj);

  virtual ~AliMuonTrackCuts();

  virtual UInt_t GetSelectionMask ( const TObject* obj );
  virtual Bool_t IsSelected ( TObject* obj );
  virtual Bool_t IsSelected ( TList* /*list */ );
  
  void SetDefaultFilterMask();
  void SetPassNumber ( Int_t passNumber ) { fPassNumber = passNumber; }
  void SetIsMC ( Bool_t isMC = kTRUE ) { fIsMC = isMC; }
  void SetAllowDefaultParams ( Bool_t allowDefaultParams = kTRUE, Int_t passNumber = -1 ) { fAllowDefaultParams = allowDefaultParams; fPassNumber = passNumber; }
  void SetCustomParamFromRun ( Int_t runNumber = -1, Int_t passNumber = -1 );
  
  /// Get pass number
  Int_t GetPassNumber () const { return fPassNumber; }

  Bool_t SetRun ( const AliInputEventHandler* eventHandler );
  
  void Print ( Option_t* option = "" ) const;
  
  Bool_t TrackPtCutMatchTrigClass ( const AliVParticle* track, const TArrayI ptCutFromClass) const;

  TVector3 GetCorrectedDCA ( const AliVParticle* track ) const;
  Double_t GetAverageMomentum ( const AliVParticle* track ) const;
  Bool_t IsThetaAbs23 ( const AliVParticle* track ) const;

  /// Apply also sharp pt cut when matching with trigger
  void ApplySharpPtCutInMatching ( Bool_t sharpPtCut = kTRUE ) { fSharpPtCut = sharpPtCut; }
  /// Get flag to apply the sharp pt cut when matching with trigger
  Bool_t IsApplySharpPtCutInMatching () const { return fSharpPtCut; }
  
  /// Get track cuts param (you're not supposed to modify its content
  const AliOADBMuonTrackCutsParam GetMuonTrackCutsParam () const { return fOADBParam; };
  AliOADBMuonTrackCutsParam* CustomParam ();

 private:
  
  Bool_t ReadParamFromOADB ( Int_t runNumber, Int_t passNumber );

  Bool_t fIsMC;               ///< Monte Carlo analysis
  Bool_t fUseCustomParam;     ///< Use custom parameters (do not search in OADB)
  Bool_t fSharpPtCut;         ///< Flag to apply sharp pt cut in track-trigger matching
  Bool_t fAllowDefaultParams; ///< Flag to allow default parameters from OADB
  Int_t fPassNumber;          ///< Pass number
  AliOADBMuonTrackCutsParam fOADBParam; ///< Track param in OADB

  ClassDef(AliMuonTrackCuts, 4); // Class for muon track filters
};
 
#endif

