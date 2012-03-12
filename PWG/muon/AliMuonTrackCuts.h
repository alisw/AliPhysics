#ifndef ALIMUONTRACKCUTS_H
#define ALIMUONTRACKCUTS_H

#include "AliAnalysisCuts.h"
#include "TArrayD.h"

class AliVParticle;
class TList;
class TVector3;

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
  AliMuonTrackCuts(const char* name, const char* title, Bool_t isESD);
  AliMuonTrackCuts(const AliMuonTrackCuts& obj);
  AliMuonTrackCuts& operator=(const AliMuonTrackCuts& obj);

  virtual ~AliMuonTrackCuts();

  virtual UInt_t GetSelectionMask ( const TObject* obj );
  virtual Bool_t IsSelected ( TObject* obj );
  virtual Bool_t IsSelected ( TList* /*list */ );
  
  void SetDefaultFilterMask();

  Bool_t SetRun(Int_t runNumber);
  void SetUseCustomParam( Bool_t useCustomParam = kTRUE, Int_t runNumber = -1 );
  void SetIsMC(Bool_t isMC = kTRUE) { fIsMC = isMC; }

  void Print ( Option_t* option = "" ) const;

  TVector3 GetCorrectedDCA ( const AliVParticle* track ) const;
  Double_t GetAverageMomentum ( const AliVParticle* track ) const;

  enum {
    kThetaAbs23,   ///< Theta_abs between 2 and 3 degrees
    kThetaAbs310,  ///< Theta_abs between 3 and 10 degrees
    kNthetaAbs     ///< Number of theta abs bins
  };

  // Parameters
  enum {
    kMeanDcaX,        ///< Average track DCA_x
    kMeanDcaY,        ///< Average track DCA_y
    kMeanDcaZ,        ///< Average track DCA_z
    kMeanPCorr23,     ///< Average momentum correction in 2-3 deg
    kMeanPCorr310,    ///< Average momentum correction in 3-10 deg
    kSigmaPdca23,     ///< Sigma_PxDCA in 2-3 deg
    kSigmaPdca310,    ///< Sigma_PxDCA in 3-10 deg
    kNSigmaPdcaCut,   ///< Cut value in units of sigma_PxDCA
    kChi2NormCut,     ///< Cut on the normalized chi2 of track
    kRelPResolution,  ///< Relative momentum resolution
    kSlopeResolution, ///< Slope resolution
    kSharpPtApt,      ///< Sharp tracker pt cut for Apt
    kSharpPtLpt,      ///< Sharp tracker pt cut for Lpt
    kSharpPtHpt,      ///< Sharp tracker pt cut for Hpt
    kNParameters      ///< Total number of parameters
  };

  void SetMeanDCA ( Double_t xAtDca, Double_t yAtDca, Double_t zAtDca = 0.);
  TVector3 GetMeanDCA () const;

  void SetMeanPCorr ( Double_t pCorrThetaAbs23, Double_t pCorrThetaAbs310 ); 
  Double_t GetMeanPCorr ( Double_t rAtAbsEnd ) const;

  void SetSigmaPdca ( Double_t sigmaThetaAbs23, Double_t sigmaThetaAbs310 );
  Double_t GetSigmaPdca ( Double_t rAtAbsEnd ) const;

  void SetNSigmaPdca ( Double_t nSigmas );
  Double_t GetNSigmaPdca () const;

  void SetChi2NormCut ( Double_t chi2normCut );
  Double_t GetChi2NormCut () const;
  
  void SetRelPResolution ( Double_t relPResolution );
  Double_t GetRelPResolution () const;
  
  void SetSlopeResolution ( Double_t slopeResolution );
  Double_t GetSlopeResolution () const;
  
  void SetSharpPtCut ( Int_t trigPtCut, Double_t ptCutValue );
  Double_t GetSharpPtCut ( Int_t trigPtCut, Bool_t warn = kTRUE ) const;

  Bool_t StreamParameters ( Int_t runNumber, Int_t maxRun );
  
  /// Apply also sharp pt cut when matching with trigger
  void ApplySharpPtCutInMatching ( Bool_t sharpPtCut = kTRUE ) { fSharpPtCut = sharpPtCut; }

 private:
  
  Int_t GetThetaAbsBin ( Double_t rAtAbsEnd ) const;
  Bool_t SetParameter ( Int_t iparam, Float_t value );
  Bool_t RunMatchesRange ( Int_t runNumber, const Char_t* objName );

  Bool_t fIsESD;            ///< Event is ESD
  Bool_t fIsMC;             ///< Monte Carlo analysis
  Bool_t fUseCustomParam;   ///< Use custom parameters (do not search in OADB)
  Bool_t fSharpPtCut;       ///< Flag to apply sharp pt cut in track-trigger matching

  TArrayD fParameters;      ///< List of parameters

  ClassDef(AliMuonTrackCuts, 2); // Class for muon track filters
};
 
#endif

