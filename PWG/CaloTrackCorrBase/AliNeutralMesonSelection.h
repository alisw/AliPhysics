#ifndef ALINEUTRALMESONSELECTION_H
#define ALINEUTRALMESONSELECTION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliNeutralMesonSelection
/// \ingroup CaloTrackCorrelationsBase
/// \brief Class that contains methods to select candidate cluster pairs to neutral meson.
///
/// Class that contains methods to select candidate pairs to neutral meson. 
/// 2 main selections, invariant mass around pi0 (also any other mass),
/// apperture angle to distinguish from combinatorial.
/// There is a 3rd cut based on the gamma correlation on phi or pt.
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//_________________________________________________________________________

// --- ROOT system ---
#include<TObject.h>
#include<TArrayD.h>
#include<TString.h>

class TLorentzVector ;
class TList   ;
class TH2F    ;

class AliNeutralMesonSelection : public TObject {
  
 public:
    
  AliNeutralMesonSelection() ; // default ctor
  
  /// Virtual destructor.
  virtual ~AliNeutralMesonSelection() { ; }   

  // General
    
  TList *  GetCreateOutputObjects();
  
  void     InitParameters();	
  
  void     Print(const Option_t * opt) const;

  Bool_t   AreNeutralMesonSelectionHistosKept()   const { return fKeepNeutralMesonHistos ; }
  void     KeepNeutralMesonSelectionHistos(Bool_t keep) { fKeepNeutralMesonHistos = keep ; }
  
  Bool_t   SelectPair(TLorentzVector particlei,  TLorentzVector particlej, Int_t calo)  ;
  
  void     SetParticle(TString particleName) ;  // Do some default settings for "Pi0" or "Eta"
  TString  GetParticle()                          const { return fParticle               ; }
  
  Int_t  GetDebug()                               const { return fDebug                  ; }
  void   SetDebug(Int_t d)                              { fDebug = d                     ; }
  
  // Asymmetry selection
    
  Float_t  GetAsymmetryCut()                      const { return fAsymmetryCut           ; }
  void     SetAsymmetryCut(Float_t asy)                 { fAsymmetryCut = asy            ; }
  
  void     SwitchOnAsymmetryCut()                       { fUseAsymmetryCut = kTRUE       ; }
  void     SwitchOffAsymmetryCut()                      { fUseAsymmetryCut = kFALSE      ; }
  
  //Opening angle selection 
  
  Double_t GetAngleMaxParam(Int_t i)              const { return fAngleMaxParam.At(i)    ; }
  void     SetAngleMaxParam(Int_t i, Double_t par)      { fAngleMaxParam.AddAt(par,i)    ; }
  
  void     SetShiftMinAngleCut(Float_t a, Float_t b)    { fShiftMinAngle[0] = a          ;
                                                          fShiftMinAngle[1] = b          ; }
  
  void     SwitchOnAngleSelection()                     { fUseAngleCut = kTRUE           ; }
  void     SwitchOffAngleSelection()                    { fUseAngleCut = kFALSE          ; }

  Bool_t   IsAngleInWindow(Float_t angle, Float_t e) const ;
  
  //Invariant mass selection
  
  Double_t GetInvMassMaxCut()                     const { return fInvMassMaxCut          ; }
  Double_t GetInvMassMinCut()                     const { return fInvMassMinCut          ; }
  
  void     SetInvMassCutRange(Double_t invmassmin, Double_t invmassmax)
            { fInvMassMaxCut =invmassmax;  fInvMassMinCut =invmassmin                    ; }	
  
  void     SetSideBandCutRanges( Double_t lmin, Double_t lmax, 
                                 Double_t rmin, Double_t rmax )
            { fLeftBandMinCut  = lmin ; fLeftBandMaxCut  = lmax ; 
              fRightBandMinCut = rmin ; fRightBandMaxCut = rmax ; }	
  
  void     SetInvMassCutMaxParameters(Float_t a, Float_t b, Float_t c)
            { fInvMassMaxCutParam[0] = a ; 
              fInvMassMaxCutParam[1] = b ; 
              fInvMassMaxCutParam[2] = c ; }	

  Double_t GetMass()                              const { return fM                      ; }
  void     SetMass(Double_t m)                          { fM = m                         ; }
  
  
  // Decay photon bit setting
  static const Int_t fgkMaxNDecayBits = 8;
  
  enum decayTypes { kPi0          = 0, kEta          = 1,
                    kPi0RightSide = 2, kEtaRightSide = 3, kEtaLeftSide  = 4,
                    kPi0LeftSide  = 5, kEtaBothSides = 6, kPi0BothSides = 7 } ;

  UInt_t    GetDecayBit()                         const { return fDecayBit               ; }
  
  void    SetDecayBit(Int_t &tag, UInt_t set) const {
    // Set bit of type set (decayTypes) in tag
    tag |= (1<<set) ;
  }
  
  void    SetDecayBit(Int_t &tag) const {
    // Set bit of type set (decayTypes) in tag
    tag |= (1<<fDecayBit) ;
  }
  
  Bool_t  CheckDecayBit(Int_t tag, UInt_t test) const {
    // Check if in tag the bit test (decayTypes) is set.
    if (tag & (1<<test) ) return  kTRUE ;
    else return kFALSE ;
  }

  Bool_t  CheckDecayBit(Int_t tag) const {
    // Check if in tag the bit test (decayTypes) is set.
    if (tag & (1<<fDecayBit) ) return  kTRUE ;
    else return kFALSE ;
  }
    
  // Histograms setters and getters
  
  virtual void SetHistoERangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoNEBins = n ;
    fHistoEMax = max ;
    fHistoEMin = min ;
  }
  
  Int_t   GetHistoNEBins()     const { return fHistoNEBins    ; }
  Float_t GetHistoEMin()       const { return fHistoEMin      ; }
  Float_t GetHistoEMax()       const { return fHistoEMax      ; }
    
  virtual void SetHistoAngleRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoNAngleBins = n ;
    fHistoAngleMax = max ;
    fHistoAngleMin = min ;
  }
  
  Int_t   GetHistoNAngleBins() const { return fHistoNAngleBins ; }
  Float_t GetHistoAngleMin()   const { return fHistoAngleMin   ; }
  Float_t GetHistoAngleMax()   const { return fHistoAngleMax   ; }
  
  virtual void SetHistoIMRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoNIMBins = n ;
    fHistoIMMax = max ;
    fHistoIMMin = min ;
  }
  
  Int_t   GetHistoNIMBins()    const { return fHistoNIMBins    ; }
  Float_t GetHistoIMMin()      const { return fHistoIMMin      ; }
  Float_t GetHistoIMMax()      const { return fHistoIMMax      ; }
  
 private:
  
  Float_t  fAsymmetryCut  ;               ///<  Asymmetry cut.
  
  Bool_t   fUseAsymmetryCut;              ///<  Use the asymmetry cut.

  Double_t fM ;                           ///<  Mass of the neutral meson.
  
  Double_t fInvMassMaxCut ;               ///<  Invariant Mass cut maximum.
  
  Double_t fInvMassMinCut ;               ///<  Invariant Masscut minimun.
  
  Double_t fInvMassMaxCutParam[3];        ///<  Variable invariant mass max cut, for pi0 in EMCAL.
  
  Double_t fLeftBandMinCut  ;             ///<  Side Band selection, min left  band cut.
  
  Double_t fLeftBandMaxCut  ;             ///<  Side Band selection, max left  band cut.
  
  Double_t fRightBandMinCut ;             ///<  Side Band selection, min right band cut.
  
  Double_t fRightBandMaxCut ;             ///<  Side Band selection, max right band cut.
  
  TArrayD  fAngleMaxParam ;               ///<  Maximum opening angle selection parameters.
  
  Bool_t   fUseAngleCut   ;               ///<  Select pairs depending on their opening angle.
  
  Float_t  fShiftMinAngle[2] ;            ///<  Correction shift for min angle from true kinematic limit, resolution effects.
  
  Bool_t   fKeepNeutralMesonHistos ;      ///<  Keep neutral meson selection histograms.
  
  TString  fParticle ;                    ///<  Meutral meson name (Pi0, Eta, +SideBand).
  
  UInt_t   fDecayBit;                     ///<  Decay type flag, set while selecting, depending on fParticle and side range. See enum decayTypes for possible bits.

  Int_t    fDebug ;                       ///< Debug level.
    
  // Histograms
  TH2F *   fhAnglePairNoCut ;             //!<! Aperture angle of decay photons, no cuts.
  
  TH2F *   fhAnglePairOpeningAngleCut ;   //!<! Aperture angle of decay photons, cut on opening angle.
  
  TH2F *   fhAnglePairAsymmetryCut ;      //!<! Aperture angle of decay photons, asymmetry cut.
  
  TH2F *   fhAnglePairAllCut ;            //!<! Aperture angle of decay photons, all cuts.
  
  
  TH2F *   fhInvMassPairNoCut ;           //!<! Invariant mass of decay photons, no cuts.
  
  TH2F *   fhInvMassPairOpeningAngleCut ; //!<! Invariant mass of decay photons, cut on opening angle.
  
  TH2F *   fhInvMassPairAsymmetryCut ;    //!<! Invariant mass of decay photons, asymmetry cut.  
  
  TH2F *   fhInvMassPairAllCut ;          //!<! Invariant mass of decay photons, all cuts.

  TH2F *   fhAsymmetryNoCut ;             //!<! Asymmetry of decay photons, no cuts.
  
  TH2F *   fhAsymmetryOpeningAngleCut ;   //!<! Asymmetry of decay photons, cut on opening angle.
  
  TH2F *   fhAsymmetryAllCut ;            //!<! Asymmetry of decay photons, all cuts.
  
  
  // Histograms binning and range
  Int_t    fHistoNEBins ;                 ///< Number of bins in pi0 E axis.
  
  Float_t  fHistoEMax ;                   ///< Maximum value of pi0 E histogram range.
  
  Float_t  fHistoEMin ;                   ///< Minimum value of pi0 E histogram range.
  
  Int_t    fHistoNAngleBins ;             ///< Number of bins in angle axis.
  
  Float_t  fHistoAngleMax ;               ///< Maximum value of angle histogram range.
  
  Float_t  fHistoAngleMin ;               ///< Minimum value of angle histogram range.
  
  
  Int_t    fHistoNIMBins ;                ///< Number of bins in Invariant Mass axis.
  
  Float_t  fHistoIMMax ;                  ///< Maximum value of Invariant Mass histogram range.
  
  Float_t  fHistoIMMin ;                  ///< Minimum value of Invariant Mass histogram range. 
  
  /// Copy constructor not implemented.
  AliNeutralMesonSelection(              const AliNeutralMesonSelection & nm) ;
  
  /// Assignment operator not implemented.
  AliNeutralMesonSelection & operator = (const AliNeutralMesonSelection & nm) ; 
  
  /// \cond CLASSIMP
  ClassDef(AliNeutralMesonSelection,8) ;
  /// \endcond

} ;

#endif //ALINEUTRALMESONSELECTION_H



