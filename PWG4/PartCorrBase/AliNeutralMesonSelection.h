#ifndef ALINEUTRALMESONSELECTION_H
#define ALINEUTRALMESONSELECTION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: AliNeutralMesonSelection.h 27413 2008-07-18 13:28:12Z gconesab $ */

//_________________________________________________________________________
// Class that contains methods to select candidate pairs to neutral meson 
// 2 main selections, invariant mass around pi0 (also any other mass),
// apperture angle to distinguish from combinatorial.
// There is a 3rd cut based on the gamma correlation on phi or pt.
//-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include<TObject.h>
#include<TArrayD.h>

class TLorentzVector ;
class TList ;
class TH2F ;

//--- ANALYSIS system ---

class AliNeutralMesonSelection : public TObject {
  
 public: 
  AliNeutralMesonSelection() ; // default ctor
  virtual ~AliNeutralMesonSelection() { ; } //virtual dtor  

  // General

  TList *  GetCreateOutputObjects();
  
  void     InitParameters();	
  
  void     Print(const Option_t * opt) const;

  Bool_t   AreNeutralMesonSelectionHistosKept()   const { return fKeepNeutralMesonHistos ; }
  void     KeepNeutralMesonSelectionHistos(Bool_t keep) { fKeepNeutralMesonHistos = keep ; }
  
  Bool_t   SelectPair(TLorentzVector particlei,  TLorentzVector particlej, TString calo)  ;
  
  void     SetParticle(TString particleName) ;  // Do some default settings for "Pi0" or "Eta"
  
  // Asymmetry selection
    
  Float_t  GetAsymmetryCut()                     const { return fAsymmetryCut            ; }
  void     SetAsymmetryCut(Float_t asy)                { fAsymmetryCut = asy             ; }	
  
  void     SwitchOnAsymmetryCut()                      { fUseAsymmetryCut = kTRUE        ; }
  void     SwitchOffAsymmetryCut()                     { fUseAsymmetryCut = kFALSE       ; }
  
  //Opening angle selection 
  
  Double_t GetAngleMaxParam(Int_t i)              const { return fAngleMaxParam.At(i)    ; }
  void     SetAngleMaxParam(Int_t i, Double_t par)      { fAngleMaxParam.AddAt(par,i)    ; }
  
  void     SetShiftMinAngleCut(Float_t a, Float_t b)    { fShiftMinAngle[0] = a          ; 
                                                          fShiftMinAngle[1] = b          ; }
  
  void     SwitchOnAngleSelection()                     { fUseAngleCut = kTRUE           ; }
  void     SwitchOffAngleSelection()                    { fUseAngleCut = kFALSE          ; }

  Bool_t   IsAngleInWindow(const Float_t angle, const Float_t e) const ;
  
  //Invariant mass selection
  
  Double_t GetInvMassMaxCut()                     const { return fInvMassMaxCut          ; }
  Double_t GetInvMassMinCut()                     const { return fInvMassMinCut          ; }
  
  void     SetInvMassCutRange(Double_t invmassmin, Double_t invmassmax)
            { fInvMassMaxCut =invmassmax;  fInvMassMinCut =invmassmin                    ; }	
  
  void     SetInvMassCutMaxParameters(Float_t a, Float_t b, Float_t c)
            { fInvMassMaxCutParam[0] = a ; 
              fInvMassMaxCutParam[1] = b ; 
              fInvMassMaxCutParam[2] = c ; }	

  Double_t GetMass()                              const { return fM                      ; }
  void     SetMass(Double_t m)                          { fM = m                         ; }
      
  //Histogrammes setters and getters
  
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
  
  Float_t  fAsymmetryCut  ;               // Asymmetry cut
  Bool_t   fUseAsymmetryCut;              // Use the asymmetry cut
  
  Double_t fM ;                           // Mass of the neutral meson
  Double_t fInvMassMaxCut ;               // Invariant Mass cut maximum
  Double_t fInvMassMinCut ;               // Invariant Masscut minimun
  Double_t fInvMassMaxCutParam[3];        // Variable invariant mass max cut, for pi0 in EMCAL
  
  TArrayD  fAngleMaxParam ;               // Max opening angle selection parameters
  Bool_t   fUseAngleCut   ;               // Select pairs depending on their opening angle
  Float_t  fShiftMinAngle[2] ;            // Correction shift for min angle from true kinematic limit, resolution effects
  
  Bool_t   fKeepNeutralMesonHistos ;      // Keep neutral meson selection histograms

  //Histograms
  TH2F *   fhAnglePairNoCut ;             //! Aperture angle of decay photons, no cuts
  TH2F *   fhAnglePairOpeningAngleCut ;   //! Aperture angle of decay photons, cut on opening angle
  TH2F *   fhAnglePairAsymmetryCut ;      //! Aperture angle of decay photons, asymmetry cut
  TH2F *   fhAnglePairAllCut ;            //! Aperture angle of decay photons, all cuts
  
  TH2F *   fhInvMassPairNoCut ;           //! Invariant mass of decay photons, no cuts
  TH2F *   fhInvMassPairOpeningAngleCut ; //! Invariant mass of decay photons, cut on opening angle
  TH2F *   fhInvMassPairAsymmetryCut ;    //! Invariant mass of decay photons, asymmetry cut  
  TH2F *   fhInvMassPairAllCut ;          //! Invariant mass of decay photons, all cuts

  TH2F *   fhAsymmetryNoCut ;             //! Asymmetry of decay photons, no cuts
  TH2F *   fhAsymmetryOpeningAngleCut ;   //! Asymmetry of decay photons, cut on opening angle
  TH2F *   fhAsymmetryAllCut ;            //! Asymmetry of decay photons, all cuts
  
  //Histograms binning and range    
  Int_t    fHistoNEBins ;                 // Number of bins in pi0 E axis
  Float_t  fHistoEMax ;                   // Maximum value of pi0 E histogram range
  Float_t  fHistoEMin ;                   // Minimum value of pi0 E histogram range
  
  Int_t    fHistoNAngleBins ;             // Number of bins in angle axis
  Float_t  fHistoAngleMax ;               // Maximum value of angle histogram range
  Float_t  fHistoAngleMin ;               // Minimum value of angle histogram range
  
  Int_t    fHistoNIMBins ;                // Number of bins in Invariant Mass axis
  Float_t  fHistoIMMax ;                  // Maximum value of Invariant Mass histogram range
  Float_t  fHistoIMMin ;                  // Minimum value of Invariant Mass histogram range  
  
  AliNeutralMesonSelection(const AliNeutralMesonSelection & g) ;               // cpy ctor
  AliNeutralMesonSelection & operator = (const AliNeutralMesonSelection & g) ; // cpy assignment
  
  ClassDef(AliNeutralMesonSelection,6)
    
} ;

#endif //ALINEUTRALMESONSELECTION_H



