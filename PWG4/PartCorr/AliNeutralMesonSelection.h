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
class Riostream ;

//--- ANALYSIS system ---
class AliLog ;

class AliNeutralMesonSelection : public TObject {

 public: 
  
  AliNeutralMesonSelection() ; // default ctor
  AliNeutralMesonSelection(const AliNeutralMesonSelection & g) ; // cpy ctor
  AliNeutralMesonSelection & operator = (const AliNeutralMesonSelection & g) ;//cpy assignment
  virtual ~AliNeutralMesonSelection() ; //virtual dtor

  TList * GetCreateOutputObjects();
  
  Double_t GetAngleMaxParam(Int_t i) const {return fAngleMaxParam.At(i) ; }
  void SetAngleMaxParam(Int_t i, Double_t par){fAngleMaxParam.AddAt(par,i) ; }
  
  Double_t GetInvMassMaxCut() const {return fInvMassMaxCut ; }
  Double_t GetInvMassMinCut() const {return fInvMassMinCut ; }
  void SetInvMassCutRange(Double_t invmassmin, Double_t invmassmax)
  {fInvMassMaxCut =invmassmax;  fInvMassMinCut =invmassmin;}	

  Double_t GetMass() const {return fM ; }
  void SetMass(Double_t m) { fM =m ; }

  Bool_t AreNeutralMesonSelectionHistosKept() const { return fKeepNeutralMesonHistos ; }
  void KeepNeutralMesonSelectionHistos(Bool_t keep) { fKeepNeutralMesonHistos = keep ; }

  void InitParameters();	
  Bool_t IsAngleInWindow(const Float_t angle, const Float_t e) const ;
  void Print(const Option_t * opt) const;

  Bool_t  SelectPair(TLorentzVector particlei,  TLorentzVector particlej)  ;
  
  private:
  Double_t fM ; //mass of the neutral meson
  Double_t   fInvMassMaxCut ;  // Invariant Mass cut maximum
  Double_t   fInvMassMinCut ;  // Invariant Masscut minimun
  TArrayD    fAngleMaxParam ; //Max opening angle selection parameters

  Bool_t  fKeepNeutralMesonHistos ; // Keep neutral meson selection histograms

  //Histograms
  TH2F * fhAnglePairNoCut  ;  //Aperture angle of decay photons, no cuts
  TH2F * fhAnglePairOpeningAngleCut   ;  //Aperture angle of decay photons, cut on opening angle
  TH2F * fhAnglePairAllCut   ;  //Aperture angle of decay photons, all cuts
  TH2F * fhInvMassPairNoCut    ;  //Invariant mass of decay photons, no cuts
  TH2F * fhInvMassPairOpeningAngleCut  ;  //Invariant mass of decay photons, cut on opening angle
  TH2F * fhInvMassPairAllCut   ;  //Invariant mass of decay photons, all cuts
  
  ClassDef(AliNeutralMesonSelection,1)
    
    } ;


#endif //ALINEUTRALMESONSELECTION_H



