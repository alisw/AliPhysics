#ifndef ALINEUTRALMESONSELECTION_H
#define ALINEUTRALMESONSELECTION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.2  2007/08/17 12:40:04  schutz
 * New analysis classes by Gustavo Conesa
 *
 * Revision 1.1.2.1  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */
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
class TParticle ;
class TList ;
class TH2F ;
class Riostream ;

//--- AliRoot system ---
class AliLog ;

class AliNeutralMesonSelection : public TObject {

 public: 
  
  AliNeutralMesonSelection() ; // default ctor
  AliNeutralMesonSelection(const AliNeutralMesonSelection & g) ; // cpy ctor
  AliNeutralMesonSelection & operator = (const AliNeutralMesonSelection & g) ;//cpy assignment
  virtual ~AliNeutralMesonSelection() ; //virtual dtor

  enum Type {kSelectPhiMinPt, kSelectPhiPtRatio, kNoSelectPhiPt};

  TList * GetCreateOutputObjects();
  
  Double_t GetAngleMaxParam(Int_t i) const {return fAngleMaxParam.At(i) ; }
  void SetAngleMaxParam(Int_t i, Double_t par){fAngleMaxParam.AddAt(par,i) ; }
  
  Double_t GetInvMassMaxCut() const {return fInvMassMaxCut ; }
  Double_t GetInvMassMinCut() const {return fInvMassMinCut ; }
  void SetInvMassCutRange(Double_t invmassmin, Double_t invmassmax)
  {fInvMassMaxCut =invmassmax;  fInvMassMinCut =invmassmin;}	
  
  Double_t GetDeltaPhiMaxCut() const {return fDeltaPhiMaxCut ; }
  Double_t GetDeltaPhiMinCut() const {return fDeltaPhiMinCut ; }
  void SetDeltaPhiCutRange(Double_t phimin, Double_t phimax)
  {fDeltaPhiMaxCut =phimax;  fDeltaPhiMinCut =phimin;}

  Double_t GetRatioMaxCut() const {return fRatioMaxCut ; }
  Double_t GetRatioMinCut() const {return fRatioMinCut ; }
  void SetRatioCutRange(Double_t ratiomin, Double_t ratiomax)
  {fRatioMaxCut = ratiomax;  fRatioMinCut = ratiomin;}

  Float_t    GetMinPt() const {return fMinPt ; }
  void SetMinPt(Float_t pt){fMinPt = pt; };

  Double_t GetMass() const {return fM ; }
  void SetMass(Double_t m) { fM =m ; }
  
  Int_t GetPhiPtSelection() const { return fSelect ; }
  void SetPhiPtSelection(Int_t ana ){  fSelect = ana ; }

  Bool_t AreNeutralMesonSelectionHistosKept() const { return fKeepNeutralMesonHistos ; }
  void KeepNeutralMesonSelectionHistos(Bool_t keep) { fKeepNeutralMesonHistos = keep ; }

  void InitParameters();	
  Bool_t IsAngleInWindow(const Float_t angle, const Float_t e) const ;
  void Print(const Option_t * opt) const;

  Bool_t  CutPtPhi(Double_t ptg, Double_t phig, Double_t pt, Double_t phi) const ;
  Bool_t  SelectPair(TParticle * photon, TLorentzVector particlei,  TLorentzVector particlej)  ;
  
  private:
  Int_t fSelect; //Pair selection depends on analysis
  Double_t fM ; //mass of the neutral meson
  Double_t   fInvMassMaxCut ;  // Invariant Mass cut maximum
  Double_t   fInvMassMinCut ;  // Invariant Masscut minimun
  TArrayD    fAngleMaxParam ; //Max opening angle selection parameters
  Double_t   fMinPt;       // Minimum pt 
  Double_t   fDeltaPhiMaxCut ;      // Leading particle - gamma Ratio cut maximum
  Double_t   fDeltaPhiMinCut ;      // Leading particle - gamma Ratio cut maximum
  Double_t   fRatioMaxCut ;    // Leading particle/gamma Ratio cut maximum
  Double_t   fRatioMinCut ;    // Leading particle/gamma Ratio cut minimum
  Bool_t  fKeepNeutralMesonHistos ; // Keep neutral meson selection histograms

  //Histograms
  TH2F * fhAnglePairNoCut  ;  //Aperture angle of decay photons, no cuts
  TH2F * fhAnglePairCorrelationCut  ;  //Aperture angle of decay photons, cut on phi/pT correlation with prompt gamma
  TH2F * fhAnglePairOpeningAngleCut   ;  //Aperture angle of decay photons, cut on opening angle
  TH2F * fhAnglePairAllCut   ;  //Aperture angle of decay photons, all cuts
  TH2F * fhInvMassPairNoCut    ;  //Invariant mass of decay photons, no cuts
  TH2F * fhInvMassPairCorrelationCut    ;   //Invariant mass of decay photons, cut on phi/pT correlation with prompt gamma
  TH2F * fhInvMassPairOpeningAngleCut  ;  //Invariant mass of decay photons, cut on opening angle
  TH2F * fhInvMassPairAllCut   ;  //Invariant mass of decay photons, all cuts
  
  ClassDef(AliNeutralMesonSelection,1)
    
    } ;


#endif //ALINEUTRALMESONSELECTION_H



