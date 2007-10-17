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
//-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include<TObject.h>
#include <TArrayD.h>

class TLorentzVector ;
class TParticle ;
class TList ;
class TH2F ;

class AliNeutralMesonSelection : public TObject {

 public: 
  
  AliNeutralMesonSelection() ; // default ctor
  AliNeutralMesonSelection(const AliNeutralMesonSelection & g) ; // cpy ctor
  AliNeutralMesonSelection & operator = (const AliNeutralMesonSelection & g) ;//cpy assignment
  virtual ~AliNeutralMesonSelection() ; //virtual dtor

  enum type_t {kSelectPhiMinPt, kSelectPhiPtRatio, kNoSelectPhiPt};

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
  
  Int_t GetPhiPtSelection(){  return fSelect ; }
  void SetPhiPtSelection(Int_t ana ){  fSelect = ana ; }

  Bool_t AreNeutralMesonSelectionHistosKept() { return fKeepNeutralMesonHistos ; }
  void KeepNeutralMesonSelectionHistos(Bool_t keep) { fKeepNeutralMesonHistos = keep ; }

  void InitParameters();	
  Bool_t IsAngleInWindow(const Float_t angle, const Float_t e);
  void Print(const Option_t * opt) const;

  Bool_t  CutPtPhi(Double_t ptg, Double_t phig, Double_t pt, Double_t phi) ;
  Bool_t  SelectPair(TParticle * photon, TLorentzVector particlei,  TLorentzVector particlej)  ;
  
  private:
  Int_t fSelect; //Pair selection depends on analysis
  Double_t fM ; //mass of the neutral meson
  Double_t   fInvMassMaxCut ;  // Invariant Mass cut maximum
  Double_t   fInvMassMinCut ;  // Invariant Masscut minimun
  TArrayD    fAngleMaxParam ; //Max opening angle selection parameters
  Double_t   fMinPt;       // Minimum pt 
  Double_t   fDeltaPhiMaxCut ;      // 
  Double_t   fDeltaPhiMinCut ;      // 
  Double_t   fRatioMaxCut ;    // Leading particle/gamma Ratio cut maximum
  Double_t   fRatioMinCut ;    // Leading particle/gamma Ratio cut minimum
  Bool_t  fKeepNeutralMesonHistos ; // Keep neutral meson selection histograms

  //Histograms
  TH2F * fhAnglePairNoCut  ; 
  TH2F * fhAnglePairCorrelationCut  ; 
  TH2F * fhAnglePairOpeningAngleCut   ; 
  TH2F * fhAnglePairAllCut   ; 
  TH2F * fhInvMassPairNoCut    ; 
  TH2F * fhInvMassPairCorrelationCut    ;
  TH2F * fhInvMassPairOpeningAngleCut  ; 
  TH2F * fhInvMassPairAllCut   ; 
  
  ClassDef(AliNeutralMesonSelection,1)
    
    } ;


#endif //ALINEUTRALMESONSELECTION_H



