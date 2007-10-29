#ifndef ALIANAGAMMACORRELATION_H
#define ALIANAGAMMACORRELATION_H
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
// Base class for  gamma correlations 
//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include <TParticle.h> 
#include <TClonesArray.h> 
#include <TH2F.h>
#include<TObject.h>
#include<TList.h>

class AliNeutralMesonSelection ;

class AliAnaGammaCorrelation : public TObject {

public: 
  
  AliAnaGammaCorrelation() ; // default ctor
  AliAnaGammaCorrelation(const AliAnaGammaCorrelation & g) ; // cpy ctor
  AliAnaGammaCorrelation & operator = (const AliAnaGammaCorrelation & g) ;//cpy assignment
  virtual ~AliAnaGammaCorrelation() ; //virtual dtor

  enum Corrtype {kParton, kHadron, kJetLeadCone, kJetFinder};

  //General methods

  AliNeutralMesonSelection * GetNeutralMesonSelection() 
  { return fNeutralMesonSelection ; }
  void SetNeutralMesonSelection(AliNeutralMesonSelection * nms) 
  { fNeutralMesonSelection = nms ; } 

  TList * GetOutputContainer() const {return fOutputContainer ;} 
  void SetOutputContainer(TList * oc) {fOutputContainer = oc ;}  

  void InitParameters();

  Int_t GetCorrelationType() const {  return fCorrelationType ; }
  void SetCorrelationType(Int_t ana ){  fCorrelationType = ana ; }

  void Print(const Option_t * opt) const;
 
  Bool_t     AreJetsOnlyInCTS() const {return fJetsOnlyInCTS ; } 
  void SetJetsOnlyInCTS(Bool_t opt){fJetsOnlyInCTS = opt; }

  virtual TList * GetCreateOutputObjects() {return fOutputContainer ;}
  virtual void MakeGammaCorrelation(TParticle * ,  TClonesArray *, TClonesArray *)  {;}

  //Gamma hadron correlations methods: kHadron
  Float_t    GetMinPtHadron() const {return fMinPtHadron ; }
  void SetMinPtHadron(Float_t pt){fMinPtHadron = pt; };
  
  Double_t GetDeltaPhiMaxCut() const {return fDeltaPhiMaxCut ; }
  Double_t GetDeltaPhiMinCut() const {return fDeltaPhiMinCut ; }
  void SetDeltaPhiCutRange(Double_t phimin, Double_t phimax)
  {fDeltaPhiMaxCut =phimax;  fDeltaPhiMinCut =phimin;}

  Double_t GetRatioMaxCut() const {return fRatioMaxCut ; }
  Double_t GetRatioMinCut() const {return fRatioMinCut ; }
  void SetRatioCutRange(Double_t ratiomin, Double_t ratiomax)
  {fRatioMaxCut = ratiomax;  fRatioMinCut = ratiomin;}
  
  private:
  
  TList * fOutputContainer; //Histograms container
  AliNeutralMesonSelection *  fNeutralMesonSelection ; //! Pointer to pair selection for pi0 identification.

  Int_t  fCorrelationType; //Type of correlation analysis
  Bool_t   fJetsOnlyInCTS ;    // Jets measured only in TPC+ITS.
  
  private:
  //Gamma hadron correlations data members kGammaHadron
  Double_t   fMinPtHadron;       // Minimum pt of hadron (kHadron)
  Double_t   fDeltaPhiMaxCut ;      // Minimum Delta Phi Gamma-Hadron/jet in leading cone
  Double_t   fDeltaPhiMinCut ;      //  Maximum Delta Phi Gamma-Hadron/ jet in leading cone
  Double_t   fRatioMaxCut ;    // Leading particle/gamma Ratio cut maximum (kLeadJetCone)
  Double_t   fRatioMinCut ;    // Leading particle/gamma Ratio cut minimum (kLeadJetCone)
  
  ClassDef(AliAnaGammaCorrelation,1)
} ;
 

#endif //ALIANAGAMMACORRELATION_H



