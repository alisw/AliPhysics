#ifndef ALIPHOSGammaJet_H
#define ALIPHOSGammaJet_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

//_________________________________________________________________________
//  Class for the analysis of gamma-jet correlations     
//                  
//*-- Author: Gustavo Conesa & Yves Schutz (IFIC, CERN)

// --- ROOT system ---
#include "TTask.h"

// --- AliRoot header files ---

class AliPHOSGammaJet : public TTask {

public: 

  AliPHOSGammaJet() ; // default ctor
  AliPHOSGammaJet(const TString inputfilename) ; //ctor 
  AliPHOSGammaJet(const AliPHOSGammaJet & gj) ; // cpy ctor
  virtual ~AliPHOSGammaJet() ; // dtor
  virtual void   Exec(Option_t * = ""); 
  void GetGammaJet(TList & particleList, TLorentzVector & gamma, Int_t & id) ; 
  void GetLeadingCharge(TList & particleList, TLorentzVector & charge, Int_t & id) ;
  void GetLeadingPi0(TList & particleList, TLorentzVector & pi0) ;
//    void GetLeadingGammaPair(TList &particleList, TLorentzVector &gammapair, Int_t & id, 
//  			   Double_t & thetacut,Double_t & ratiocut1, Double_t & ratiocut2,
//  			   Double_t & invmasscut1,Double_t & invmasscut2);
  void Pi0Decay(Double_t mPi0, TLorentzVector &p0, 
		TLorentzVector &p1, TLorentzVector &p2, Double_t &angle) ; 
private: 

  ClassDef(AliPHOSGammaJet,1)
} ; 
#endif //ALIPHOSGammaJet_H
