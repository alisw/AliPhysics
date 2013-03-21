#ifndef ALIANARANDOMTRIGGER_H
#define ALIANARANDOMTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
// Gerenate a random trigger, input for other analysis
// Set flat energy distribution over acceptance of EMCAL, PHOS or CTS
// Be careful, correlate only with Min Bias events this random trigger particle
//
//
//-- Author: Gustavo Conesa (LPSC-Grenoble)

// Root system
class TH2F; 
#include <TRandom3.h>

// Analysis system
#include "AliAnaCaloTrackCorrBaseClass.h"
 
class AliAnaRandomTrigger : public AliAnaCaloTrackCorrBaseClass {
  
 public: 
  AliAnaRandomTrigger() ; // default ctor
  virtual ~AliAnaRandomTrigger() { ; } //virtual dtor

  Bool_t       ExcludeDeadBadRegions(const Float_t eta, const Float_t phi);
  
  TObjString * GetAnalysisCuts();
  
  TList      * GetCreateOutputObjects();
    
  void         InitParameters();
    
  void         MakeAnalysisFillAOD()  ;
  
  void         MakeAnalysisFillHistograms() ; 
  
  void         Print(const Option_t * opt) const;
  
  void         SetDetector(TString detector)       { fDetector  = detector ; }  
  
  void         SetEtaCut(Float_t min, Float_t max) { fEtaCut[0] = min ; fEtaCut[1] = max;}
  
  void         SetPhiCut(Float_t min, Float_t max) { fPhiCut[0] = min ; fPhiCut[1] = max;} // radians
  
  void         SetNumberOfRandomParticles(Int_t n) { fNRandom   = n   ; }
  
 private:
  
  TString    fDetector ; // Detector : EMCAL, PHOS, CTS
  Float_t    fEtaCut[2]; // Eta acceptance
  Float_t    fPhiCut[2]; // Phi acceptance, radians
  TRandom3   fRandom   ; // Random generator
  Int_t      fNRandom  ; // Number of random particles per event
  
  //Constrol histograms 
  TH1F     * fhE;        //! E distribution
  TH1F     * fhPt;       //! pT distribution
  TH2F     * fhPhi;      //! phi distribution vs pT, negative
  TH2F     * fhEta;      //! eta distribution vs pT, negative
  TH2F     * fhEtaPhi;   //! eta vs phi distribution of positive charge  


  AliAnaRandomTrigger(              const AliAnaRandomTrigger & r) ; // cpy ctor
  AliAnaRandomTrigger & operator = (const AliAnaRandomTrigger & r) ; //cpy assignment
  
  ClassDef(AliAnaRandomTrigger,2)

} ;


#endif //ALIANARANDOMTRIGGER_H



