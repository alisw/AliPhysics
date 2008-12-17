#ifndef ALIITSQACHECKER_H
#define ALIITSQACHECKER_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//
//  Checks the quality assurance. 
//  By comparing with reference data
//  INFN Torino
//  W.Ferrarese  P.Cerello  Mag 2008
//


// --- ROOT system ---
class TFile ; 
class TH2F ;  

// --- AliRoot header files ---
#include "AliQA.h"
#include "AliQACheckerBase.h"

class AliITSQASPDChecker;
class AliITSQASDDChecker;
class AliITSQASSDChecker;
class AliITSLoader ; 

class AliITSQAChecker: public AliQACheckerBase {

friend class AliITSQASPDChecker;
friend class AliITSQASDDChecker;
friend class AliITSQASSDChecker;

public:
  AliITSQAChecker(Bool_t kMode = kFALSE, Short_t subDet = 0, Short_t ldc = 0) ;         // ctor
  AliITSQAChecker(const AliITSQAChecker& qac);
  AliITSQAChecker& operator=(const AliITSQAChecker& qac);  
  virtual ~AliITSQAChecker() {;} // dtor
  void SetMode(Bool_t kMode) { fkOnline = kMode; }
  void SetSubDet(Short_t subdet) { fDet = subdet; }
  void SetLDC(Short_t ldc) { fLDC = ldc; }
  Bool_t GetMode() { return fkOnline; }
  Short_t GetSubDet() { return fDet; }
  Short_t GetLDC() { return fLDC; }
  virtual void SetTaskOffset(Int_t SPDOffset, Int_t SDDOffset, Int_t SSDOffset);
  virtual void SetDetTaskOffset(Int_t subdet=0,Int_t offset=0);

protected:

  virtual void SetSPDTaskOffset(Int_t SPDOffset){fSPDOffset = SPDOffset;} ;
  virtual void SetSDDTaskOffset(Int_t SDDOffset){fSDDOffset = SDDOffset;} ;
  virtual void SetSSDTaskOffset(Int_t SSDOffset){fSSDOffset = SSDOffset;} ;

  virtual Double_t Check(AliQA::ALITASK_t /*index*/){return 0.5;}
  virtual Double_t Check(AliQA::ALITASK_t index, TObjArray * list ) ;
  Double_t Check(AliQA::ALITASK_t, TNtupleD*) {AliFatal("Not implemented\n"); return 0;}
private:

  Bool_t  fkOnline;
  Short_t fDet;  
  Short_t fLDC;
  Int_t fSPDOffset; //starting point for the QACheck list
  Int_t fSDDOffset;
  Int_t fSSDOffset;

  AliITSQASPDChecker *fSPDChecker;  // SPD Checker
  AliITSQASDDChecker *fSDDChecker;  // SDD Checker
  AliITSQASSDChecker *fSSDChecker;  // SSD Checker

  ClassDef(AliITSQAChecker,3)  // description 

};

#endif // AliITSQAChecker_H
