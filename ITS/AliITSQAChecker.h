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
#include "AliQAv1.h"
#include "AliQACheckerBase.h"

class AliITSQASPDChecker;
class AliITSQASDDChecker;
class AliITSQASSDChecker;
class AliITSLoader ; 

class AliITSQAChecker: public AliQACheckerBase {

public:
  AliITSQAChecker(Bool_t kMode = kFALSE, Short_t subDet = 0, Short_t ldc = 0) ;         // ctor
  AliITSQAChecker(const AliITSQAChecker& qac);
  AliITSQAChecker& operator=(const AliITSQAChecker& qac);  
  virtual ~AliITSQAChecker();// dtor
  void SetMode(Bool_t kMode) { fkOnline = kMode; }
  void SetSubDet(Short_t subdet) { fDet = subdet; }
  void SetLDC(Short_t ldc) { fLDC = ldc; }
  Bool_t GetMode() { return fkOnline; }
  Short_t GetSubDet() { return fDet; }
  Short_t GetLDC() { return fLDC; }
  virtual void SetTaskOffset(Int_t SPDOffset, Int_t SDDOffset, Int_t SSDOffset);
  virtual void SetHisto(Int_t SPDhisto, Int_t SDDhisto, Int_t SSDhisto);
  virtual void SetDetTaskOffset(Int_t subdet=0,Int_t offset=0);
  virtual void InitQACheckerLimits();
  virtual void CreateStepForBit(Double_t histonumb,Double_t *steprange);
  virtual void SetQA(AliQAv1::ALITASK_t index, Double_t * value) const;
  virtual void SetDetHisto(Int_t subdet=0,Int_t histo=0);

  virtual Int_t GetSPDHisto(){return fSPDHisto;} ;
  virtual Int_t GetSDDHisto(){return fSDDHisto;} ;
  virtual Int_t GetSSDHisto(){return fSSDHisto;} ;


protected:
  virtual void Check(Double_t * test, AliQAv1::ALITASK_t index, TObjArray ** list, const AliDetectorRecoParam * recoParam) ;
  virtual void SetSPDTaskOffset(Int_t SPDOffset){fSPDOffset = SPDOffset;} ;
  virtual void SetSDDTaskOffset(Int_t SDDOffset){fSDDOffset = SDDOffset;} ;
  virtual void SetSSDTaskOffset(Int_t SSDOffset){fSSDOffset = SSDOffset;} ;

  virtual void SetSPDHisto(Int_t SPDhisto){fSPDHisto = SPDhisto;} ;
  virtual void SetSDDHisto(Int_t SDDhisto){fSDDHisto = SDDhisto;} ;
  virtual void SetSSDHisto(Int_t SSDhisto){fSSDHisto = SSDhisto;} ;

private:

  Bool_t  fkOnline;
  Short_t fDet;  
  Short_t fLDC;

  Int_t fSPDOffset; //starting point for the QACheck list
  Int_t fSDDOffset;
  Int_t fSSDOffset;

  Int_t fSPDHisto;
  Int_t fSDDHisto;
  Int_t fSSDHisto;

  AliITSQASPDChecker *fSPDChecker;  // SPD Checker
  AliITSQASDDChecker *fSDDChecker;  // SDD Checker
  AliITSQASSDChecker *fSSDChecker;  // SSD Checker

  ClassDef(AliITSQAChecker,3)  // description 

};

#endif // AliITSQAChecker_H
