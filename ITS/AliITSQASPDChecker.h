#ifndef ALIITSQASPDCHECKER_H
#define ALIITSQASPDCHECKER_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//
//  Checks the quality assurance. 
//  By comparing with reference data
//  INFN Torino
//  P. Cerello - apr 2008
//


// --- ROOT system ---
class TFile ; 
class TH2F ;  

// --- AliRoot header files ---
#include "AliQAv1.h"
#include "AliQACheckerBase.h"
#include "AliITSQAChecker.h"
class AliITSLoader ; 

class AliITSQASPDChecker: public TObject {

public:
  AliITSQASPDChecker():fSubDetOffset(0),fStepBitSPD(NULL),fLowSPDValue(NULL),fHighSPDValue(NULL) {;}          // ctor
  AliITSQASPDChecker& operator = (const AliITSQASPDChecker& qac) ; //operator =
  virtual ~AliITSQASPDChecker() {if(fStepBitSPD) delete[] fStepBitSPD ;if(fLowSPDValue)delete[]fLowSPDValue;if(fHighSPDValue) delete[]fHighSPDValue;} // dtor
  Double_t Check(AliQAv1::ALITASK_t index, TObjArray * list);
  void SetTaskOffset(Int_t TaskOffset);

  void SetStepBit(Double_t *steprange);
  Double_t *GetStepBit(){return fStepBitSPD;};
  void SetSPDLimits(Float_t *lowvalue, Float_t * highvalue);
private:
  
  AliITSQASPDChecker(const AliITSQASPDChecker& qac):TObject(),fSubDetOffset(qac.fSubDetOffset),fStepBitSPD(qac.fStepBitSPD),fLowSPDValue(qac.fLowSPDValue),fHighSPDValue(qac.fHighSPDValue){;}  // cpy ctor   
  Int_t fSubDetOffset;            // checking operation starting point
  Double_t *fStepBitSPD;
  Float_t *fLowSPDValue;
  Float_t *fHighSPDValue;

  ClassDef(AliITSQASPDChecker,2)  // description 

};

#endif // AliITSQASPDChecker_H
