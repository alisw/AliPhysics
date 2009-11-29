#ifndef ALIITSQASDDCHECKER_H
#define ALIITSQASDDCHECKER_H
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
#include "AliITSCalibrationSDD.h"
class AliITSLoader ; 

class AliITSQASDDChecker: public TObject {

public:
  AliITSQASDDChecker():fSubDetOffset(0),fStepBitSDD(NULL),fLowSDDValue(NULL),fHighSDDValue(NULL),fCalibration(NULL) {;}          // ctor
  AliITSQASDDChecker& operator = (const AliITSQASDDChecker& qac) ; //operator =
  virtual ~AliITSQASDDChecker(); /*{if(fStepBitSDD) delete[] fStepBitSDD ;if(fLowSDDValue)delete[]fLowSDDValue;if(fHighSDDValue) delete[]fHighSDDValue;if(fCalibration)delete fCalibration;} */// dtor
  virtual Double_t Check(AliQAv1::ALITASK_t index, TObjArray * list);
  virtual void SetTaskOffset(Int_t taskoffset);
  virtual void SetStepBit(Double_t *steprange);
  virtual Double_t *GetStepBit(){return fStepBitSDD;};
  virtual void SetSDDLimits(Float_t *lowvalue, Float_t * highvalue);
private:
  AliITSQASDDChecker(const AliITSQASDDChecker& qac):TObject(),fSubDetOffset(qac.fSubDetOffset),fStepBitSDD(qac.fStepBitSDD),fLowSDDValue(qac.fLowSDDValue),fHighSDDValue(qac.fHighSDDValue),fCalibration(qac.fCalibration) {;} // cpy ctor   
  Int_t fSubDetOffset;            // checking operation starting point
  Double_t *fStepBitSDD;
  Float_t *fLowSDDValue;
  Float_t *fHighSDDValue;
  TObjArray *fCalibration;

  static const Int_t fgknSDDmodules = 260; // number of SDD modules
  static const Int_t fgkmodoffset = 240;   // number of SPD modules
  ClassDef(AliITSQASDDChecker,2)  // description 

};

#endif // AliITSQASDDChecker_H
