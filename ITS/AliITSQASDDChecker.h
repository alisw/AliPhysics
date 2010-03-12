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

#include "AliQAv1.h"
// --- ROOT system ---
class TFile;
class TH2F; 

// --- AliRoot header files ---

class AliQACheckerBase;
class AliITSQAChecker;
class AliITSCalibrationSDD;
class AliITSLoader;
class TSystem;
class AliQAManager;
class AliLog;
class TF1;
class TCanvas;

class AliITSQASDDChecker: public TObject {

public:
  AliITSQASDDChecker():fSubDetOffset(0),fStepBitSDD(NULL),fLowSDDValue(NULL),fHighSDDValue(NULL),fCalibration(NULL) {;}          // ctor
  AliITSQASDDChecker& operator = (const AliITSQASDDChecker& qac) ; //operator =
  virtual ~AliITSQASDDChecker(); /*{if(fStepBitSDD) delete[] fStepBitSDD ;if(fLowSDDValue)delete[]fLowSDDValue;if(fHighSDDValue) delete[]fHighSDDValue;if(fCalibration)delete fCalibration;} */// dtor
  virtual Double_t Check(AliQAv1::ALITASK_t index, const TObjArray * list, const AliDetectorRecoParam * recoParam);
  virtual void SetTaskOffset(Int_t taskoffset);
  virtual void SetStepBit(const Double_t *steprange);
  virtual Double_t *GetStepBit(){return fStepBitSDD;};
  virtual void SetSDDLimits(const Float_t *lowvalue, const Float_t * highvalue);
private:
  AliITSQASDDChecker(const AliITSQASDDChecker& qac):TObject(),fSubDetOffset(qac.fSubDetOffset),fStepBitSDD(qac.fStepBitSDD),fLowSDDValue(qac.fLowSDDValue),fHighSDDValue(qac.fHighSDDValue),fCalibration(qac.fCalibration) {;} // cpy ctor   
  Int_t fSubDetOffset;            // checking operation starting point
  Double_t *fStepBitSDD;          //step size for each QAbit(kINFO, kWARNING,kERROR,kFATAL)
  Float_t *fLowSDDValue;          //low value of each QA bit range 
  Float_t *fHighSDDValue;         //High value of each QA bit range
  TObjArray *fCalibration;        //TObjArray with Calibration SDD Objects

  static const Int_t fgknSDDmodules = 260; // number of SDD modules
  static const Int_t fgkmodoffset = 240;   // number of SPD modules
  ClassDef(AliITSQASDDChecker,2)  // description 

};

#endif // AliITSQASDDChecker_H
