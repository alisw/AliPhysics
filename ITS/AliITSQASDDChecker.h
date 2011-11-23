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
//#include "TPad.h"
// --- ROOT system ---
class TFile;
class TH2F; 

// --- AliRoot header files ---

class AliQACheckerBase;
class AliQAChecker;
class AliITSQAChecker;
class AliITSCalibrationSDD;
class AliITSLoader;
class TSystem;
class AliQAManager;
class AliLog;
class TF1;
class TCanvas;
class AliQAv1;

class AliITSQASDDChecker: public TObject{

public:

  AliITSQASDDChecker();  // ctor
  virtual ~AliITSQASDDChecker(); // dtor
  virtual Double_t Check(AliQAv1::ALITASK_t index, const TObjArray * list, const AliDetectorRecoParam * recoParam);
  virtual void SetTaskOffset(Int_t taskoffset);
  virtual void SetStepBit(const Double_t *steprange);
  virtual Double_t *GetStepBit(){return fStepBitSDD;};
  virtual void SetSDDLimits(const Float_t *lowvalue, const Float_t * highvalue);
  void SetEventSpecieForCheck(Int_t esforcheck=0){ fESforCheck=esforcheck;}
  Int_t GetEventSpecieForCheck(){return  fESforCheck;}

  virtual Bool_t   MakeSDDImage( TObjArray ** list, AliQAv1::TASKINDEX_t task, AliQAv1::MODE_t mode) ; 
  Bool_t MakeSDDRawsImage(TObjArray ** list, AliQAv1::TASKINDEX_t task, AliQAv1::MODE_t mode );//{AliInfo("The method for raw image has been called\n");}
  Bool_t MakeSDDRecPointsImage(TObjArray ** list, AliQAv1::TASKINDEX_t task, AliQAv1::MODE_t mode);//{AliInfo("The method for recpoint image has been called\n");}


 private:
  
  AliITSQASDDChecker(const AliITSQASDDChecker& qac); 
  AliITSQASDDChecker& operator = (const AliITSQASDDChecker& qac) ; 
    Int_t fSubDetOffset;            // checking operation starting point
    Double_t *fStepBitSDD;          //step size for each QAbit(kINFO, kWARNING,kERROR,kFATAL)
    Float_t *fLowSDDValue;          //low value of each QA bit range 
    Float_t *fHighSDDValue;         //High value of each QA bit range
    TObjArray *fCalibration;        //TObjArray with Calibration SDD Objects
        
    
    Float_t fThresholdForRelativeOccupancy;  // ThresholdForRelativeOccupancy (by module)
    Float_t fThresholdForRecToRawRatio; // ThresholdForRecToRawRatio (by module)
    
    TCanvas **    fImage          ; //[AliRecoParam::kNSpecies] 

    Int_t fESforCheck; //eventspecie of the list to check

    static const Int_t fgknSDDmodules = 260; // number of SDD modules
    static const Int_t fgkmodoffset = 240;   // number of SPD modules

    ClassDef(AliITSQASDDChecker,6)  // description 
      
};

#endif // AliITSQASDDChecker_H
