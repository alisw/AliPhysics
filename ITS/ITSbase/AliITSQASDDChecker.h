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
class TH1;
class TH2;
class TCanvas;
class TPaveText;
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
  Int_t GetEventSpecieForCheck() const {return  fESforCheck;}

  virtual Bool_t   MakeSDDImage( TObjArray ** list, AliQAv1::TASKINDEX_t task, AliQAv1::MODE_t mode) ; 
  Bool_t DrawHistos(TObjArray ** list, AliQAv1::TASKINDEX_t task, AliQAv1::MODE_t mode );

  enum {kActiveMod,kFilledMod,kActiveWing,kFilledWing,
	kExcludedMod,kEmptyMod,kExcludedWing,kEmptyWing,
	kExcludedButFilledMod,kActiveButEmptyMod,
	kExcludedButFilledWing,kActiveButEmptyWing,
	kNumOfSDDCheckerCounters};

 private:
  
  void FillCounters(Int_t counters[kNumOfSDDCheckerCounters][3], TH1* hmodule, TH2* hlay3, TH2* hlay4);
  Int_t CheckCounters(Int_t counters[kNumOfSDDCheckerCounters][3], Int_t jl, Int_t neventsraw, TPaveText *ptext);

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
    TPaveText *    fPaveText[AliRecoParam::kNSpecies]         ; //[AliRecoParam::kNSpecies] 

    Int_t fESforCheck; //eventspecie of the list to check

    static const Int_t fgknSDDmodules = 260; // number of SDD modules
    static const Int_t fgkmodoffset = 240;   // number of SPD modules

    ClassDef(AliITSQASDDChecker,8)  // description 
      
};

#endif // AliITSQASDDChecker_H
