#ifndef ALIFMDCALIBFAKER_H
#define ALIFMDCALIBFAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */

//____________________________________________________________________
//
//  Class to make fake calibration parameters 
//
#ifndef ROOT_TTask
# include <TTask.h>
#endif
#include "AliFMDParameters.h"	   // ALIFMDPARAMETERS_H

class AliFMDCalibFaker : public TTask
{
public:
  enum EWhat {
    kZeroSuppression =  1, 
    kSampleRate,
    kPedestal,
    kPulseGain,
    kDeadMap,
    kAltroMap
  };
  enum {
    kAll             = (1<<kZeroSuppression|1<<kSampleRate|1<<kPedestal|
			1<<kPulseGain|1<<kDeadMap|1<<kAltroMap)
  };
  AliFMDCalibFaker(Int_t mask=kAll, const char* loc="local://cdb");
  virtual ~AliFMDCalibFaker() {}
  void AddCalib(EWhat w) { SETBIT(fMask, w); }
  void RemoveCalib(EWhat w) { SETBIT(fMask, w); }
  void SetCalib(Int_t mask) { fMask = mask; }
  void SetGainSeed(Float_t g) { fGain = g; }
  void SetThresholdFactor(Float_t t) { fThresholdFactor = t; }
  void SetPedestalRange(Float_t min, Float_t max) 
  {
    fPedestalMin = min;
    fPedestalMax = (max < min ? min : max);
  }
  void SetRunRange(Int_t min, Int_t max) 
  {
    fRunMin = min;
    fRunMax = (max < min ? min : max);
  }
  void SetDeadChance(Float_t chance) { fDeadChance = chance; }
  void SetRate(UShort_t rate) { fRate = rate; }
  void SetZeroThreshold(UShort_t t) { fZeroThreshold = t; }
  void SetDefaultStorage(const char* url) { SetTitle(url); }
  void Exec(Option_t* option="");
protected:
  virtual AliFMDCalibZeroSuppression* MakeZeroSuppression();
  virtual AliFMDCalibSampleRate*      MakeSampleRate();
  virtual AliFMDCalibPedestal*        MakePedestal();
  virtual AliFMDCalibGain*            MakePulseGain();
  virtual AliFMDCalibDeadMap*         MakeDeadMap();
  virtual AliFMDAltroMapping*         MakeAltroMap();

  Long_t   fMask;            // What to write 
  Float_t  fGain;            // Gain
  Float_t  fThresholdFactor; // Threshold factor
  Float_t  fThreshold;       // Threshold
  Float_t  fPedestalMin;     // Min pedestal
  Float_t  fPedestalMax;     // Max pedestal
  Float_t  fDeadChance;      // Chance of dead channel
  UShort_t fRate;            // Sample rate 
  UShort_t fZeroThreshold;   // Zero suppression threshold
  Int_t    fRunMin;
  Int_t    fRunMax;
  
  ClassDef(AliFMDCalibFaker,0)
};

#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
// EOF
//

