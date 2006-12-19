#ifndef ALIFMDCALIBFAKER_H
#define ALIFMDCALIBFAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
/** @file    AliFMDCalibFaker.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 18:29:36 2006
    @brief   Make fake calibration data 
    @ingroup FMD_util
*/
//____________________________________________________________________
//
//  Class to make fake calibration parameters 
// Used for populating inital data base
// Utility class 
//
#ifndef ROOT_TTask
# include <TTask.h>
#endif
#include "AliFMDParameters.h"	   // ALIFMDPARAMETERS_H

/** @class AliFMDCalibFaker
    @brief Class to make fake calibration parameters. 
    @code 
    AliFMDCalibFaker f(0);
    f.AddCalib(kPedestal);
    f.AddCalib(kGain);
    f.Exec();
    @endcode 
    @ingroup FMD_util
 */
class AliFMDCalibFaker : public TTask
{
public:
  /** What to make */
  enum EWhat {
    /** Zero suppressio threshold */
    kZeroSuppression =  1, 
    /** Sample rate */
    kSampleRate,
    /** Pedestals */
    kPedestal,
    /** Gain */
    kPulseGain,
    /** Dead map */
    kDeadMap,
    /** Hardware map */
    kAltroMap, 
    /** Strip range */
    kStripRange
  };
  enum {
    /** All parameters */
    kAll             = (1<<kZeroSuppression|1<<kSampleRate|1<<kPedestal|
			1<<kPulseGain|1<<kDeadMap|1<<kAltroMap|1<<kStripRange)
  };
  /** Constructor 
      @param mask Bit mask of what to make 
      @param loc  Where to store the results */
  AliFMDCalibFaker(Int_t mask=kAll, const char* loc="local://cdb");
  /** Destructor */
  virtual ~AliFMDCalibFaker() {}
  /** Add a parameter to output 
      @param w Bit of parameter */
  void AddCalib(EWhat w) { SETBIT(fMask, w); }
  /** Remove a parameter from output 
      @param w Bit of parameter */
  void RemoveCalib(EWhat w) { SETBIT(fMask, w); }
  /** Set the bit mask of what to make 
      @param mask bit mask */
  void SetCalib(Int_t mask) { fMask = mask; }
  /** Set seed for random gain.  The gain is distributed flatly from
      90 to 110 percent of the seed. 
      @param g Seed for gain */
  void SetGainSeed(Float_t g) { fGain = g; }
  /** Set the threshold factor.  The actual threshold is the gain seed
      times the factor 
      @param t Factor */
  void SetThresholdFactor(Float_t t) { fThresholdFactor = t; }
  /** Set the limits for the random pedestal.  The pedestal values are
      randomly distributed in the range 
      @param min Minimum of range 
      @param max Maximum of range */
  void SetPedestalRange(Float_t min, Float_t max) 
  {
    fPedestalMin = min;
    fPedestalMax = (max < min ? min : max);
  }
  /** Set run validty range 
      @param min Minimum run number
      @param max Maximum run number */
  void SetRunRange(Int_t min, Int_t max) 
  {
    fRunMin = min;
    fRunMax = (max < min ? min : max);
  }
  /** Set the likelyness that a strip is dead. 
      @param chance Chance of dead channel. */
  void SetDeadChance(Float_t chance) { fDeadChance = chance; }
  /** Set Sample rate 
      @param rate Rate */
  void SetRate(UShort_t rate) { fRate = rate; }
  /** Set the zero suppression threshold 
      @param t Threshold (in ADC counts) */
  void SetZeroThreshold(UShort_t t) { fZeroThreshold = t; }
  /** Set strip validty range 
      @param min Minimum strip number
      @param max Maximum strip number */
  void SetStripRange(UShort_t min, UShort_t max) 
  {
    fStripMin = min;
    fStripMax = (max < min ? min : max);
  }
  /** Set the default output storage.  It must be a CDB URL. 
      @param url CDB URL. */
  void SetDefaultStorage(const char* url) { SetTitle(url); }
  /** Make the fake calibration parameters 
      @param option Not used */
  void Exec(Option_t* option="");
protected:
  /** Make zero suppression parameters 
      @return Map of zero suppression */
  virtual AliFMDCalibZeroSuppression* MakeZeroSuppression() const;
  /** Make sample rate parameters 
      @return Map of sample rate */
  virtual AliFMDCalibSampleRate*      MakeSampleRate() const;
  /** Make pedestal parameters 
      @return Map of pedestal */
  virtual AliFMDCalibPedestal*        MakePedestal() const;
  /** Make gain parameters 
      @return Map of gain */
  virtual AliFMDCalibGain*            MakePulseGain() const;
  /** Make dead channel parameters 
      @return Map of dead channel */
  virtual AliFMDCalibDeadMap*         MakeDeadMap() const;
  /** Make a hardware map
      @return hardware map */
  virtual AliFMDAltroMapping*         MakeAltroMap() const;
  /** Make a strip range
      @return strip range map */
  virtual AliFMDCalibStripRange*      MakeStripRange() const;
  virtual Float_t MakeNoise(Char_t ring, UShort_t str) const;
  
  Long_t   fMask;            // What to write 
  Float_t  fGain;            // Gain
  Float_t  fThresholdFactor; // Threshold factor
  Float_t  fThreshold;       // Threshold
  Float_t  fPedestalMin;     // Min pedestal
  Float_t  fPedestalMax;     // Max pedestal
  Float_t  fDeadChance;      // Chance of dead channel
  UShort_t fRate;            // Sample rate 
  UShort_t fZeroThreshold;   // Zero suppression threshold
  Int_t    fRunMin;          // Min run number
  Int_t    fRunMax;          // Max run number
  UShort_t fStripMin;        // Min strip read out
  UShort_t fStripMax;        // Max strip read out
  
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

