#ifndef ALIFMDCALIBSAMPLERATE_H
#define ALIFMDCALIBSAMPLERATE_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
#ifndef ROOT_TObject
# include <TObject.h>
#endif
#ifndef ROOT_TArrayI
# include <TArrayI.h>
#endif
//____________________________________________________________________
//
// Gain value and width for each strip in the FMD
//
class AliFMDCalibSampleRate : public TObject
{
public:
  AliFMDCalibSampleRate();
  AliFMDCalibSampleRate(const AliFMDCalibSampleRate& o);
  AliFMDCalibSampleRate& operator=(const AliFMDCalibSampleRate& o);
  void Set(UShort_t ddl, UShort_t rate);
  UShort_t Rate(UShort_t ddl) const;
protected:
  TArrayI fRates; // Sample rates 
  ClassDef(AliFMDCalibSampleRate,1); // Sample rates 
};

#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//


