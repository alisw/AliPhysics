#ifndef ALIFMDCALIBGAIN_H
#define ALIFMDCALIBGAIN_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
#ifndef ALIFMDFLOATMAP_H
# include <AliFMDFloatMap.h>
#endif
//____________________________________________________________________
//
// Gain value and width for each strip in the FMD
//
class AliFMDCalibGain : public TObject 
{
public:
  AliFMDCalibGain();
  ~AliFMDCalibGain() {}
  AliFMDCalibGain(const AliFMDCalibGain& o);
  AliFMDCalibGain& operator=(const AliFMDCalibGain& o);
  void Set(UShort_t det, Char_t ring, UShort_t sec, UShort_t str, Float_t val);
  void Set(Float_t thres) { fThreshold = thres; }
  Float_t Value(UShort_t det, Char_t ring, UShort_t sec, UShort_t str);
  Float_t Threshold() const { return fThreshold; }
private:
  AliFMDFloatMap fValue;
  Float_t        fThreshold;
  ClassDef(AliFMDCalibGain, 1) // Gain data for the FMD 
};


#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//


