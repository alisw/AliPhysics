#ifndef ALIFMDCALIBPEDESTAL_H
#define ALIFMDCALIBPEDESTAL_H
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
// Pedestal value and width for each strip in the FMD
//
class AliFMDCalibPedestal : public TObject 
{
public:
  AliFMDCalibPedestal();
  ~AliFMDCalibPedestal() {}
  AliFMDCalibPedestal(const AliFMDCalibPedestal& o);
  AliFMDCalibPedestal& operator=(const AliFMDCalibPedestal& o);
  void Set(UShort_t det, Char_t ring, UShort_t sec, UShort_t str, 
	   Float_t ped, Float_t pedW);
  Float_t Value(UShort_t det, Char_t ring, UShort_t sec, UShort_t str);
  Float_t Width(UShort_t det, Char_t ring, UShort_t sec, UShort_t str);
private:
  AliFMDFloatMap fValue;
  AliFMDFloatMap fWidth;
  ClassDef(AliFMDCalibPedestal, 1) // Pedestal data for the FMD 
};


#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//


