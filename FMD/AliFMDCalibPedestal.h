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
/** Pedestal value and width for each strip in the FMD 
    @ingroup FMD_base
*/
class AliFMDCalibPedestal : public TObject 
{
public:
  /** CTOR */
  AliFMDCalibPedestal();
  /** DTOR */
  ~AliFMDCalibPedestal() {}
  /** Copy ctor 
      @param o Object to copy from  */
  AliFMDCalibPedestal(const AliFMDCalibPedestal& o);
  /** Assignment 
      @param o Object to assign from
      @return Reference to this object   */
  AliFMDCalibPedestal& operator=(const AliFMDCalibPedestal& o);
  /** Set the values for a strip. 
      @param det  Detector 
      @param ring Ring 
      @param sec  Sector 
      @param str  Strip
      @param ped  Value of pedestal 
      @param pedW Width of pedestal */
  void Set(UShort_t det, Char_t ring, UShort_t sec, UShort_t str, 
	   Float_t ped, Float_t pedW);
  /** Get pedestal for a strip. 
      @param det  Detector 
      @param ring Ring 
      @param sec  Sector 
      @param str  Strip
      @param val  Value of gain 
      @return Pedestal for strip */  
  Float_t Value(UShort_t det, Char_t ring, UShort_t sec, UShort_t str);
  /** Get pedestal width for a strip. 
      @param det  Detector 
      @param ring Ring 
      @param sec  Sector 
      @param str  Strip
      @param val  Value of gain 
      @return Pedestal width for strip */  
  Float_t Width(UShort_t det, Char_t ring, UShort_t sec, UShort_t str);
private:
  AliFMDFloatMap fValue; /** Pedestal */
  AliFMDFloatMap fWidth; /** Pedestal width */
  ClassDef(AliFMDCalibPedestal, 1) // Pedestal data for the FMD 
};


#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//


