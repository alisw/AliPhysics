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
/** Gain value and width for each strip in the FMD
    @ingroup FMD_base
*/
class AliFMDCalibSampleRate : public TObject
{
public:
  /** CTOR */
  AliFMDCalibSampleRate();
  /** Copy CTOR
      @param o Object to copy from  */
  AliFMDCalibSampleRate(const AliFMDCalibSampleRate& o);
  /** Assignment operator 
      @param o Object to assign from 
      @return Reference to assign from  */
  AliFMDCalibSampleRate& operator=(const AliFMDCalibSampleRate& o);
  /** Set sample for a DDL
      @param ddl   DDL (detector)
      @param rate  Sample rate */
  void Set(UShort_t ddl, UShort_t rate);
  /** Get sample rate for a detector 
      @param ddl Detector (DDL) identifier
      @return Sample rate */
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


