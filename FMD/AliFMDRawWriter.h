#ifndef ALIFMDRAWWRITER_H
#define ALIFMDRAWWRITER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
/* $Id$ */
//____________________________________________________________________
// 
// Class to writer ADC values to a Raw File
//
#ifndef ROOT_TTask
# include <TTask.h>
#endif

//____________________________________________________________________
class AliFMD;
class AliAltroBuffer;
class TArrayI;


//____________________________________________________________________
class AliFMDRawWriter : public TTask 
{
public:
  AliFMDRawWriter(AliFMD* fmd);

  virtual void Exec(Option_t* option="");
  void SetSampleRate(UShort_t sampleRate=1) { fSampleRate = sampleRate; }
  void SetChannelsPerAltro(UShort_t size=128) { fChannelsPerAltro = size; }
  void SetThreshold(UShort_t t=0) { fThreshold = t; }
protected:
  virtual void WriteChannel(AliAltroBuffer* altro, 
			    UShort_t strip, UShort_t sector, Char_t ring, 
			    const TArrayI& data);
  AliFMD*       fFMD;              //! Pointer to detector description 
  UShort_t      fSampleRate;       // The sample rate (0 -> inferred from data)
  UShort_t      fChannelsPerAltro; // Number of pre-amp. channels/adc channel 
  UShort_t      fThreshold;        // Threshold for zero-suppression
  
  ClassDef(AliFMDRawWriter, 0) // Write FMD raw data to a DDL file
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
