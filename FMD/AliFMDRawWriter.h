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


//____________________________________________________________________
class AliFMDRawWriter : public TTask 
{
public:
  AliFMDRawWriter(AliFMD* fmd);

  virtual void Exec(Option_t* option="");
  void SetSampleRate(UShort_t sampleRate=1) { fSampleRate = sampleRate; }
  void SetChannelsPerAltro(UShort_t size=128) { fChannelsPerAltro = size; }
protected:
  AliFMD*       fFMD;        //! Pointer to detector description 
  UShort_t      fSampleRate; // The sample rate (if 0, inferred from data)
  UShort_t      fChannelsPerAltro;
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
