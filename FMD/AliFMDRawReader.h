#ifndef ALIFMDRAWREADER_H
#define ALIFMDRAWREADER_H
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
// Class to read ADC values from a AliRawReader object. 
// 
#ifndef ROOT_TTask
# include <TTask.h>
#endif

//____________________________________________________________________
class AliRawReader;
class TTree;


//____________________________________________________________________
class AliFMDRawReader : public TTask 
{
public:
  AliFMDRawReader(AliRawReader* reader, TTree* array);

  virtual void Exec(Option_t* option="");
protected:
  TTree*        fTree;       //! Pointer to tree to read into 
  AliRawReader* fReader;     //! Pointer to raw reader 
  UShort_t      fSampleRate; // The sample rate (if 0, inferred from data)
  ClassDef(AliFMDRawReader, 0) // Read FMD raw data into a cache 
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
