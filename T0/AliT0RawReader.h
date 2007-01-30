#ifndef ALIT0RAWREADER_H
#define ALIT0RAWREADER_H
 
#include <TTask.h>
#include <Riostream.h>
#include "AliRawReader.h"
 
class AliT0RawReader : public TTask {
  public :

  AliT0RawReader(AliRawReader *rawReader) ;

  virtual  ~AliT0RawReader();
  AliT0RawReader(const AliT0RawReader& o)
    : TTask(o)
  {}
  AliT0RawReader& operator=(const AliT0RawReader&) { return *this; }


  Bool_t  Next(); //read next raw digit
  Int_t            GetPosition();
  UInt_t         GetNextWord();
  Int_t GetData(Int_t channel, Int_t hit) {return fAllData[channel][hit];}
  
  protected :

  AliRawReader*    fRawReader;    // object for reading the raw data

  UChar_t*         fData;         // raw data
  Int_t            fPosition;     // current (32 bit) position in fData
  Int_t fAllData[110][5];
  
 ClassDef(AliT0RawReader, 0) //class for reading T0 Raw data
};

typedef AliT0RawReader AliSTARTRawReader; // for backward compatibility
 
#endif
