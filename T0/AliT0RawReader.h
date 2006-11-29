#ifndef ALIT0RAWREADER_H
#define ALIT0RAWREADER_H
 
#include <TTask.h>
#include <Riostream.h>
#include "TTree.h"
#include "AliT0digit.h"
#include "AliRawReader.h"
 
class AliT0RawReader : public TTask {
  public :

  AliT0RawReader(AliRawReader *rawReader, TTree* tree) ;

  virtual  ~AliT0RawReader();


  Bool_t  Next(); //read next raw digit
  Int_t            GetPosition();
  // void UnpackTime(Int_t outTime, Int_t outCh);
  UInt_t         GetNextWord();
  
  protected :

  AliT0digit* fDigits;
  TTree*        fTree;
  AliRawReader*    fRawReader;    // object for reading the raw data

  UChar_t*         fData;         // raw data
  Int_t            fPosition;     // current (32 bit) position in fData

  
 ClassDef(AliT0RawReader, 0) //class for reading T0 Raw data
};

typedef AliT0RawReader AliSTARTRawReader; // for backward compatibility
 
#endif
