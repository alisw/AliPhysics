#ifndef ALISTARTRAWREADER_H
#define ALISTARTRAWREADER_H
 
#include <TTask.h>
#include <Riostream.h>
#include "TTree.h"
#include "AliSTARTdigit.h"
#include "AliRawReader.h"
 
class AliSTARTRawReader : public TTask {
  public :

  AliSTARTRawReader(AliRawReader *rawReader, TTree* tree) ;

  virtual  ~AliSTARTRawReader();


  Bool_t  Next(); //read next raw digit
  
  protected :

    UInt_t           fData;         // data read for file
  AliSTARTdigit* fDigits;
  TTree*        fTree;
  AliRawReader*    fRawReader;    // object for reading the raw data
  
 ClassDef(AliSTARTRawReader, 0) //class for reading START Raw data
};
 
#endif
