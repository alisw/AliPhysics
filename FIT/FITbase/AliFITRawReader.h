#ifndef ALIFITRAWREADER_H
#define ALIFITRAWREADER_H
/***************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Alla Maevskaya INR RAS alla@inr.ru
 *
 * See cxx source for full Copyright notice                               
 ***************************************************************************/

 
#include <TTask.h>
#include "AliRawReader.h"

  
class AliFITRawReader : public TTask {
  public :

  AliFITRawReader(AliRawReader *rawReader) ;

  virtual  ~AliFITRawReader();
  AliFITRawReader(const AliFITRawReader& o): TTask(o),
    fRawReader(0),
    fData(NULL),
    fPosition(0),
    fBunchID(0),
    fPrintout(kFALSE)
    { for ( Int_t k=0; k<1000; k++)   fAllData[k] = -1;}
  
  AliFITRawReader& operator=(const AliFITRawReader&) { return *this; }


  Bool_t  Next(); //read next raw digit
  Int_t            GetPosition();
  UInt_t         GetNextWord();
  Int_t GetData( Int_t channel) {return fAllData[channel];}


  enum EFITRawReaderError {
    kIncorrectDataSize = 1,
    kWrongDRMHeader = 2,
    kWrongDRMTrailer = 3,
    kWrongTRMHeader = 4,
    kWrongTRMTrailer = 5,
    kWrongChain0Header = 6,
    kWrongChain0Trailer = 7,
    kWrongChain1Header = 8,
    kWrongChain1Trailer = 9,
    kIncorrectLUT = 10
  };

   Int_t GetTRMBunchID() {return fBunchID;};

  void SetPrintout(Bool_t pp ) {fPrintout = pp;}
  UInt_t GetChannel(Int_t iTRM, Int_t iTDC, Int_t iChain, Int_t ichannel);     
  protected :

  AliRawReader*    fRawReader;    // object for reading the raw data
  UChar_t*         fData;         // raw data
  Int_t            fPosition;     // current (32 bit) position in fData
  Int_t            fBunchID;       //bunchID from TRM chain header
  Bool_t           fPrintout;      // advanced printout
  Int_t            fAllData[1000]; // container for raw data
  
 ClassDef(AliFITRawReader,1) //class for reading FIT Raw data
};

#endif
