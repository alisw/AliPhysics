#ifndef ALIT0RAWREADER_H
#define ALIT0RAWREADER_H
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
#include "AliT0Parameters.h"

  
class AliT0RawReader : public TTask {
  public :

  AliT0RawReader(AliRawReader *rawReader,Bool_t isOnline = kFALSE) ;
//  AliT0RawReader(AliRawReader *rawReader,Bool_t isOnline = kTRUE) ;

  virtual  ~AliT0RawReader();
  AliT0RawReader(const AliT0RawReader& o): TTask(o),
    fRawReader(0),
    fData(NULL),
    fPosition(0),
    fParam(0),
    fIsOnline(kFALSE),
    fBunchID(0),
    fPrintout(kFALSE)
    {}
  
  AliT0RawReader& operator=(const AliT0RawReader&) { return *this; }


  Bool_t  Next(); //read next raw digit
  Int_t            GetPosition();
  UInt_t         GetNextWord();
  Int_t GetData( Int_t channel, Int_t hit) {return fAllData[channel][hit];}


  enum ET0RawReaderError {
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

  Bool_t IsOnlineMode () {return fIsOnline;}
  Int_t GetTRMBunchID() {return fBunchID;};

  void SetPrintout(Bool_t pp ) {fPrintout = pp;}
     
  protected :

  AliRawReader*    fRawReader;    // object for reading the raw data
  UChar_t*         fData;         // raw data
  Int_t            fPosition;     // current (32 bit) position in fData
  AliT0Parameters *fParam;       // instanse of  Parameters class
  Bool_t           fIsOnline;     // for case online DA usage
  Int_t            fBunchID;       //bunchID from TRM chain header
  Bool_t           fPrintout;      // advanced printout
  Int_t            fAllData[250][5]; // container for raw data
  
 ClassDef(AliT0RawReader,5) //class for reading T0 Raw data
};

typedef AliT0RawReader AliSTARTRawReader; // for backward compatibility
 
#endif
