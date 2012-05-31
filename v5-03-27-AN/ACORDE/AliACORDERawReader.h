#ifndef ALIACORDERAWREADER_H
#define ALIACORDERAWREADER_H
/***************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 ***************************************************************************/
// Authors:
//	Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch>
//	Arturo Fernandez Tellez   <afernan@mail.cern.ch>
 
#include <TNamed.h>
#include "AliRawReader.h"

  
class AliACORDERawReader : public TNamed {
  public :

  AliACORDERawReader(AliRawReader *rawReader,Bool_t isOnline = kFALSE) ;
//AliACORDERawReader(AliRawReader *rawReader,Bool_t isOnline = kTRUE) ;

  virtual  ~AliACORDERawReader();
  AliACORDERawReader(const AliACORDERawReader& o): TNamed(o),fRawReader(0),fData(NULL),fPosition(0),fIsOnline(kFALSE),fDataSize(0)
{
	fWord[0] = fWord[1] = fWord[2] = fWord[3] = 0;
}
  
  AliACORDERawReader& operator=(const AliACORDERawReader&) { return *this; }


  Bool_t  Next(); //read next raw digit
  Int_t            GetPosition();
  UInt_t         GetNextWord();
  Int_t GetData( Int_t channel, Int_t hit) {return fAllData[channel][hit];}


  enum EACORDERawReaderError {
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
     
  protected :

  AliRawReader*    fRawReader;    // object for reading the raw data
  UChar_t*         fData;         // raw data
  Int_t            fPosition;     // current (32 bit) position in fData
  Bool_t           fIsOnline;     // for case online DA usage
  UInt_t          fWord[4];      // data vector
  Int_t           fDataSize;     // data size
  Int_t            fAllData[110][50]; // container for raw data

enum EACORDERawStreamError {
      kRawDataSizeErr = 1
  };
  
 ClassDef(AliACORDERawReader,3) //class for reading ACORDE Raw data
};

typedef AliACORDERawReader AliSTARTRawReader; // for backward compatibility
 
#endif
