#ifndef ALIMUONRAWSTREAM_H
#define ALITMUONAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup rec
/// \class AliMUONRawStream
/// \brief Class for reading MUON raw digits
///
///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to MUON digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliRawReader;


class AliMUONRawStream: public TObject {
  public :
    AliMUONRawStream(AliRawReader* rawReader);
    AliMUONRawStream(const AliMUONRawStream& stream);
    AliMUONRawStream& operator = (const AliMUONRawStream& stream);
    virtual ~AliMUONRawStream();

    virtual Bool_t   Next();


  protected :

    AliRawReader*    fRawReader;    // object for reading the raw data

 
    static const Int_t fgkDataMax = 10000000; // size of array for uncompressed raw data
    UShort_t*        fData;         // uncompressed raw data
    Int_t            fDataSize;     // actual size of the uncompressed raw data
    Int_t            fPosition;     // current position in fData
    Int_t            fCount;        // counter of words to be read for current trailer


    ClassDef(AliMUONRawStream, 0)    // base class for reading MUON raw digits
};

#endif
