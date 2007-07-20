#ifndef ALIHMPIDRAWSTREAM_H
#define ALIHMPIDRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
///
/// This is a class for reading raw data digits for HMPID.
/// The data format is taken from the document provided by Paolo Martinengo.
///
/// cvetan.cheshkov@cern.ch 19/07/2007
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliRawReader;

class AliHMPIDRawStream: public TObject {
  public :
    AliHMPIDRawStream(AliRawReader* rawReader);
    virtual ~AliHMPIDRawStream();

    virtual void             Reset();
    virtual Bool_t           Next();
    void                     Init();

    Int_t GetDDLNumber()  const { return fDDLNumber; }  // Provide current DDL number
    Short_t GetCharge(Int_t row, Int_t dilogic, Int_t pad) const; // Provide the charge observed in certain row,dilogic,pad channel

      enum EHMPIDRawStreamError {
      kRawDataSizeErr = 1,
      kRowMarkerErr = 2,
      kWrongRowErr = 3,
      kWrongDilogicErr = 4,
      kWrongPadErr = 5,
      kEoEFlagErr = 6,
      kEoESizeErr = 7,
      kEoEDILOGICErr = 8,
      kEoERowErr = 9,
      kBadSegWordErr = 10,
      kWrongSegErr = 11
    };

    enum {
      kNRows       = 25, // Number of rows (starting from 1 !)
      kNDILOGICAdd = 11, // Number of DILOGIC addresses in a row (starting from 1 !)
      kNPadAdd     = 48, // Number of
      kNRowsPerSegment = 8 // Number of rows per segment
    };

  private :

    AliHMPIDRawStream& operator = (const AliHMPIDRawStream& stream);
    AliHMPIDRawStream(const AliHMPIDRawStream& stream);

    UInt_t           GetNextWord();

    Short_t          fCharge[kNRows][kNDILOGICAdd][kNPadAdd]; // Array for charge values for all channels in one DDL

    Int_t            fDDLNumber;    // index of current DDL number

    AliRawReader*    fRawReader;    // object for reading the raw data

    UChar_t*         fData;         // raw data

    Int_t            fPosition;     // current position in fData
    //    Int_t            fCount;        // counter of words to be read for current trailer

    ClassDef(AliHMPIDRawStream, 0)  // base class for reading HMPID raw digits
};

#endif
