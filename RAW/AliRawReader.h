#ifndef ALIRAWREADER_H
#define ALIRAWREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include "AliMiniHeader.h"


class AliRawReader: public TObject {
  public :
    AliRawReader();
    AliRawReader(const AliRawReader& rawReader);
    AliRawReader& operator = (const AliRawReader& rawReader);
    virtual ~AliRawReader() {};

    void             Select(Int_t detectorID, 
			    Int_t minDDLID = -1, Int_t maxDDLID = -1);

    virtual UInt_t   GetType() const = 0;
    virtual UInt_t   GetRunNumber() const = 0;
    virtual const UInt_t* GetEventId() const = 0;
    virtual const UInt_t* GetTriggerPattern() const = 0;
    virtual const UInt_t* GetDetectorPattern() const = 0;
    virtual const UInt_t* GetAttributes() const = 0;
    virtual UInt_t   GetGDCId() const = 0;

    Int_t            GetDataSize() const {return fMiniHeader->fSize;};
    Int_t            GetDetectorID() const {return fMiniHeader->fDetectorID;};
    Int_t            GetDDLID() const {return fMiniHeader->fDDLID;};
    Int_t            GetVersion() const {return fMiniHeader->fVersion;};
    Bool_t           IsCompressed() const {return fMiniHeader->fCompressionFlag != 0;};

    virtual Bool_t   ReadMiniHeader() = 0;
    virtual Bool_t   ReadNextData(UChar_t*& data) = 0;
    virtual Bool_t   ReadNextInt(UInt_t& data);
    virtual Bool_t   ReadNextShort(UShort_t& data);
    virtual Bool_t   ReadNextChar(UChar_t& data);

    virtual Bool_t   Reset() = 0;

    enum {kErrMagic=1, kErrNoMiniHeader=2, kErrMiniMagic=3, 
	  kErrSize=4, kErrOutOfBounds=5};
    virtual Int_t    CheckData() const;
    Int_t            GetErrorCode() {return fErrorCode;};

  protected :
    Bool_t           IsSelected() const;

    Bool_t           CheckMiniHeader(AliMiniHeader* miniHeader = NULL) const;
    virtual Bool_t   ReadNext(UChar_t* data, Int_t size) = 0;

    AliMiniHeader*   fMiniHeader;  // current mini header
    Int_t            fCount;       // counter of bytes to be read for current DDL

    Int_t            fSelectDetectorID;  // id of selected detector (<0 = no selection)
    Int_t            fSelectMinDDLID;    // minimal index of selected DDLs (<0 = no selection)
    Int_t            fSelectMaxDDLID;    // maximal index of selected DDLs (<0 = no selection)

    Int_t            fErrorCode;         // code of last error

    ClassDef(AliRawReader, 0) // base class for reading raw digits
};

#endif
