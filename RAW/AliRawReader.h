#ifndef ALIRAWREADER_H
#define ALIRAWREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>


struct AliMiniHeader {
  UInt_t    fSize;
  UChar_t   fDetectorID;
  UChar_t   fMagicWord[3];
  UChar_t   fVersion;
  UChar_t   fCompressionFlag;
  UShort_t  fDDLID;
};

class AliRawReader: public TObject {
  public :
    AliRawReader();

    void             Select(Int_t detectorID, 
			    Int_t minDDLID = -1, Int_t maxDDLID = -1);

    virtual UInt_t   GetType() = 0;
    virtual UInt_t   GetRunNumber() = 0;
    virtual const UInt_t* GetEventId() = 0;
    virtual const UInt_t* GetTriggerPattern() = 0;
    virtual const UInt_t* GetDetectorPattern() = 0;
    virtual const UInt_t* GetAttributes() = 0;
    virtual UInt_t   GetGDCId() = 0;

    inline Int_t     GetDataSize() const {return fMiniHeader->fSize;};
    inline Int_t     GetDetectorID() const {return fMiniHeader->fDetectorID;};
    inline Int_t     GetDDLID() const {return fMiniHeader->fDDLID;};
    inline Int_t     GetVersion() const {return fMiniHeader->fVersion;};
    inline Bool_t    IsCompressed() const {return fMiniHeader->fCompressionFlag != 0;};

    virtual Bool_t   ReadMiniHeader() = 0;
    virtual Bool_t   ReadNextData(UChar_t*& data) = 0;
    virtual Bool_t   ReadNextInt(UInt_t& data);
    virtual Bool_t   ReadNextShort(UShort_t& data);
    virtual Bool_t   ReadNextChar(UChar_t& data);

    virtual Bool_t   Reset() = 0;

  protected :
    Bool_t           IsSelected();

    Bool_t           CheckMiniHeader();
    virtual Bool_t   ReadNext(UChar_t* data, Int_t size) = 0;

    AliMiniHeader*   fMiniHeader;  // current mini header
    Int_t            fCount;       // counter of bytes to be read for current DDL

    Int_t            fSelectDetectorID;  // id of selected detector (<0 = no selection)
    Int_t            fSelectMinDDLID;    // minimal index of selected DDLs (<0 = no selection)
    Int_t            fSelectMaxDDLID;    // maximal index of selected DDLs (<0 = no selection)

    ClassDef(AliRawReader, 0) // base class for reading raw digits
};

#endif
