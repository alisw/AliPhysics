#ifndef ALIRAWREADER_H
#define ALIRAWREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This is the base class for reading raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TArrayI.h>
#include "AliRawDataHeader.h"

class AliRawEventHeaderBase;

class AliRawReader: public TObject {
  public :
    AliRawReader();
    AliRawReader(const AliRawReader& rawReader);
    AliRawReader& operator = (const AliRawReader& rawReader);
    virtual ~AliRawReader();

    void             Select(Int_t detectorID, 
			    Int_t minDDLID = -1, Int_t maxDDLID = -1);
    void             Select(const char *detectorName, 
			    Int_t minDDLID = -1, Int_t maxDDLID = -1);
    void             SelectEquipment(Int_t equipmentType, 
				     Int_t minEquipmentId = -1, 
				     Int_t maxEquipmentId = -1);
    void             SkipInvalid(Bool_t skip = kTRUE)
      {fSkipInvalid = skip;};
    void             SelectEvents(Int_t type);
    virtual void     RequireHeader(Bool_t required)
      {fRequireHeader = required;};

    virtual const AliRawEventHeaderBase* GetEventHeader() const {return NULL;};

    virtual UInt_t   GetType() const = 0;
    virtual UInt_t   GetRunNumber() const = 0;
    virtual const UInt_t* GetEventId() const = 0;
    virtual const UInt_t* GetTriggerPattern() const = 0;
    virtual const UInt_t* GetDetectorPattern() const = 0;
    virtual const UInt_t* GetAttributes() const = 0;
    virtual const UInt_t* GetSubEventAttributes() const = 0;
    virtual UInt_t   GetLDCId() const = 0;
    virtual UInt_t   GetGDCId() const = 0;

    virtual Int_t    GetEquipmentSize() const = 0;
    virtual Int_t    GetEquipmentType() const = 0;
    virtual Int_t    GetEquipmentId() const = 0;
    Int_t            GetMappedEquipmentId() const;
    Bool_t           LoadEquipmentIdsMap(const char *fileName);
    virtual const UInt_t* GetEquipmentAttributes() const = 0;
    virtual Int_t    GetEquipmentElementSize() const = 0;
    virtual Int_t    GetEquipmentHeaderSize() const = 0;

    Int_t            GetDetectorID() const;
    Int_t            GetDDLID() const;

    Int_t            GetDataSize() const 
      {if (fHeader) {
	if (fHeader->fSize != 0xFFFFFFFF) return fHeader->fSize - sizeof(AliRawDataHeader); 
	else return GetEquipmentSize() - GetEquipmentHeaderSize() - sizeof(AliRawDataHeader);
      } else return GetEquipmentSize() - GetEquipmentHeaderSize();};

    Int_t            GetVersion() const 
      {if (fHeader) return fHeader->fVersion; else return -1;};
    Bool_t           IsValid() const 
      {if (fHeader) return fHeader->TestAttribute(0); 
      else return kFALSE;};
    Bool_t           IsCompressed() const 
      {if (fHeader) return fHeader->TestAttribute(1); 
      else return kFALSE;};
    Bool_t           TestBlockAttribute(Int_t index) const
      {if (fHeader) return fHeader->TestAttribute(index); 
      else return kFALSE;};
    UChar_t          GetBlockAttributes() const 
      {if (fHeader) return fHeader->GetAttributes(); 
      else return 0;};
    UInt_t           GetStatusBits() const
      {if (fHeader) return fHeader->GetStatus(); 
      else return 0;};
    const AliRawDataHeader* GetDataHeader() const
      {return fHeader;}

    virtual Bool_t   ReadHeader() = 0;
    virtual Bool_t   ReadNextData(UChar_t*& data) = 0;
    virtual Bool_t   ReadNextInt(UInt_t& data);
    virtual Bool_t   ReadNextShort(UShort_t& data);
    virtual Bool_t   ReadNextChar(UChar_t& data);
    virtual Bool_t   ReadNext(UChar_t* data, Int_t size) = 0;

    virtual Bool_t   Reset() = 0;

    virtual Bool_t   NextEvent() = 0;
    virtual Bool_t   RewindEvents() = 0;

    enum {kErrMagic=1, kErrNoDataHeader=2, 
	  kErrSize=4, kErrOutOfBounds=8};
    virtual Int_t    CheckData() const;
    Int_t            GetErrorCode() const {return fErrorCode;};

    void             DumpData(Int_t limit = -1);

  protected :
    Bool_t           IsSelected() const;
    Bool_t           IsEventSelected() const;

    TArrayI         *fEquipmentIdsIn;       // array of equipment Ids to be mapped
    TArrayI         *fEquipmentIdsOut;      // array of mapped equipment Ids

    Bool_t           fRequireHeader;        // if false, data without header is accepted

    AliRawDataHeader* fHeader;              // current data header
    Int_t            fCount;                // counter of bytes to be read for current DDL

    Int_t            fSelectEquipmentType;  // type of selected equipment (<0 = no selection)
    Int_t            fSelectMinEquipmentId; // minimal index of selected equipment (<0 = no selection)
    Int_t            fSelectMaxEquipmentId; // maximal index of selected equipment (<0 = no selection)
    Bool_t           fSkipInvalid;          // skip invalid data
    Int_t            fSelectEventType;      // type of selected events (<0 = no selection)

    Int_t            fErrorCode;            // code of last error

    ClassDef(AliRawReader, 0) // base class for reading raw digits
};

#endif
