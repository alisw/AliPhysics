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
    void             SelectEquipment(Int_t equipmentType, 
				     Int_t minEquipmentId = -1, 
				     Int_t maxEquipmentId = -1);

    virtual UInt_t   GetType() const = 0;
    virtual UInt_t   GetRunNumber() const = 0;
    virtual const UInt_t* GetEventId() const = 0;
    virtual const UInt_t* GetTriggerPattern() const = 0;
    virtual const UInt_t* GetDetectorPattern() const = 0;
    virtual const UInt_t* GetAttributes() const = 0;
    virtual UInt_t   GetLDCId() const = 0;
    virtual UInt_t   GetGDCId() const = 0;

    virtual Int_t    GetEquipmentSize() const = 0;
    virtual Int_t    GetEquipmentType() const = 0;
    virtual Int_t    GetEquipmentId() const = 0;
    virtual const UInt_t* GetEquipmentAttributes() const = 0;
    virtual Int_t    GetEquipmentElementSize() const = 0;

    Int_t            GetDataSize() const 
      {if (fMiniHeader) return fMiniHeader->fSize; 
      else return GetEquipmentSize();};

    Int_t            GetDetectorID() const 
      {if (fMiniHeader) return fMiniHeader->fDetectorID; else return -1;};
    Int_t            GetDDLID() const 
      {if (fMiniHeader) return fMiniHeader->fDDLID; else return -1;};
    Int_t            GetVersion() const 
      {if (fMiniHeader) return fMiniHeader->fVersion; else return -1;};
    Bool_t           IsCompressed() const 
      {if (fMiniHeader) return fMiniHeader->fCompressionFlag != 0; 
      else return kFALSE;};

    virtual Bool_t   ReadHeader() = 0;
    virtual Bool_t   ReadNextData(UChar_t*& data) = 0;
    virtual Bool_t   ReadNextInt(UInt_t& data);
    virtual Bool_t   ReadNextShort(UShort_t& data);
    virtual Bool_t   ReadNextChar(UChar_t& data);

    virtual Bool_t   Reset() = 0;

    enum {kErrMagic=1, kErrNoMiniHeader=2, kErrMiniMagic=4, 
	  kErrSize=8, kErrOutOfBounds=16};
    virtual Int_t    CheckData() const;
    Int_t            GetErrorCode() {return fErrorCode;};

    void             DumpData(Int_t limit = -1);

  protected :
    Bool_t           IsSelected() const;

    Bool_t           CheckMiniHeader(AliMiniHeader* miniHeader = NULL) const;
    virtual Bool_t   ReadNext(UChar_t* data, Int_t size) = 0;

    AliMiniHeader*   fMiniHeader;  // current mini header
    Int_t            fCount;       // counter of bytes to be read for current DDL

    Int_t            fSelectDetectorID;  // id of selected detector (<0 = no selection)
    Int_t            fSelectMinDDLID;    // minimal index of selected DDLs (<0 = no selection)
    Int_t            fSelectMaxDDLID;    // maximal index of selected DDLs (<0 = no selection)
    Int_t            fSelectEquipmentType;  // type of selected equipment (<0 = no selection)
    Int_t            fSelectMinEquipmentId; // minimal index of selected equipment (<0 = no selection)
    Int_t            fSelectMaxEquipmentId; // maximal index of selected equipment (<0 = no selection)

    Int_t            fErrorCode;         // code of last error

    ClassDef(AliRawReader, 0) // base class for reading raw digits
};

#endif
