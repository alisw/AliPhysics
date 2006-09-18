#ifndef ALIRAWREADERFILE_H
#define ALIRAWREADERFILE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This is a class for reading raw data files.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliRawReader.h"
#ifdef __CINT__
class fstream;
#else
#include <Riostream.h>
#endif
#include <TString.h>


class AliRawReaderFile: public AliRawReader {
  public :
    AliRawReaderFile(Int_t eventNumber = -1);
    AliRawReaderFile(const char* dirName, Int_t eventNumber = -1);
    virtual ~AliRawReaderFile();

    virtual void     RequireHeader(Bool_t required);

    virtual UInt_t   GetType() const {return 0;};
    virtual UInt_t   GetRunNumber() const {return 0;};
    virtual const UInt_t* GetEventId() const {return 0;};
    virtual const UInt_t* GetTriggerPattern() const {return 0;};
    virtual const UInt_t* GetDetectorPattern() const {return 0;};
    virtual const UInt_t* GetAttributes() const {return 0;};
    virtual const UInt_t* GetSubEventAttributes() const {return 0;};
    virtual UInt_t   GetLDCId() const {return 0;};
    virtual UInt_t   GetGDCId() const {return 0;};

    virtual Int_t    GetEquipmentSize() const {return 0;};
    virtual Int_t    GetEquipmentType() const {return 0;};
    virtual Int_t    GetEquipmentId() const {return fEquipmentId;};
    virtual const UInt_t* GetEquipmentAttributes() const {return NULL;};
    virtual Int_t    GetEquipmentElementSize() const {return 0;};
    virtual Int_t    GetEquipmentHeaderSize() const {return 0;};

    virtual Bool_t   ReadHeader();
    virtual Bool_t   ReadNextData(UChar_t*& data);
    virtual Bool_t   ReadNext(UChar_t* data, Int_t size);

    virtual Bool_t   Reset();

    virtual Bool_t   NextEvent();
    virtual Bool_t   RewindEvents();

  protected :
    TString          GetDirName() const;
    void*            OpenDirectory();
    Bool_t           OpenNextFile();

    Int_t            fEventIndex;  // index of the event
    TString          fDirName;     // name of the input directory
    void*            fDirectory;   // pointer to the input directory
    fstream*         fStream;      // stream of raw digits
    Int_t            fEquipmentId; // equipment ID from file name
    UChar_t*         fBuffer;      // buffer for payload
    Int_t            fBufferSize;  // size of fBuffer in bytes

  private :
    AliRawReaderFile(const AliRawReaderFile& rawReader);
    AliRawReaderFile& operator = (const AliRawReaderFile& rawReader);

    ClassDef(AliRawReaderFile, 0) // class for reading raw digits from a file
};

#endif
