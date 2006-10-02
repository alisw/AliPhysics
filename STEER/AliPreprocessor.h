#ifndef ALI_PREPROCESSOR_H
#define ALI_PREPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// This class is the CDBPreProcessor interface,
// supposed to be implemented by any detector
// interested in immediate processing of data 
// which is retrieved from DCS, DAQ or HLT.
//

#include <TNamed.h>

class TList;
class TMap;

class AliCDBMetaData;
class AliShuttleInterface;

class AliPreprocessor : public TNamed
{
  public:
    enum { kDAQ, kDCS, kHLT };

    AliPreprocessor(const char* detector, AliShuttleInterface* shuttle);
    virtual ~AliPreprocessor();

    virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual UInt_t Process(TMap* dcsAliasMap) = 0;

  protected:
    UInt_t Store(const char* pathLevel2, const char* pathLevel3, TObject* object,
    		AliCDBMetaData* metaData, Int_t validityStart = 0, Bool_t validityInfinite = kFALSE);
    UInt_t StoreReferenceData(const char* pathLevel2, const char* pathLevel3, TObject* object,
    		AliCDBMetaData* metaData);
    const char* GetFile(Int_t system, const char* id, const char* source);
    TList* GetFileSources(Int_t system, const char* id);
    void Log(const char* message);

    Int_t fRun;         // current run
    UInt_t fStartTime;  // starttime of current run
    UInt_t fEndTime;    // endtime of current run

  private:
    AliPreprocessor(const AliPreprocessor & source);
    AliPreprocessor & operator=(const AliPreprocessor & source);
    AliShuttleInterface* fShuttle;   // link to Shuttle

    ClassDef(AliPreprocessor, 0);
};

#endif
