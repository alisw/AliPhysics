#ifndef ALI_SHUTTLE_INTERFACE_H
#define ALI_SHUTTLE_INTERFACE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// abstract interface class to AliShuttle
//

#include <TObject.h>

class TList;
class AliCDBMetaData;
class AliPreprocessor;
class AliCDBPath;

class AliShuttleInterface : public TObject
{
  public:
    enum { kDAQ = 0, kDCS, kHLT };
    static const char* fkSystemNames[3];  // names of the systems providing data to the shuttle

    virtual UInt_t Store(const AliCDBPath& path, TObject* object, AliCDBMetaData* metaData,
    				Int_t validityStart = 0, Bool_t validityInfinite = kFALSE) = 0;
    virtual UInt_t StoreReferenceData(const AliCDBPath& path, TObject* object, AliCDBMetaData* metaData) = 0;
    virtual const char* GetFile(Int_t system, const char* detector, const char* id, const char* source) = 0;
    virtual TList* GetFileSources(Int_t system, const char* detector, const char* id) = 0;
    virtual void Log(const char* detector, const char* message) = 0;

    virtual void RegisterPreprocessor(AliPreprocessor* preprocessor) = 0;

  private:
    ClassDef(AliShuttleInterface, 0);
};

#endif
