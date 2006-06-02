#ifndef ALI_TEST_SHUTTLE_H
#define ALI_TEST_SHUTTLE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// test implementation of the AliShuttleInterface, to be used for local tests of preprocessors
//

#include <AliShuttleInterface.h>

class TMap;
class TList;
class AliCDBMetaData;

class AliTestShuttle : public AliShuttleInterface
{
  public:
    AliTestShuttle(TMap* inputFiles);
    virtual ~AliTestShuttle();

    virtual Int_t Store(const char* detector, TObject* object, AliCDBMetaData* metaData);
    virtual const char* GetFile(Int_t system, const char* detector, const char* id, const char* source);
    virtual TList* GetFileSources(Int_t system, const char* detector, const char* id);
    virtual void Log(const char* detector, const char* message);

  protected:
    TMap* fInputFiles;   // files for GetFile, GetFileSources

  private:
    ClassDef(AliTestShuttle, 0);
};

#endif
