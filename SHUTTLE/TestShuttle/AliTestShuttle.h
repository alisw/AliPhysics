#ifndef ALI_TEST_SHUTTLE_H
#define ALI_TEST_SHUTTLE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// test implementation of the AliShuttleInterface, to be used for local tests of preprocessors
//

#include <AliShuttleInterface.h>
#include <TString.h>

class TMap;
class TList;
class AliCDBMetaData;
class AliCDBPath;

class AliTestShuttle : public AliShuttleInterface
{
  public:
    AliTestShuttle(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual ~AliTestShuttle();

    void AddInputFile(Int_t system, const char* detector, const char* id, const char* source, const char* fileName);
    void SetDCSInput(TMap* dcsAliasMap) { fDcsAliasMap = dcsAliasMap; }
    void AddInputRunParameter(const char* key, const char* value);

    void Process();

    // AliShuttleInterface functions
    virtual UInt_t Store(const AliCDBPath& path, TObject* object, AliCDBMetaData* metaData,
        				Int_t validityStart = 0, Bool_t validityInfinite = kFALSE);
    virtual UInt_t StoreReferenceData(const AliCDBPath& path, TObject* object, AliCDBMetaData* metaData);
    virtual const char* GetFile(Int_t system, const char* detector, const char* id, const char* source);
    virtual TList* GetFileSources(Int_t system, const char* detector, const char* id);
    virtual const char* GetRunParameter(const char* key);
    virtual void Log(const char* detector, const char* message);

    static TString GetMainCDB () {return fgkMainCDB;}
    static void SetMainCDB (TString mainCDB) {fgkMainCDB = mainCDB;}
    static TString GetMainRefStorage() {return fgkMainRefStorage;}
    static void SetMainRefStorage (TString mainRefStorage) {fgkMainRefStorage = mainRefStorage;}

    static void SetShuttleTempDir (const char* tmpDir);
    static const char* GetShuttleTempDir() {return fgkShuttleTempDir.Data();}
    static void SetShuttleLogDir (const char* logDir);
    static const char* GetShuttleLogDir() {return fgkShuttleLogDir.Data();}

    virtual void RegisterPreprocessor(AliPreprocessor* preprocessor);

  protected:
    static TString 	    fgkMainCDB;		// URI of the main (Grid) CDB storage
    static TString 	    fgkMainRefStorage;	// URI of the main (Grid) REFERENCE storage
    static TString 	    fgkShuttleTempDir;	// base path of SHUTTLE temp folder
    static TString 	    fgkShuttleLogDir;	// path of SHUTTLE log folder

    Int_t fRun;         // run that is simulated with the AliTestShuttle
    UInt_t fStartTime;  // starttime that is simulated with the AliTestShuttle
    UInt_t fEndTime;    // endtime that is simulated with the AliTestShuttle

    TMap* fInputFiles;      // files for GetFile, GetFileSources
    TMap* fRunParameters;    // run parameters
    TObjArray* fPreprocessors; // list of preprocessors that are to be tested
    TMap* fDcsAliasMap; // DCS data for testing

  private:
    ClassDef(AliTestShuttle, 0);
};

#endif
