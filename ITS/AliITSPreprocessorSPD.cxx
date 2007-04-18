///////////////////////////////////////////////
//  Author: Henrik Tydesjo                   //
//  Preprocessor Class for the SPD           //
//                                           //
///////////////////////////////////////////////

#include "AliITSPreprocessorSPD.h"
#include "AliITSCalibrationSPD.h"
#include "AliITSOnlineCalibrationSPD.h"
#include "AliITSOnlineCalibrationSPDhandler.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliCDBMetaData.h"
#include "AliShuttleInterface.h"
#include "AliLog.h"
#include <TTimeStamp.h>
#include <TObjString.h>
#include <TSystem.h>

ClassImp(AliITSPreprocessorSPD)

//______________________________________________________________________________________________
AliITSPreprocessorSPD::AliITSPreprocessorSPD(AliShuttleInterface* shuttle) :
  AliPreprocessor("SPD", shuttle)
{
  // constructor
}

//______________________________________________________________________________________________
AliITSPreprocessorSPD::~AliITSPreprocessorSPD()
{
  // destructor
}

//______________________________________________________________________________________________
void AliITSPreprocessorSPD::Initialize(Int_t run, UInt_t startTime,
	UInt_t endTime)
{
  // initialize
  AliPreprocessor::Initialize(run, startTime, endTime);

  AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run,
	       TTimeStamp(startTime).AsString(),
	       TTimeStamp(endTime).AsString()));
}

//______________________________________________________________________________________________
UInt_t AliITSPreprocessorSPD::Process(TMap* /*dcsAliasMap*/)
{
  // Do the actual preprocessing


  // *** CHECK RUN TYPE ***

  TString runType = GetRunType();
  if(runType != "SPD_STANDALONE_CALIBRATION") {
    Log("Nothing to do");
    return 0;
  }


  // *** REFERENCE DATA ***

  // Store the container files as reference data (one file for each equipment)
  for (UInt_t eq=0; eq<20; eq++) {
    Char_t id[20];
    sprintf(id,"SPD_reference_%d",eq);
    TList* list = GetFileSources(kDAQ,id); // (the id should be unique, so always 1 file)
    if (list) {
      TObjString* fileNameEntry = (TObjString*) list->First();
      if (fileNameEntry!=NULL) {
	TString fileName = GetFile(kDAQ, id, fileNameEntry->GetString().Data());
	Char_t refCAT[10];
	sprintf(refCAT,"SPDref_eq_%d.root",eq);
	if (!StoreReferenceFile(fileName.Data(),refCAT)) {
	  return 1;
	}
      }
    }
  }


  // *** NOISY DATA ***

  // Read old calibration
  AliCDBEntry* cdbEntry = GetFromOCDB("Calib", "CalibSPD");
  TObjArray* spdEntry;
  if(cdbEntry) {
    spdEntry = (TObjArray*)cdbEntry->GetObject();
    if(!spdEntry) return 1;
  }
  else {
    Log("Old calibration not found in database. This is required for further processing.");
    return 1;
  }

  TString pwd = gSystem->pwd();  // remember the working directory, to cd back later
  TString tempDir = Form("%s",AliShuttleInterface::GetShuttleTempDir());

  // Retrieve and unpack tared calibration files from FXS
  TList* list = GetFileSources(kDAQ,"SPD_noisy");
  if (list) {
    UInt_t index = 0;
    while (list->At(index)!=NULL) {
      TObjString* fileNameEntry = (TObjString*) list->At(index);
      TString fileName = GetFile(kDAQ, "SPD_noisy", fileNameEntry->GetString().Data());
      gSystem->cd(tempDir.Data());
      Char_t command[100];
      sprintf(command,"tar -xf %s",fileName.Data());
      gSystem->Exec(command);
      index++;
    }
  }

  gSystem->cd(pwd.Data());

  // Update the database entries if needed
  UInt_t nrUpdatedMods = 0;
  AliITSOnlineCalibrationSPDhandler* handler = new AliITSOnlineCalibrationSPDhandler();
  Char_t fileLoc[100];
  sprintf(fileLoc,"%s",AliShuttleInterface::GetShuttleTempDir());
  handler->SetFileLocation(fileLoc);
  for (Int_t module=0; module<240; module++) {
    handler->SetModuleNr(module);
    if (handler->ReadFromFile()) {
      ((AliITSCalibrationSPD*) spdEntry->At(module)) -> SetNrNoisy( handler->GetNrNoisy() );
      ((AliITSCalibrationSPD*) spdEntry->At(module)) -> SetNoisyList( handler->GetNoisyArray() );
      nrUpdatedMods++;
    }
  }
  delete handler;
  if (nrUpdatedMods>0) {
    Log(Form("Noisy lists for %d modules will be updated and stored...",nrUpdatedMods));
    // Store the cdb entry
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0);
    metaData.SetResponsible("Henrik Tydesjo");
    metaData.SetComment("Preprocessor test for SPD.");  
    if (!Store("Calib", "CalibSPD", spdEntry, &metaData, 0, kTRUE)) {
      return 1;
    }
    //    delete spdEntry;
    Log("Database updated.");
  }

  return 0; // 0 means success



}

