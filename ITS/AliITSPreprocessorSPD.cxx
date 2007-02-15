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
UInt_t AliITSPreprocessorSPD::Process(TMap* dcsAliasMap)
{
  // Do the actual preprocessing

  UInt_t result = 1;


//  // *** REFERENCE DATA ***

//  // Store the container files as reference data (one for each equipment)
//  for (UInt_t eq=0; eq<20; eq++) {
//    Char_t id[20];
//    sprintf(id,"SPDref%d",eq);
//    TList* list = GetFileSources(kDAQ,id); // (the id is actually unique, so always 1 file)
//    if (list) {
//      TObjString* fileNameEntry = (TObjString*) list->First();
//      Char_t* fileName = (Char_t*) GetFile(kDAQ, id, fileNameEntry->GetString().Data());
//      AliITSOnlineSPDscan *scan = new AliITSOnlineSPDscan((Char_t*)fileName);
//      TObjArray* arr = scan->GetAsTObjArray();
//      Char_t refCAT[10];
//      sprintf(refCAT,"Ref%d",eq);
//      StoreReferenceData("SHUTTLE", refCAT, arr, &metaData, 0, 0);
//    }
//  }



  // *** NOISY DATA ***

  // Read old calibration
  TString runStr = GetRunParameter("run");
  Int_t run = atoi(runStr);
  AliCDBStorage* storage = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB());
  if(!storage) return 0;
  AliCDBEntry* cdbEntry = storage->Get("ITS/Calib/SPDNoisy", run);
  TObjArray* spdEntry;
  if(cdbEntry) {
    spdEntry = (TObjArray*)cdbEntry->GetObject();
    if(!spdEntry) return 0;
  }
  else {
    Log(Form("Old calibration not found in database. A fresh one will be created."));
    // Create fresh calibration
    spdEntry = new TObjArray(240);
    for(Int_t module=0;module<240;module++){
      AliITSCalibrationSPD* cal = new AliITSCalibrationSPD();
      spdEntry->Add(cal);
    }
    spdEntry->SetOwner(kTRUE);
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


  // Update the database entries
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
    Log(Form("Noisy lists for %d modules updated.",nrUpdatedMods));

    // Store the cdb entry
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0);
    metaData.SetResponsible("Henrik Tydesjo");
    metaData.SetComment("Preprocessor test for SPD.");  
    result = Store("Calib", "SPDNoisy", spdEntry, &metaData, 0, kTRUE);
    delete spdEntry;
    Log("Database updated.");
  }

  return result;



}

