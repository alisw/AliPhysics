///////////////////////////////////////////////
//  Author: Henrik Tydesjo                   //
//  Preprocessor Class for the SPD           //
//                                           //
///////////////////////////////////////////////

#include "AliITSPreprocessorSPD.h"
#include "AliITSCalibrationSPD.h"
#include "AliITSOnlineCalibrationSPDhandler.h"
#include "AliCDBEntry.h"
#include "AliCDBMetaData.h"
#include "AliShuttleInterface.h"
#include "AliLog.h"
#include <TTimeStamp.h>
#include <TObjString.h>
#include <TSystem.h>
#include <TList.h>

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


  // *** GET RUN TYPE ***

  TString runType = GetRunType();

  UInt_t nrEqForScan  = 0;
  UInt_t nrEqForPhysN = 0;
  UInt_t nrEqForPhysD = 0;

  // ******************************************************************************************** //
  // *** REFERENCE DATA *** //

  // Standalone runs:
  if (runType == "DAQ_MIN_TH_SCAN" ||
      runType == "DAQ_MEAN_TH_SCAN" ||
      runType == "DAQ_UNIFORMITY_SCAN" ||
      runType == "DAQ_NOISY_PIX_SCAN" ||
      runType == "DAQ_PIX_DELAY_SCAN" ||
      runType == "DAQ_FO_UNIF_SCAN") {
    // Store the scan container files as reference data (0 or 1 file for each equipment)
    for (UInt_t eq=0; eq<20; eq++) {
      TString id = Form("SPD_ref_scan_%d",eq);
      TList* list = GetFileSources(kDAQ,id.Data()); // (the id should be unique, so always 1 file)
      if (list) {
	TObjString* fileNameEntry = (TObjString*) list->First();
	if (fileNameEntry!=NULL) {
	  nrEqForScan++;
	  TString fileName = GetFile(kDAQ, id, fileNameEntry->GetString().Data());
	  TString refCAT = Form("SPD_ref_scan_eq_%d",eq);
	  if (!StoreReferenceFile(fileName.Data(),refCAT.Data())) {
	    Log(Form("Failed to store reference file %s.",fileName.Data()));
	    return 1;
	  }
	}
      }
    }
  }

  // Physics runs (online monitoring):
  if (runType == "PHYSICS") {
    // Store the phys "per run" container files as reference data (0 or 1 file for each equipment)
    for (UInt_t eq=0; eq<20; eq++) {
      TString id = Form("SPD_ref_phys_%d",eq);
      TList* list = GetFileSources(kDAQ,id.Data()); // (the id should be unique, so always 1 file)
      if (list) {
	TObjString* fileNameEntry = (TObjString*) list->First();
	if (fileNameEntry!=NULL) {
	  nrEqForPhysN++;
	  TString fileName = GetFile(kDAQ, id, fileNameEntry->GetString().Data());
	  TString refCAT = Form("SPD_ref_phys_eq_%d",eq);
	  if (!StoreReferenceFile(fileName.Data(),refCAT.Data())) {
	    Log(Form("Failed to store reference file %s.",fileName.Data()));
	    return 1;
	  }
	}
      }
    }
    // Store the phys "dead" container files as reference data (0 or 1 file for each equipment)
    for (UInt_t eq=0; eq<20; eq++) {
      TString id = Form("SPD_ref_phys_dead_%d",eq);
      TList* list = GetFileSources(kDAQ,id.Data()); // (the id should be unique, so always 1 file)
      if (list) {
	TObjString* fileNameEntry = (TObjString*) list->First();
	if (fileNameEntry!=NULL) {
	  nrEqForPhysD++;
	  TString fileName = GetFile(kDAQ, id, fileNameEntry->GetString().Data());
	  TString refCAT = Form("SPD_ref_phys_dead_eq_%d",eq);
	  if (!StoreReferenceFile(fileName.Data(),refCAT.Data())) {
	    Log(Form("Failed to store reference file %s.",fileName.Data()));
	    return 1;
	  }
	}
      }
    }
  }

  // ******************************************************************************************** //


  // *** NOISY AND DEAD DATA *** //

  // Standalone runs:
  if (runType == "DAQ_NOISY_PIX_SCAN") {
    // Retrieve and unpack tared calibration files from FXS
    TList* list = GetFileSources(kDAQ,"SPD_scan_noisy");
    if (list) {
      UInt_t index = 0;
      while (list->At(index)!=NULL) {
	TObjString* fileNameEntry = (TObjString*) list->At(index);
	TString fileName = GetFile(kDAQ, "SPD_scan_noisy", fileNameEntry->GetString().Data());
	TString command = Form("tar -xf %s",fileName.Data());
	gSystem->Exec(command.Data());
	index++;
      }
    }
    // Create new database entries
    TObjArray* spdEntryNoisy = new TObjArray(240);
    spdEntryNoisy->SetOwner(kTRUE);
    for(UInt_t module=0; module<240; module++){
      AliITSCalibrationSPD* calObj = new AliITSCalibrationSPD();
      spdEntryNoisy->Add(calObj);
    }
    // Add noisy from the copied FXS files
    AliITSOnlineCalibrationSPDhandler* handler = new AliITSOnlineCalibrationSPDhandler();
    TString fileLoc = ".";
    handler->SetFileLocation(fileLoc.Data());
    handler->ReadNoisyFromFiles();
    for (Int_t module=0; module<240; module++) {
      ((AliITSCalibrationSPD*) spdEntryNoisy->At(module)) -> SetNrBad( handler->GetNrNoisy(module) );
      ((AliITSCalibrationSPD*) spdEntryNoisy->At(module)) -> SetBadList( handler->GetNoisyArray(module) );
    }
    delete handler;
    // Store the new calibration objects in OCDB
    Log("Noisy lists (scan) will be stored...");
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0);
    metaData.SetResponsible("Henrik Tydesjo");
    metaData.SetComment("Created by SPD PreProcessor");
    // validity for this run until infinity
    if (!Store("Calib", "SPDNoisy", spdEntryNoisy, &metaData, 0, kTRUE)) {
      Log("Failed to store calibration data.");
      return 1;
    }
    Log("Database updated.");
    delete spdEntryNoisy;
  }

  // Physics runs (online monitoring):
  else if (runType == "PHYSICS") {

    // Noisy pixels:
    // Read noisy from previous calibration
    AliCDBEntry* cdbEntry = GetFromOCDB("Calib", "SPDNoisy");
    TObjArray* spdEntryNoisy;
    if(cdbEntry) {
      spdEntryNoisy = (TObjArray*)cdbEntry->GetObject();
      if(!spdEntryNoisy) return 1;
    }
    else {
      Log("Old calibration not found in database. This is required for further processing.");
      return 1;
    }
    AliITSOnlineCalibrationSPDhandler* handOld = new AliITSOnlineCalibrationSPDhandler();
    handOld->ReadNoisyFromCalibObj(spdEntryNoisy);
    // Retrieve and unpack tared calibration files from FXS
    TList* listN = GetFileSources(kDAQ,"SPD_phys_noisy");
    if (listN) {
      UInt_t index = 0;
      while (listN->At(index)!=NULL) {
	TObjString* fileNameEntry = (TObjString*) listN->At(index);
	TString fileName = GetFile(kDAQ, "SPD_phys_noisy", fileNameEntry->GetString().Data());
	TString command = Form("tar -xf %s",fileName.Data());
	gSystem->Exec(command.Data());
	index++;
      }
    }
    AliITSOnlineCalibrationSPDhandler* handNew = new AliITSOnlineCalibrationSPDhandler();
    handNew->SetFileLocation(".");
    handNew->ReadNoisyFromFiles();
    // add the new list to the old one
    UInt_t nrNewNoisy = handOld->AddNoisyFrom(handNew);
    // If new noisy pixels were found: Update calibration objects
    if (nrNewNoisy>0) {
      for (Int_t module=0; module<240; module++) {
	((AliITSCalibrationSPD*) spdEntryNoisy->At(module)) -> SetNrBad( handOld->GetNrNoisy(module) );
	((AliITSCalibrationSPD*) spdEntryNoisy->At(module)) -> SetBadList( handOld->GetNoisyArray(module) );
      }
      // Store the new calibration objects in OCDB
      Log("Noisy lists (phys) will be stored...");
      AliCDBMetaData metaData;
      metaData.SetBeamPeriod(0);
      metaData.SetResponsible("Henrik Tydesjo");
      metaData.SetComment("Created by SPD PreProcessor");  
      // validity for this run only
      if (!Store("Calib", "SPDNoisy", spdEntryNoisy, &metaData, 0, kFALSE)) {
	Log("Failed to store calibration data.");
	return 1;
      }
      Log("Database updated.");
    }
    delete handNew;

    // Dead pixels:
    // Retrieve and unpack tared calibration files from FXS
    TList* listD = GetFileSources(kDAQ,"SPD_phys_dead");
    if (listD) {
      UInt_t index = 0;
      while (listD->At(index)!=NULL) {
	TObjString* fileNameEntry = (TObjString*) listD->At(index);
	TString fileName = GetFile(kDAQ, "SPD_phys_dead", fileNameEntry->GetString().Data());
	TString command = Form("tar -xf %s",fileName.Data());
	gSystem->Exec(command.Data());
	index++;
      }
    }
    // Create new database entries
    TObjArray* spdEntryDead = new TObjArray(240);
    spdEntryDead->SetOwner(kTRUE);
    for(UInt_t module=0; module<240; module++){
      AliITSCalibrationSPD* calObj = new AliITSCalibrationSPD();
      spdEntryDead->Add(calObj);
    }
    // Add dead from the copied FXS files
    handOld->SetFileLocation(".");
    handOld->ReadDeadFromFiles();
    for (Int_t module=0; module<240; module++) {
      ((AliITSCalibrationSPD*) spdEntryDead->At(module)) -> SetNrBad( handOld->GetNrDead(module) );
      ((AliITSCalibrationSPD*) spdEntryDead->At(module)) -> SetBadList( handOld->GetDeadArray(module) );
    }
    delete handOld;
    // Store the new calibration objects in OCDB
    Log("Dead lists (phys) will be stored...");
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0);
    metaData.SetResponsible("Henrik Tydesjo");
    metaData.SetComment("Created by SPD PreProcessor");
    // validity for this run only
    if (!Store("Calib", "SPDDead", spdEntryDead, &metaData, 0, kFALSE)) {
      Log("Failed to store calibration data.");
      return 1;
    }
    Log("Database updated.");
    delete spdEntryDead;

  }


  return 0; // 0 means success

}

