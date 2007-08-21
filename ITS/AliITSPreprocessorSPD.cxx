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



  // *** REFERENCE DATA *** //

  // Standalone runs:
  if (runType == "DAQ_MIN_TH_SCAN" ||
      runType == "DAQ_MEAN_TH_SCAN" ||
      runType == "DAQ_UNIFORMITY_SCAN" ||
      runType == "DAQ_NOISY_PIX_SCAN" ||
      runType == "DAQ_PIX_DELAY_SCAN" ||
      runType == "DAQ_FO_UNIF_SCAN") {
    // Store the scan container files as reference data (one file for each equipment)
    for (UInt_t eq=0; eq<20; eq++) {
      TString id = Form("SPD_reference_%d",eq);
      TList* list = GetFileSources(kDAQ,id.Data()); // (the id should be unique, so always 1 file)
      if (list) {
	TObjString* fileNameEntry = (TObjString*) list->First();
	if (fileNameEntry!=NULL) {
	  TString fileName = GetFile(kDAQ, id, fileNameEntry->GetString().Data());
	  TString refCAT = Form("SPDref_eq_%d.root",eq);
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
    // *** code to be written *** //
  }



  // *** NOISY AND DEAD DATA *** //

  if (runType == "DAQ_NOISY_PIX_SCAN" || runType == "PHYSICS") {
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

    // Standalone runs:
    if (runType == "DAQ_NOISY_PIX_SCAN") {
      UInt_t nrUpdatedMods = 0;
      // Retrieve and unpack tared calibration files from FXS
      TList* list = GetFileSources(kDAQ,"SPD_noisy");
      if (list) {
	UInt_t index = 0;
	while (list->At(index)!=NULL) {
	  TObjString* fileNameEntry = (TObjString*) list->At(index);
	  TString fileName = GetFile(kDAQ, "SPD_noisy", fileNameEntry->GetString().Data());
	  TString command = Form("tar -xf %s",fileName.Data());
	  gSystem->Exec(command.Data());
	  index++;
	}
      }
      // Update the database entries for the modules that were scanned
      AliITSOnlineCalibrationSPDhandler* handler = new AliITSOnlineCalibrationSPDhandler();
      TString fileLoc = ".";
      handler->SetFileLocation(fileLoc.Data());
      for (Int_t module=0; module<240; module++) {
	if (handler->ReadFromFile(module)) {
	  ((AliITSCalibrationSPD*) spdEntry->At(module)) -> SetNrNoisy( handler->GetNrNoisy(module) );
	  ((AliITSCalibrationSPD*) spdEntry->At(module)) -> SetNoisyList( handler->GetNoisyArray(module) );
	  nrUpdatedMods++;
	}
      }
      delete handler;
      // Store the new calibration objects (if any modifications were made) in OCDB
      if (nrUpdatedMods>0) {
	Log(Form("Noisy lists for %d modules will be updated and stored...",nrUpdatedMods));
	AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("Henrik Tydesjo");
	metaData.SetComment("Preprocessor test for SPD.");  
	if (!Store("Calib", "CalibSPD", spdEntry, &metaData, 0, kTRUE)) {
	  Log("Failed to store calibration data.");
	  return 1;
	}
	Log("Database updated.");
      }
    }

    // Physics runs (online monitoring):
    if (runType == "PHYSICS") {
      // *** code to be written *** //
    }


  }



  return 0; // 0 means success

}

