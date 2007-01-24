///////////////////////////////////////////////
//  Author: Henrik Tydesjo                   //
//  Preprocessor Class for the SPD           //
//                                           //
///////////////////////////////////////////////

#include "AliITSPreprocessorSPD.h"
#include "AliITSCalibrationSPD.h"
#include "AliITSOnlineCalibrationSPD.h"
#include "AliITSOnlineCalibrationSPDhandler.h"
#include "AliCDBMetaData.h"
#include "AliLog.h"
#include <TTimeStamp.h>
#include <TObjString.h>

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

  // DEAD/NOISY DATA
  // Create the cdb entry, at first filled with empty cal objects (one for each module)
  TObjArray *spdEntry = new TObjArray(240);
  for(Int_t module=0;module<240;module++){
    AliITSCalibrationSPD* cal = new AliITSCalibrationSPD();
    spdEntry->Add(cal);
  }
  spdEntry->SetOwner(kTRUE);
  // Get all the files corresponding to the different modules and fill the dead/noisy lists
  AliITSOnlineCalibrationSPDhandler* handler = new AliITSOnlineCalibrationSPDhandler();
  for (Int_t module=0; module<240; module++) {
    Char_t id[20];
    sprintf(id,"SPD_%d",module);
    TList* list = GetFileSources(kDAQ,id); // (the id is actually unique, so always 1 file)
    if (list) {
      TObjString* fileNameEntry = (TObjString*) list->First();
      Char_t* fileName = (Char_t*) GetFile(kDAQ, id, fileNameEntry->GetString().Data());
      handler->SetModuleNr(module);
      handler->ReadFromFile(fileName);
      ((AliITSCalibrationSPD*) spdEntry->At(module)) -> SetNrDead( handler->GetNrDead() );
      ((AliITSCalibrationSPD*) spdEntry->At(module)) -> SetDeadList( handler->GetDeadArray() );
      ((AliITSCalibrationSPD*) spdEntry->At(module)) -> SetNrNoisy( handler->GetNrNoisy() );
      ((AliITSCalibrationSPD*) spdEntry->At(module)) -> SetNoisyList( handler->GetNoisyArray() );
    }
  }
  delete handler;
  // Store the cdb entry
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Henrik Tydesjo");
  metaData.SetComment("Preprocessor test for SPD.");  
  UInt_t result = Store("SHUTTLE", "Calib", spdEntry, &metaData, 0, 0);
  delete spdEntry;
  
//  // REFERENCE DATA
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

  return result;
}

