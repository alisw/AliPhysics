// --- ROOT system
#include <TFile.h>
#include <TTimeStamp.h>

#include "AliZDCPreprocessor.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliZDCDataDCS.h"
#include "AliZDCCalibData.h"

//
// This class is an example for a simple preprocessor.
// It takes data from DCS and passes it to the class AliZDCDataDCS, which
// reformats its. This class is then written to the CDB.
//

ClassImp(AliZDCPreprocessor)

//______________________________________________________________________________________________
AliZDCPreprocessor::AliZDCPreprocessor(const char* detector, AliShuttleInterface* shuttle) :
  AliPreprocessor(detector, shuttle),
  fData(0)
{
  // constructor
}

//______________________________________________________________________________________________
AliZDCPreprocessor::~AliZDCPreprocessor()
{
  // destructor
}

//______________________________________________________________________________________________
void AliZDCPreprocessor::Initialize(Int_t run, UInt_t startTime,
	UInt_t endTime)
{
  // Creates AliZDCDataDCS object

  AliPreprocessor::Initialize(run, startTime, endTime);

	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run,
		TTimeStamp(startTime).AsString(),
		TTimeStamp(endTime).AsString()));

	fData = new AliZDCDataDCS(fRun, fStartTime, fEndTime);
}

//______________________________________________________________________________________________
UInt_t AliZDCPreprocessor::Process(TMap* dcsAliasMap)
{
  // Fills data into a AliZDCDataDCS object
  if(!dcsAliasMap) return 0;

  // The processing of the DCS input data is forwarded to AliZDCDataDCS
  Float_t DCSValues[26];
  fData->ProcessData(*dcsAliasMap, DCSValues);
  dcsAliasMap->Print("");
  //
  AliZDCCalibData *calibdata;
  calibdata->SetDCSCalibData(DCSValues);

  const char* PedfileName = GetFile(kDAQ, "PEDESTALS", "LDC0");
  if(PedfileName){
    AliInfo(Form("File %s connected to analyze pedestal events", PedfileName));
    //TFile *file = TFile::Open(PedfileName);
    AliCDBEntry *entry = AliCDBManager::Instance()->Get("ZDC/Calib/Data/");
    AliZDCCalibData *calibpeddata = (AliZDCCalibData*) entry->GetObject();
    calibpeddata->Print("");
    calibdata->SetMeanPed(calibpeddata->GetMeanPed());
    calibdata->SetMeanPedWidth(calibpeddata->GetMeanPedWidth());
    calibdata->SetOOTPed(calibpeddata->GetOOTPed());
    calibdata->SetOOTPedWidth(calibpeddata->GetOOTPedWidth());
    calibdata->SetPedCorrCoeff(calibpeddata->GetPedCorrCoeff());
  }
  else AliInfo(Form("File %s not found", PedfileName));

  const char* EMDfileName = GetFile(kDAQ, "MUTUALEMD", "LDC0");
  if(EMDfileName){
    AliInfo(Form("File %s connected to analyze EM dissociation events", EMDfileName));
    //TFile *file = TFile::Open(EMDfileName);
    AliCDBEntry *entry = AliCDBManager::Instance()->Get("ZDC/Calib/Data/");
    AliZDCCalibData *calibEMDdata = (AliZDCCalibData*) entry->GetObject();
    calibdata->SetEnCalib(calibEMDdata->GetEnCalib());
  }
  else AliInfo(Form("File %s not found", EMDfileName));
   
  calibdata->Print("");
  
  //Now we have to store the final CDB file
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Chiara");
  metaData.SetComment("This preprocessor fills an AliZDCDataDCS object.");

  UInt_t result = Store(fData, &metaData);
  delete fData;
  fData = 0;

  return result;
}

