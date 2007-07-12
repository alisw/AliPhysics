/**************************************************************************
 * Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#include "AliTPCPreprocessor.h"
#include "AliShuttleInterface.h"

#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliTPCSensorTempArray.h"
#include "AliTPCROC.h"
#include "AliTPCCalROC.h"
#include "AliTPCCalPad.h"
#include "AliTPCCalibPedestal.h"
#include "TFile.h"
#include "TTree.h"
#include "TEnv.h"

#include <TTimeStamp.h>

const Int_t kValCutTemp = 100;               // discard temperatures > 100 degrees
const Int_t kDiffCutTemp = 5;	             // discard temperature differences > 5 degrees
const TString kPedestalRunType = "PEDESTAL_RUN";  // pedestal run identifier

//
// This class is the SHUTTLE preprocessor for the TPC detector.
// It contains several components, this far the part containing
// temperatures is implemented
//

ClassImp(AliTPCPreprocessor)

//______________________________________________________________________________________________
AliTPCPreprocessor::AliTPCPreprocessor(AliShuttleInterface* shuttle) :
  AliPreprocessor("TPC",shuttle),
  fConfEnv(0), fTemp(0), fPressure(0), fConfigOK(kTRUE), fROC(0)
{
  // constructor
  fROC = AliTPCROC::Instance();
}
//______________________________________________________________________________________________
// AliTPCPreprocessor::AliTPCPreprocessor(const AliTPCPreprocessor& org) :
//   AliPreprocessor(org),
//   fConfEnv(0), fTemp(0), fPressure(0), fConfigOK(kTRUE)
// {
//   // copy constructor not implemented
//   //   -- missing underlying copy constructor in AliPreprocessor
//
//   Fatal("AliTPCPreprocessor", "copy constructor not implemented");
//
// //  fTemp = new AliTPCSensorTempArray(*(org.fTemp));
// }

//______________________________________________________________________________________________
AliTPCPreprocessor::~AliTPCPreprocessor()
{
  // destructor

  delete fTemp;
  delete fPressure;
}
//______________________________________________________________________________________________
AliTPCPreprocessor& AliTPCPreprocessor::operator = (const AliTPCPreprocessor& )
{
  Fatal("operator =", "assignment operator not implemented");
  return *this;
}


//______________________________________________________________________________________________
void AliTPCPreprocessor::Initialize(Int_t run, UInt_t startTime,
	UInt_t endTime)
{
  // Creates AliTestDataDCS object

  AliPreprocessor::Initialize(run, startTime, endTime);

	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run,
		TTimeStamp(startTime).AsString(),
		TTimeStamp(endTime).AsString()));

  // Preprocessor configuration

	AliCDBEntry* entry = GetFromOCDB("Config", "Temperature");
        if (entry) fConfEnv = (TEnv*) entry->GetObject();
        if ( fConfEnv==0 ) {
           AliWarning(Form("Preprocessor Config OCDB entry missing.\n"));
           Log("AliTPCPreprocsessor: Preprocessor Config OCDB entry missing.\n");
        }

  // Temperature sensors

        TTree *confTree = 0;
	entry = GetFromOCDB("Config", "Temperature");
        if (entry) confTree = (TTree*) entry->GetObject();
        if ( confTree==0 ) {
           AliError(Form("Temperature Config OCDB entry missing.\n"));
           Log("AliTPCPreprocsessor: Temperature Config OCDB entry missing.\n");
	   fConfigOK = kFALSE;
	   return;
        }
        fTemp = new AliTPCSensorTempArray(fStartTime, fEndTime, confTree);
	fTemp->SetValCut(kValCutTemp);
	fTemp->SetDiffCut(kDiffCutTemp);

  // Pressure sensors

        confTree=0;
	entry=0;
	entry = GetFromOCDB("Config", "Pressure");
        if (entry) confTree = (TTree*) entry->GetObject();
        if ( confTree==0 ) {
           AliError(Form("Pressure Config OCDB entry missing.\n"));
           Log("AliTPCPreprocsessor: Pressure Config OCDB entry missing.\n");
	   fConfigOK = kFALSE;
	   return;
        }
	fPressure = new AliDCSSensorArray(fStartTime, fEndTime, confTree);

}

//______________________________________________________________________________________________
UInt_t AliTPCPreprocessor::Process(TMap* dcsAliasMap)
{
  // Fills data into TPC calibrations objects

  // Amanda servers provide information directly through dcsAliasMap

  if (!dcsAliasMap) return 9;
  if (dcsAliasMap->GetEntries() == 0 ) return 9;
  if (!fConfigOK) return 9;

  TString runType = GetRunType();

  // Temperature sensors are processed by AliTPCCalTemp


  UInt_t tempResult = MapTemperature(dcsAliasMap);
  UInt_t result=tempResult;

  // Pressure sensors

  UInt_t pressureResult = MapPressure(dcsAliasMap);
  result += pressureResult;

  // Other calibration information will be retrieved through FXS files
  //  examples:
  //    TList* fileSourcesDAQ = GetFile(AliShuttleInterface::kDAQ, "pedestals");
  //    const char* fileNamePed = GetFile(AliShuttleInterface::kDAQ, "pedestals", "LDC1");
  //
  //    TList* fileSourcesHLT = GetFile(AliShuttleInterface::kHLT, "calib");
  //    const char* fileNameHLT = GetFile(AliShuttleInterface::kHLT, "calib", "LDC1");


  if(runType == kPedestalRunType) {
    Int_t pedestalSource = AliShuttleInterface::kDAQ;
    TString source = fConfEnv->GetValue("Pedestal","DAQ");
    source.ToUpper();
    if ( source == "HLT" ) pedestalSource = AliShuttleInterface::kHLT;
    UInt_t pedestalResult = ExtractPedestals(pedestalSource);
    result += pedestalResult;

  }


  return result;
}
//______________________________________________________________________________________________
UInt_t AliTPCPreprocessor::MapTemperature(TMap* dcsAliasMap)
{

   // extract DCS temperature maps. Perform fits to save space

  UInt_t result=0;
  TMap *map = fTemp->ExtractDCS(dcsAliasMap);
  if (map) {
    fTemp->MakeSplineFit(map);
    AliInfo(Form("Temperature values extracted, fits performed.\n"));
  } else {
    AliError(Form("No temperature map extracted.\n"));
    Log("AliTPCPreprocsessor: no temperature map extracted. \n");
    result=9;
  }
  delete map;
  // Now store the final CDB file

  if ( result == 0 ) {
        AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("Haavard Helstrup");
	metaData.SetComment("Preprocessor AliTPC data base entries.");

	Bool_t storeOK = Store("Calib", "Temperature", fTemp, &metaData, 0, kFALSE);
        if ( !storeOK )  result=1;

   }

   return result;

}
//______________________________________________________________________________________________

UInt_t AliTPCPreprocessor::MapPressure(TMap* dcsAliasMap)
{

   // extract DCS temperature maps. Perform fits to save space

  UInt_t result=0;
  TMap *map = fPressure->ExtractDCS(dcsAliasMap);
  if (map) {
    fPressure->MakeSplineFit(map);
    AliInfo(Form("Pressure values extracted, fits performed.\n"));
  } else {
    AliError(Form("No atmospheric pressure map extracted.\n"));
    Log("AliTPCPreprocsessor: no atmospheric pressure map extracted. \n");
    result=9;
  }
  delete map;
  // Now store the final CDB file

  if ( result == 0 ) {
        AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("Haavard Helstrup");
	metaData.SetComment("Preprocessor AliTPC data base entries.");

	Bool_t storeOK = Store("Calib", "Pressure", fPressure, &metaData, 0, 0);
        if ( !storeOK ) result=1;

   }

   return result;

}


//______________________________________________________________________________________________

UInt_t AliTPCPreprocessor::ExtractPedestals(Int_t sourceFXS)
{
 //
 //  Read pedestal file from file exchage server
 //  Keep original entry from OCDB in case no new pedestals are available
 //
 AliTPCCalPad *calPadPed=0;
 AliCDBEntry* entry = GetFromOCDB("Calib", "Pedestals");
 if (entry) calPadPed = (AliTPCCalPad*)entry->GetObject();
 if ( calPadPed==NULL ) {
     AliWarning(Form("No previous TPC pedestal entry available.\n"));
     Log("AliTPCPreprocsessor: No previous TPC pedestal entry available.\n");
     calPadPed = new AliTPCCalPad("PedestalsMean","PedestalsMean");
 }

 AliTPCCalPad *calPadRMS=0;
 entry = GetFromOCDB("Calib", "Noise");
 if (entry) calPadRMS = (AliTPCCalPad*)entry->GetObject();
 if ( calPadRMS==NULL ) {
     AliWarning(Form("No previous TPC noise entry available.\n"));
     Log("AliTPCPreprocsessor: No previous TPC noise entry available.\n");
     calPadRMS = new AliTPCCalPad("PedestalsRMS","PedestalsRMS");
 }


 UInt_t result=0;

 Int_t nSectors = fROC->GetNSectors();
 TList* list = GetFileSources(sourceFXS,"pedestals");
 if (list) {

//  loop through all files from LDCs

    UInt_t index = 0;
    while (list->At(index)!=NULL) {
     TObjString* fileNameEntry = (TObjString*) list->At(index);
     if (fileNameEntry!=NULL) {
        TString fileName = GetFile(sourceFXS, "pedestals",
	                                 fileNameEntry->GetString().Data());
        TFile *f = TFile::Open(fileName);
        AliTPCCalibPedestal *calPed;
	f->GetObject("AliTPCCalibPedestal",calPed);

        //  replace entries for the sectors available in the present file

        for (Int_t sector=0; sector<nSectors; sector++) {
           AliTPCCalROC *rocPed=calPed->GetCalRocPedestal(sector, kFALSE);
           if ( rocPed )  calPadPed->SetCalROC(rocPed,sector);
           AliTPCCalROC *rocRMS=calPed->GetCalRocRMS(sector, kFALSE);
           if ( rocRMS )  calPadRMS->SetCalROC(rocRMS,sector);
        }
      }
     ++index;
    }  // while(list)
//
//  Store updated pedestal entry to OCDB
//
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0);
    metaData.SetResponsible("Haavard Helstrup");
    metaData.SetComment("Preprocessor AliTPC data base entries.");

    Bool_t storeOK = Store("Calib", "Pedestals", calPadPed, &metaData, 0, kTRUE);
    if ( !storeOK ) ++result;
    storeOK = Store("Calib", "PadNoise", calPadRMS, &metaData, 0, kTRUE);
    if ( !storeOK ) ++result;

  }

  return result;
}

//______________________________________________________________________________________________


