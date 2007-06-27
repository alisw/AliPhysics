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

#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliTPCSensorTempArray.h"
#include "AliTPCDBPressure.h"

#include <TTimeStamp.h>

const Int_t kValCutTemp = 100;         // discard temperatures > 100 degrees
const Int_t kDiffCutTemp = 5;	   // discard temperature differences > 5 degrees

//
// This class is the SHUTTLE preprocessor for the TPC detector.
// It contains several components, this far the part containing
// temperatures is implemented
//

ClassImp(AliTPCPreprocessor)

//______________________________________________________________________________________________
AliTPCPreprocessor::AliTPCPreprocessor(AliShuttleInterface* shuttle) :
  AliPreprocessor("TPC",shuttle),
  fTemp(0), fPressure(0), fConfigOK(kTRUE)
{
  // constructor
}
//______________________________________________________________________________________________
// AliTPCPreprocessor::AliTPCPreprocessor(const AliTPCPreprocessor& org) :
//   AliPreprocessor(org),
//   fTemp(0), fPressure(0), fConfigOK(kTRUE)
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

  // Temperature sensors

        TTree *confTree = 0;
	AliCDBEntry* entry = GetFromOCDB("Config", "Temperature");
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

   if (!dcsAliasMap) return 9;
   if (!fConfigOK) return 9;

  // Amanda servers provide information directly through dcsAliasMap

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

	result = Store("Calib", "Temperature", fTemp, &metaData, 0, 0);
        if ( result == 1 ) {                  
          result = 0;
	} else {
	  result = 1;
	}                      // revert to new return code conventions
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

	result = Store("Calib", "Pressure", fPressure, &metaData, 0, 0);
        if ( result == 1 ) {                  
          result = 0;
	} else {
	  result = 1;
	}                      // revert to new return code conventions
   }

   return result;
}
