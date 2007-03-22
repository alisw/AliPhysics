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

#include <TTimeStamp.h>

//
// This class is the SHUTTLE preprocessor for the TPC detector.
// It contains several components, this far the part containing 
// temperatures is implemented
//

ClassImp(AliTPCPreprocessor)

//______________________________________________________________________________________________
AliTPCPreprocessor::AliTPCPreprocessor(const char* detector, AliShuttleInterface* shuttle) :
  AliPreprocessor(detector, shuttle),
  fTemp(0)
{
  // constructor
}

//______________________________________________________________________________________________
AliTPCPreprocessor::~AliTPCPreprocessor()
{
  // destructor
  
  delete fTemp;
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

        fTemp = new AliTPCSensorTempArray(fStartTime, fEndTime);
}

//______________________________________________________________________________________________
UInt_t AliTPCPreprocessor::Process(TMap* dcsAliasMap)
{
  // Fills data into TPC calibrations objects

  if (!dcsAliasMap) return 0;

  // Amanda servers provide information directly through dcsAliasMap

  // Temperature sensors are processed by AliTPCCalTemp
  
  UInt_t tempResult = MapTemperature(dcsAliasMap);
  UInt_t result=tempResult;
  
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

  TMap *map = fTemp->ExtractDCS(dcsAliasMap);
  if (map) {
    fTemp->MakeSplineFit(map);
    AliInfo(Form("Temperature values extracted, fits performed.\n"));
  } else {
    AliError(Form("No temperature map extracted.\n"));
    Log("AliTPCPreprocsessor: no temperature map extracted. \n");
  }
  delete map;
  // Now store the final CDB file
  
  AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("Haavard Helstrup");
	metaData.SetComment("Preprocessor AliTPC data base entries.");

	UInt_t result = Store("TPC/Calib/Temperature", "Data", fTemp, &metaData, 0, 0);
	delete fTemp;
	fTemp = 0;

   return result;
}
