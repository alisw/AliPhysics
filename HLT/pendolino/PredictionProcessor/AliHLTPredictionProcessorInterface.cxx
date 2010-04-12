// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Sebastian Bablok <Sebastian.Bablok@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTPredictionProcessorInterface.cxx
    @author Sebastian Bablok
    @date   
    @brief  
*/

#include "AliHLTPredictionProcessorInterface.h"
#include "AliHLTPendolino.h"
#include "TObjArray.h"
#include "AliDCSValue.h"


ClassImp(AliHLTPredictionProcessorInterface)

AliHLTPredictionProcessorInterface::AliHLTPredictionProcessorInterface(
			const char* detector, AliHLTPendolino* pendolino) 
  : AliPreprocessor(detector, reinterpret_cast<AliShuttleInterface*>(pendolino))
  , fpPend(pendolino)
  , fPredict(kFALSE)
  , fStartTime(0)
  , fEndTime(0)
{
}


AliHLTPredictionProcessorInterface::~AliHLTPredictionProcessorInterface() {

}

Int_t AliHLTPredictionProcessorInterface::GetRunNumber() {
	return fpPend->GetRunNumber();
}

void AliHLTPredictionProcessorInterface::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  // initializes AliHLTPredictionProcessorHLT
  if (GetRunNumber()!=run) {
    Log(Form("run number argument %d differs from pendolino run number %d", run, GetRunNumber()));
  }
  fStartTime = startTime;
  fEndTime = endTime;
}

Bool_t AliHLTPredictionProcessorInterface::includeAliCDBEntryInList(
            const TString& entryPath) {

    return fpPend->IncludeAliCDBEntryInList(entryPath);
}

Bool_t AliHLTPredictionProcessorInterface::GetSensorValue(TMap* dcsAliasMap,
							  const char* stringId, Float_t *value) const
{
  // extracts the sensor value
  // return last value read from sensor specified by stringId
  
  TObjArray* valueSet;
  TPair* pair = (TPair*)dcsAliasMap->FindObject(stringId);
  if (pair) {
    valueSet = (TObjArray*)pair->Value();
    Int_t nentriesDCS = valueSet->GetEntriesFast() - 1;
    if(nentriesDCS>=0){
      AliDCSValue *val = (AliDCSValue *)valueSet->At(nentriesDCS);
      *value=val->GetFloat();
      return kTRUE;
    }
  }
  return kFALSE;
}
