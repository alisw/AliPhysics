// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

// @file   AliHLTModulePreprocessor.cxx
// @author Matthias Richter
// @date   2008-01-22
// @brief  Base class for HLT module preprocessors
// 

#include <cstdlib>
#include <cassert>
#include "AliHLTModulePreprocessor.h"
#include "AliHLTShuttleInterface.h"
#include "TObjString.h"
#include "TString.h"
#include "TMap.h"
#include "TObject.h"
#include "TObjArray.h"


ClassImp(AliHLTModulePreprocessor)

AliHLTModulePreprocessor::AliHLTModulePreprocessor() 
  :
  fpInterface(NULL),
  fActiveDetectors(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTModulePreprocessor::~AliHLTModulePreprocessor() 
{
  // see header file for function documentation
}

// TODO: map this constants to AliHLTDAQ class
#ifndef MFT_UPGRADE
const Int_t AliHLTModulePreprocessor::kNDetectors = 22;
#else
const Int_t AliHLTModulePreprocessor::kNDetectors = 21;
#endif

const char* AliHLTModulePreprocessor::fgkDetectorName[kNDetectors] = 
{
  "ITSSPD",
  "ITSSDD",
  "ITSSSD",
  "TPC",
  "TRD",
  "TOF",
  "HMPID",
  "PHOS",
  "CPV",
  "PMD",
  "MUONTRK",
  "MUONTRG",
  "FMD",
  "T0",
  "VZERO", // Name to be changed to V0 ?
  "ZDC",
  "ACORDE",
  "TRG",
  "EMCAL",
  "DAQ_TEST",
  "HLT"
  #ifndef MFT_UPGRADE
  , "MFT"
  #endif
};

void AliHLTModulePreprocessor::SetShuttleInterface(AliHLTShuttleInterface* pInterface)
{
  assert(fpInterface==NULL || fpInterface==pInterface || pInterface==NULL);
  fpInterface=pInterface;
}

Int_t AliHLTModulePreprocessor::GetRun() const
{
  // see header file for function documentation

  assert(fpInterface);
  if (!fpInterface) return 0;
  return fpInterface->PreprocessorGetRun();
}

UInt_t AliHLTModulePreprocessor::GetStartTime() const
{
  // see header file for function documentation

  assert(fpInterface);
  if (!fpInterface) return 0;
  return fpInterface->PreprocessorGetStartTime();
}

UInt_t AliHLTModulePreprocessor::GetEndTime() const
{
  // see header file for function documentation

  assert(fpInterface);
  if (!fpInterface) return 0;
  return fpInterface->PreprocessorGetEndTime();
}

Bool_t AliHLTModulePreprocessor::Store(const char* pathLevel2, const char* pathLevel3, TObject* object,
				 AliCDBMetaData* metaData, Int_t validityStart, Bool_t validityInfinite)
{
  // see header file for function documentation

  assert(fpInterface);
  if (!fpInterface) return 0;
  return fpInterface->PreprocessorStore(pathLevel2, pathLevel3, object, metaData, validityStart, validityInfinite);
}

Bool_t AliHLTModulePreprocessor::StoreReferenceData(const char* pathLevel2, const char* pathLevel3, TObject* object,
					      AliCDBMetaData* metaData)
{
  // see header file for function documentation

  assert(fpInterface);
  if (!fpInterface) return 0;
  return fpInterface->PreprocessorStoreReferenceData(pathLevel2, pathLevel3, object, metaData);
}

Bool_t AliHLTModulePreprocessor::StoreReferenceFile(const char* localFile, const char* gridFileName)
{
  // see header file for function documentation

  assert(fpInterface);
  if (!fpInterface) return 0;
  return fpInterface->PreprocessorStoreReferenceFile(localFile, gridFileName);
}

Bool_t AliHLTModulePreprocessor::StoreRunMetadataFile(const char* localFile, const char* gridFileName)
{
  // see header file for function documentation

  assert(fpInterface);
  if (!fpInterface) return 0;
  return fpInterface->PreprocessorStoreRunMetadataFile(localFile, gridFileName);
}
    
const char* AliHLTModulePreprocessor::GetFile(Int_t system, const char* id, const char* source)
{
  // see header file for function documentation

  assert(fpInterface);
  if (!fpInterface) return 0;
  return fpInterface->PreprocessorGetFile(system, id, source);
}

TList* AliHLTModulePreprocessor::GetFileSources(Int_t system, const char* id)
{
  // see header file for function documentation

  assert(fpInterface);
  if (!fpInterface) return 0;
  return fpInterface->PreprocessorGetFileSources(system, id);
}

TList* AliHLTModulePreprocessor::GetFileIDs(Int_t system, const char* source)
{
  // see header file for function documentation

  assert(fpInterface);
  if (!fpInterface) return 0;
  return fpInterface->PreprocessorGetFileIDs(system, source);
}

const char* AliHLTModulePreprocessor::GetRunParameter(const char* param)
{
  // see header file for function documentation

  assert(fpInterface);
  if (!fpInterface) return 0;
  return fpInterface->PreprocessorGetRunParameter(param);
}

AliCDBEntry* AliHLTModulePreprocessor::GetFromOCDB(const char* pathLevel2, const char* pathLevel3)
{
  // see header file for function documentation

  assert(fpInterface);
  if (!fpInterface) return 0;
  return fpInterface->PreprocessorGetFromOCDB(pathLevel2, pathLevel3);
}

const char* AliHLTModulePreprocessor::GetRunType()
{
  // see header file for function documentation

  assert(fpInterface);
  if (!fpInterface) return 0;
  return fpInterface->PreprocessorGetRunType();
}

void AliHLTModulePreprocessor::Log(const char* message)
{
  // see header file for function documentation

  assert(fpInterface);
  if (!fpInterface) return;
  fpInterface->PreprocessorLog(message);
}

Int_t AliHLTModulePreprocessor::DetectorBitMask(const char *detectorName)
{
  // Return the detector index
  // corresponding to a given
  // detector name
  TString detStr = detectorName;

  Int_t iDet;
  for(iDet = 0; iDet < kNDetectors; iDet++) {
    if (detStr.CompareTo(fgkDetectorName[iDet],TString::kIgnoreCase) == 0)
      break;
  }
  if (iDet == kNDetectors) 
    {
      TString errormessage;
      errormessage.Form("Invalid detector name: %s !",detectorName);
      Log(errormessage.Data());
      return -1;
    }

  Int_t detectorbitmask = 0;
  if(iDet > 32)
    {
      TString errormessage2;
      errormessage2.Form("Invalid detector bit position in detectorMask: %d !", iDet);
      Log(errormessage2.Data());
      return -1;
    }
  
  detectorbitmask = (1 << iDet);
  return detectorbitmask;
}

Bool_t AliHLTModulePreprocessor::GetDetectorStatus(Int_t detectorbitmask)
{
  // see header file for function documentation
  // retrieve list of active detectors from previous run.
  fActiveDetectors = atoi(AliHLTModulePreprocessor::GetRunParameter("detectorMask"));
 
  if((fActiveDetectors & detectorbitmask) != 0)
    {
      return 1;
    }
  else
    {
      return 0;
    }
}

// function copied from AliDAQ.cxx
const char *AliHLTModulePreprocessor::DetectorName(Int_t detectorID)
{
  // Returns the name of particular
  // detector identified by its index
  if (detectorID < 0 || detectorID >= kNDetectors) 
    {
      TString errormessage;
      errormessage.Form("Invalid detector index: %d (%d -> %d) !",detectorID,0,kNDetectors-1);
      Log(errormessage.Data());
      return "";
    }
  return fgkDetectorName[detectorID];
}

TObject* AliHLTModulePreprocessor::GetFromMap(TMap* dcsAliasMap, const char* stringId) const
{
  /// extract object from the alias map
  /// return the last object from the value set
  TObject* object=dcsAliasMap->FindObject(stringId);
  if (!object) return NULL;
  TPair* pair = dynamic_cast<TPair*>(object);
  if (pair && pair->Value()) {
    TObjArray* valueSet = dynamic_cast<TObjArray*>(pair->Value());
    if (!valueSet) return NULL;
    Int_t nentriesDCS = valueSet->GetEntriesFast() - 1;
    if(nentriesDCS>=0 && valueSet->At(nentriesDCS)){
      return valueSet->At(nentriesDCS);
    }
  }
  return NULL;
}
