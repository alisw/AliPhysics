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

/**
 * @file   AliHLTModulePreprocessor.cxx
 * @author Matthias Richter
 * @date   2008-01-22
 * @brief  Base class for HLT module preprocessors
 */

#include <cassert>
#include "AliHLTModulePreprocessor.h"
#include "AliHLTShuttleInterface.h"

ClassImp(AliHLTModulePreprocessor)

AliHLTModulePreprocessor::AliHLTModulePreprocessor() 
  :
  fpInterface(NULL)
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

void AliHLTModulePreprocessor::SetShuttleInterface(AliHLTShuttleInterface* pInterface)
{
  assert(fpInterface==NULL || fpInterface==pInterface || pInterface==NULL);
  fpInterface=pInterface;
}

Int_t AliHLTModulePreprocessor::GetRun()
{
  // see header file for function documentation

  assert(fpInterface);
  if (!fpInterface) return 0;
  return fpInterface->PreprocessorGetRun();
}

UInt_t AliHLTModulePreprocessor::GetStartTime()
{
  // see header file for function documentation

  assert(fpInterface);
  if (!fpInterface) return 0;
  return fpInterface->PreprocessorGetStartTime();
}

UInt_t AliHLTModulePreprocessor::GetEndTime()
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
