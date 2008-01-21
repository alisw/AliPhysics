// $Id: AliHLTModulePreprocessor.cxx 23039 2007-12-13 20:53:02Z richterm $

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
#include "AliHLTPreprocessor.h"

ClassImp(AliHLTModulePreprocessor)

AliHLTModulePreprocessor::AliHLTModulePreprocessor() 
  :
  fpContainer(NULL)
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

void AliHLTModulePreprocessor::SetContainer(AliHLTPreprocessor* pContainer)
{
  assert(fpContainer==NULL || fpContainer==pContainer || pContainer==NULL);
  fpContainer=pContainer;
}

Int_t AliHLTModulePreprocessor::GetRun()
{
  // see header file for function documentation

  assert(fpContainer);
  if (!fpContainer) return 0;
  return fpContainer->GetRun();
}

UInt_t AliHLTModulePreprocessor::GetStartTime()
{
  // see header file for function documentation

  assert(fpContainer);
  if (!fpContainer) return 0;
  return fpContainer->GetStartTime();
}

UInt_t AliHLTModulePreprocessor::GetEndTime()
{
  // see header file for function documentation

  assert(fpContainer);
  if (!fpContainer) return 0;
  return fpContainer->GetEndTime();
}

Bool_t AliHLTModulePreprocessor::Store(const char* pathLevel2, const char* pathLevel3, TObject* object,
				 AliCDBMetaData* metaData, Int_t validityStart, Bool_t validityInfinite)
{
  // see header file for function documentation

  assert(fpContainer);
  if (!fpContainer) return 0;
  return fpContainer->Store(pathLevel2, pathLevel3, object, metaData, validityStart, validityInfinite);
}

Bool_t AliHLTModulePreprocessor::StoreReferenceData(const char* pathLevel2, const char* pathLevel3, TObject* object,
					      AliCDBMetaData* metaData)
{
  // see header file for function documentation

  assert(fpContainer);
  if (!fpContainer) return 0;
  return fpContainer->StoreReferenceData(pathLevel2, pathLevel3, object, metaData);
}

Bool_t AliHLTModulePreprocessor::StoreReferenceFile(const char* localFile, const char* gridFileName)
{
  // see header file for function documentation

  assert(fpContainer);
  if (!fpContainer) return 0;
  return fpContainer->StoreReferenceFile(localFile, gridFileName);
}

Bool_t AliHLTModulePreprocessor::StoreRunMetadataFile(const char* localFile, const char* gridFileName)
{
  // see header file for function documentation

  assert(fpContainer);
  if (!fpContainer) return 0;
  return fpContainer->StoreRunMetadataFile(localFile, gridFileName);
}
    
const char* AliHLTModulePreprocessor::GetFile(Int_t system, const char* id, const char* source)
{
  // see header file for function documentation

  assert(fpContainer);
  if (!fpContainer) return 0;
  return fpContainer->GetFile(system, id, source);
}

TList* AliHLTModulePreprocessor::GetFileSources(Int_t system, const char* id)
{
  // see header file for function documentation

  assert(fpContainer);
  if (!fpContainer) return 0;
  return fpContainer->GetFileSources(system, id);
}

TList* AliHLTModulePreprocessor::GetFileIDs(Int_t system, const char* source)
{
  // see header file for function documentation

  assert(fpContainer);
  if (!fpContainer) return 0;
  return fpContainer->GetFileIDs(system, source);
}

const char* AliHLTModulePreprocessor::GetRunParameter(const char* param)
{
  // see header file for function documentation

  assert(fpContainer);
  if (!fpContainer) return 0;
  return fpContainer->GetRunParameter(param);
}

AliCDBEntry* AliHLTModulePreprocessor::GetFromOCDB(const char* pathLevel2, const char* pathLevel3)
{
  // see header file for function documentation

  assert(fpContainer);
  if (!fpContainer) return 0;
  return fpContainer->GetFromOCDB(pathLevel2, pathLevel3);
}

const char* AliHLTModulePreprocessor::GetRunType()
{
  // see header file for function documentation

  assert(fpContainer);
  if (!fpContainer) return 0;
  return fpContainer->GetRunType();
}

void AliHLTModulePreprocessor::Log(const char* message)
{
  // see header file for function documentation

  assert(fpContainer);
  if (!fpContainer) return;
  fpContainer->Log(message);
}
