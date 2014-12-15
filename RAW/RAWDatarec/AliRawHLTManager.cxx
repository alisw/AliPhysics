// $Id: AliRawHLTManager.cxx 23039 2007-12-13 20:53:02Z richterm $

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

/** @file   AliRawHLTManager.cxx
    @author Matthias Richter
    @date   
    @brief  dynamic generation of HLT RAW readers and streams
*/

#include "AliRawHLTManager.h"
#include "AliLog.h"
#include "TSystem.h"
#include "TClass.h"

ClassImp(AliRawHLTManager)

AliRawHLTManager::AliRawHLTManager()
{
  // The class gives dynamic access to creater methods for HLT RAW readers and
  // streams without any library dependencies to HLT libraries.
  //
  // The AliRawReaderHLT allows the redirection of input from the HLT DDL links
  // to the detector equipment ids. To access the data, the AliRawReaderHLT
  // needs a valid RAW reader (parent).
  // usage:
  // AliRawReader* pHLTReader=AliRawHLTManager::CreateRawReaderHLT(pParent, "TPC");
}

AliRawHLTManager::~AliRawHLTManager()
{
  // destructor
}

int AliRawHLTManager::fLibraryStatus=kUnloaded;
AliRawReaderHLTCreateInstance_t AliRawHLTManager::fFctCreateRawReaderHLT=NULL;
void* AliRawHLTManager::fFctCreateRawStream=NULL;

AliRawReader* AliRawHLTManager::CreateRawReaderHLT(AliRawReader* pParent, const char* detectors)
{
  // see header file for class documentation
  if (fLibraryStatus==kUnloaded) fLibraryStatus=LoadLibrary();
  if (fLibraryStatus==kUnavailable) return NULL;

  if (!fFctCreateRawReaderHLT) {
    AliErrorClass("internal error, library loaded but function entry not known");
    return NULL;
  }
  return ((AliRawReaderHLTCreateInstance_t)fFctCreateRawReaderHLT)(pParent, detectors);
}

TObject* AliRawHLTManager::CreateRawStream(const char* /*className*/)
{
  // see header file for class documentation

  // not yet implemented
  return NULL;
}

int AliRawHLTManager::LoadLibrary()
{
  // see header file for class documentation
  int iResult=kUnavailable;
  if (fLibraryStatus!=kUnloaded) return fLibraryStatus;

  // strictly speaken we do not need a trial counter as gSystem->Load only returns 1 if the
  // library has been loaded. If it was already loaded we get 0 
  int iTrials=0;
  do {
    fFctCreateRawReaderHLT=(AliRawReaderHLTCreateInstance_t)gSystem->DynFindSymbol(ALIHLTREC_LIBRARY, ALIRAWREADERHLT_CREATE_INSTANCE);
  } while (fFctCreateRawReaderHLT==NULL && gSystem->Load(ALIHLTREC_LIBRARY)==0 && iTrials++<1);
  if (fFctCreateRawReaderHLT) {
    iResult=kLoaded;
  } else {
    AliErrorClass(Form("can not find library/entry %s/%s", ALIHLTREC_LIBRARY, ALIRAWREADERHLT_CREATE_INSTANCE));
    iResult=kUnavailable;
  }
  return iResult;
}
