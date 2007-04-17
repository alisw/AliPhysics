// $Id$

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Matthias Richter <Matthias.Richter@ift.uib.no>                *
 *          for The ALICE Off-line Project.                               *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTOfflineInterface.cxx
    @author Matthias Richter
    @date   
    @brief  the HLT interface to AliRoot
*/

#include "AliHLTOfflineInterface.h"
#include "AliHLTLogging.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOfflineInterface)

AliHLTOfflineInterface::AliHLTOfflineInterface()
  :
  fpRunLoader(NULL),
  fpRawReader(NULL),
  fpESD(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

TList AliHLTOfflineInterface::fgList;
TObjLink* AliHLTOfflineInterface::fgCurrentLnk=NULL;
AliRunLoader* AliHLTOfflineInterface::fgpRunLoader=NULL;
AliRawReader* AliHLTOfflineInterface::fgpRawReader=NULL;

AliHLTOfflineInterface::AliHLTOfflineInterface(AliRunLoader* pRunLoader, AliRawReader* pRawReader)
  :
  fpRunLoader(pRunLoader),
  fpRawReader(pRawReader),
  fpESD(NULL)
{
}

AliHLTOfflineInterface::AliHLTOfflineInterface(const AliHLTOfflineInterface&)
  :
  TObject(),
  fpRunLoader(NULL),
  fpRawReader(NULL),
  fpESD(NULL)
{
  // see header file for class documentation
  //HLTFatal("copy constructor untested");
}

AliHLTOfflineInterface& AliHLTOfflineInterface::operator=(const AliHLTOfflineInterface&)
{ 
  // see header file for class documentation
  //HLTFatal("assignment operator untested");
  return *this;
}

AliHLTOfflineInterface::~AliHLTOfflineInterface()
{
}

AliRunLoader* AliHLTOfflineInterface::GetRunLoader() const
{
  return fpRunLoader!=NULL?fpRunLoader:fgpRunLoader;
}

AliRawReader* AliHLTOfflineInterface::GetRawReader() const
{
  return fpRawReader!=NULL?fpRawReader:fgpRawReader;
}

int AliHLTOfflineInterface::SetESD(Int_t eventNo, AliESD* pESD)
{
  fpESD=pESD;
  return 0;
}

AliESD* AliHLTOfflineInterface::GetESD() const
{
  return fpESD;
}

int AliHLTOfflineInterface::SetParams(AliRunLoader* runLoader, AliRawReader* rawReader)
{
  // see header file for class documentation
  int iResult=0;
  if (fpRunLoader!=NULL && fpRunLoader!=runLoader) {
    //HLTWarning("overriding previous instance of Run Loader %p with %p", fpRunLoader, runLoader);
  }
  fpRunLoader=runLoader;
  if (fpRawReader!=NULL && fpRawReader!=rawReader) {
    //HLTWarning("overriding previous instance of RawReader %p with %p", fpRawReader, rawReader);
  }
  fpRawReader=rawReader;
  return iResult;
}

int AliHLTOfflineInterface::Reset()
{
  // see header file for class documentation
  int iResult=0;
  fpRunLoader=NULL;
  fpRawReader=NULL;
  fpESD=NULL;
  return iResult;
}

int AliHLTOfflineInterface::SetParamsToComponents(AliRunLoader* runLoader, AliRawReader* rawReader)
{
  // see header file for class documentation
  AliHLTLogging log;
  int iResult=0;
  int count=0;
  fgpRunLoader=runLoader;
  fgpRawReader=rawReader;
  fgCurrentLnk=fgList.FirstLink();
  while (fgCurrentLnk!=NULL) {
    AliHLTOfflineInterface* pComponent=reinterpret_cast<AliHLTOfflineInterface*>(fgCurrentLnk->GetObject());
    int iLocal=0;
    if (pComponent) iLocal=pComponent->SetParams(runLoader, rawReader);
    if (iLocal<0) {
      log.LoggingVarargs(kHLTLogWarning, "AliHLTOfflineInterface","SetParamsToComponents", __FILE__, __LINE__,
			"parameter initialization failed for component %p with result %d", pComponent, iLocal);
      if (iResult>=0) iResult=iLocal;
    }
    count++;
    fgCurrentLnk=fgCurrentLnk->Next();
  }
  if (iResult>=0) {
//       log.LoggingVarargs(kHLTLogInfo, "AliHLTOfflineInterface","SetParamsToComponents", __FILE__, __LINE__,
// 			"parameters set to %d offline interface component(s)", count);
  }
  return iResult;
}

int AliHLTOfflineInterface::ResetComponents()
{
  // see header file for class documentation
  int iResult=0;
  fgCurrentLnk=fgList.FirstLink();
  while (fgCurrentLnk!=NULL) {
    AliHLTOfflineInterface* pComponent=reinterpret_cast<AliHLTOfflineInterface*>(fgCurrentLnk->GetObject());
    int iLocal=0;
    if (pComponent) iLocal=pComponent->Reset();
    if (iLocal<0) {
      if (iResult>=0) iResult=iLocal;
    }
    fgCurrentLnk=fgCurrentLnk->Next();
  }
  return iResult;
}

int AliHLTOfflineInterface::FillComponentESDs(int eventNo, AliRunLoader* runLoader, AliESD* esd)
{
  // see header file for class documentation
  int iResult=0;
  fgCurrentLnk=fgList.FirstLink();
  while (fgCurrentLnk!=NULL) {
    AliHLTOfflineInterface* pComponent=reinterpret_cast<AliHLTOfflineInterface*>(fgCurrentLnk->GetObject());
    int iLocal=0;
    if (pComponent) {
      pComponent->SetESD(eventNo, esd);
      if (pComponent->GetRunLoader()!=runLoader) {
	//HLTWarning("runLoader missmatch: component %p was reconstructed with runLoader %p, but got %p now", pComponent, pComponent->GetRunLoader(), runLoader);
      }
      iLocal=pComponent->FillESD(eventNo, runLoader, esd);
    }
    if (iLocal<0) {
      if (iResult>=0) iResult=iLocal;
    }
    fgCurrentLnk=fgCurrentLnk->Next();
  }
  return iResult;
}

int AliHLTOfflineInterface::Register(AliHLTOfflineInterface* me)
{
  // see header file for function documentation
  int iResult=0;
  if (fgList.FindObject(me->GetName())==NULL) {
    fgList.Add(me);
  } else {
  }
  return iResult;
}

int AliHLTOfflineInterface::Unregister(AliHLTOfflineInterface* me)
{
  // see header file for function documentation
  int iResult=0;
  if (fgList.FindObject(me)!=NULL) {
    fgList.Remove(me);
  } else {
  }
  return iResult;
}
