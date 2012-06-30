// $Id$

//**************************************************************************
//* This file is property of and copyright by the                          * 
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
//**************************************************************************/

/// @file   AliHLTOfflineInterface.cxx
/// @author Matthias Richter
/// @date   
/// @brief  the HLT interface to AliRoot
///

#include "AliHLTOfflineInterface.h"
#include "AliHLTLogging.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOfflineInterface)

AliHLTOfflineInterface::AliHLTOfflineInterface()
  :
  fpRunLoader(NULL),
  fpRawReader(NULL),
  fpESD(NULL),
  fpNext(NULL)
{
  // The class implements the basic interface to the AliRoot objects during
  // reconstructions.
  // It serves as a base class for offline source and sink interface components
  // and provides access methods for the AliRunLoader, AliRawReader and AliESDEvent
  // objects. The AliRunLoader and the AliRawReader are fixed during one run,
  // while the AliESDEvent object will be changed from event to event.<br>
  // \em Note: The digits and clusters trees are not available through this
  // interface class as they are completetly detector (AliLoader) dependend.
}

AliHLTOfflineInterface* AliHLTOfflineInterface::fgAnchor=NULL;
AliHLTOfflineInterface* AliHLTOfflineInterface::fgCurrent=NULL;
int AliHLTOfflineInterface::fgCount=0;
AliRunLoader* AliHLTOfflineInterface::fgpRunLoader=NULL;
AliRawReader* AliHLTOfflineInterface::fgpRawReader=NULL;

AliHLTOfflineInterface::AliHLTOfflineInterface(AliRunLoader* pRunLoader, AliRawReader* pRawReader)
  :
  fpRunLoader(pRunLoader),
  fpRawReader(pRawReader),
  fpESD(NULL),
  fpNext(NULL)
{
  // constructor
}

AliHLTOfflineInterface::~AliHLTOfflineInterface()
{
  // destructor
}

AliRunLoader* AliHLTOfflineInterface::GetRunLoader() const
{
  // set RawReader pointer
  return fpRunLoader!=NULL?fpRunLoader:fgpRunLoader;
}

AliRawReader* AliHLTOfflineInterface::GetRawReader() const
{
  // get RawReader pointer
  return fpRawReader!=NULL?fpRawReader:fgpRawReader;
}

int AliHLTOfflineInterface::SetESD(Int_t /*eventNo*/, AliESDEvent* pESD)
{
  // set ESD pointer
  fpESD=pESD;
  return 0;
}

AliESDEvent* AliHLTOfflineInterface::GetESD() const
{
  // get ESD pointer
  return fpESD;
}

int AliHLTOfflineInterface::SetParams(AliRunLoader* runLoader, AliRawReader* rawReader)
{
  // set parameters of the interface
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
  // pass parameters to registered components
  AliHLTLogging log;
  int iResult=0;
  int count=0;
  fgpRunLoader=runLoader;
  fgpRawReader=rawReader;
  AliHLTOfflineInterface* pCurrent=fgAnchor;
  while (pCurrent!=NULL) {
    int iLocal=0;
    if (pCurrent) iLocal=pCurrent->SetParams(runLoader, rawReader);
    if (iLocal<0) {
      log.LoggingVarargs(kHLTLogWarning, "AliHLTOfflineInterface","SetParamsToComponents", __FILE__, __LINE__,
			"parameter initialization failed for component %p with result %d", pCurrent, iLocal);
      if (iResult>=0) iResult=iLocal;
    }
    count++;
    pCurrent=pCurrent->fpNext;
  }
  if (iResult>=0) {
//       log.LoggingVarargs(kHLTLogInfo, "AliHLTOfflineInterface","SetParamsToComponents", __FILE__, __LINE__,
// 			"parameters set to %d offline interface component(s)", count);
  }
  return iResult;
}

int AliHLTOfflineInterface::ResetComponents()
{
  // loop over registered components and call Reset
  int iResult=0;
  AliHLTOfflineInterface* pCurrent=fgAnchor;
  while (pCurrent!=NULL) {
    int iLocal=0;
    if (pCurrent) iLocal=pCurrent->Reset();
    if (iLocal<0) {
      if (iResult>=0) iResult=iLocal;
    }
    pCurrent=pCurrent->fpNext;
  }
  return iResult;
}

int AliHLTOfflineInterface::FillComponentESDs(int eventNo, AliRunLoader* runLoader, AliESDEvent* esd)
{
  // loop ove registered components and call FillESD function
  int iResult=0;
  AliHLTOfflineInterface* pCurrent=fgAnchor;
  while (pCurrent!=NULL) {
    int iLocal=0;
    if (pCurrent) {
      pCurrent->SetESD(eventNo, esd);
      if (pCurrent->GetRunLoader()!=runLoader) {
	//HLTWarning("runLoader mismatch: component %p was reconstructed with runLoader %p, but got %p now", pCurrent, pCurrent->GetRunLoader(), runLoader);
      }
      iLocal=pCurrent->FillESD(eventNo, runLoader, esd);
    }
    if (iLocal<0) {
      if (iResult>=0) iResult=iLocal;
    }
    pCurrent=pCurrent->fpNext;
  }
  return iResult;
}

int AliHLTOfflineInterface::Register(AliHLTOfflineInterface* me)
{
  // register a component in the list of offline interface components
  int iResult=0;
  if (fgAnchor==NULL) {
    fgAnchor=me;
  } else {
    me->fpNext=fgAnchor;
    fgAnchor=me;
  }
  fgCount++;
  return iResult;
}

int AliHLTOfflineInterface::Unregister(AliHLTOfflineInterface* me)
{
  // remove a component from the list
  int iResult=0;
  fgCurrent=NULL;
  AliHLTOfflineInterface* prev=NULL;
  AliHLTOfflineInterface* handler=fgAnchor;
  while (handler!=NULL && handler!=me) {
    prev=handler;
    handler=handler->fpNext;
  }
  if (handler) {
    if (prev==NULL) {
      fgAnchor=handler->fpNext;
    } else {
      prev->fpNext=handler->fpNext;
    }
    fgCount--;
  }
  return iResult;
}
