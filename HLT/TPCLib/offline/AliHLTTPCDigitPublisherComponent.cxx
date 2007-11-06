// @(#) $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTPCDigitPublisherComponent.cxx
    @author Matthias Richter
    @date   
    @brief  TPC digit publisher component (input from offline).
*/

#include "AliHLTTPCDigitPublisherComponent.h"
#include "AliRunLoader.h"
#include "AliLog.h"
#include "TTree.h"
#include "AliHLTTPCDigitData.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCFileHandler.h"

/** global instance for agent registration */
AliHLTTPCDigitPublisherComponent gAliHLTTPCDigitPublisherComponent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCDigitPublisherComponent)

AliHLTTPCDigitPublisherComponent::AliHLTTPCDigitPublisherComponent()
  :
  fMaxSize(200000), // just a number to start with
  fMinSlice(-1),
  fMinPart(-1)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCFileHandler* AliHLTTPCDigitPublisherComponent::fpFileHandler=NULL;
int AliHLTTPCDigitPublisherComponent::fFileHandlerInstances=0;
int AliHLTTPCDigitPublisherComponent::fCurrEvent=-1;

AliHLTTPCDigitPublisherComponent::~AliHLTTPCDigitPublisherComponent()
{
  // see header file for class documentation
  if (fpFileHandler!=NULL && fFileHandlerInstances<=0) {
    HLTWarning("improper state, de-initialization missing");
    DoDeinit();
  }
}

const char* AliHLTTPCDigitPublisherComponent::GetComponentID()
{
  // see header file for class documentation
  return "TPCDigitPublisher";
}

AliHLTComponentDataType AliHLTTPCDigitPublisherComponent::GetOutputDataType()
{
  return AliHLTTPCDefinitions::fgkUnpackedRawDataType;
}

void AliHLTTPCDigitPublisherComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  constBase=fMaxSize;
  inputMultiplier=1;
}

AliHLTComponent* AliHLTTPCDigitPublisherComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTPCDigitPublisherComponent;
}

int AliHLTTPCDigitPublisherComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;

  // scan arguments
  TString argument="";
  int bMissingParam=0;
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

    // -slice
    if (argument.CompareTo("-slice")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      TString parameter(argv[i]);
      parameter.Remove(TString::kLeading, ' '); // remove all blanks
      if (parameter.IsDigit()) {
	fMinSlice=parameter.Atoi();
      } else {
	HLTError("wrong parameter for argument %s, number expected", argument.Data());
	iResult=-EINVAL;
      }
      // -partition
    } else if (argument.CompareTo("-partition")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      TString parameter(argv[i]);
      parameter.Remove(TString::kLeading, ' '); // remove all blanks
      if (parameter.IsDigit()) {
	fMinPart=parameter.Atoi();
      } else {
	HLTError("wrong parameter for argument %s, number expected", argument.Data());
	iResult=-EINVAL;
      }
    } else {
      HLTError("unknown argument %s", argument.Data());
      iResult=-EINVAL;
    }
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }

  if (fMinSlice<0) {
    HLTError("slice no required");
    iResult=-EINVAL;
  }

  if (fMinPart<0) {
    HLTError("partition (patch) no required");
    iResult=-EINVAL;
  }

  if (iResult<0) return iResult;

  // fetch runLoader instance from interface
  AliRunLoader* pRunLoader=GetRunLoader();
  if (pRunLoader) {
    if (fpFileHandler==NULL) {
      fpFileHandler=new AliHLTTPCFileHandler;
      fCurrEvent=-1;
      fFileHandlerInstances=1;
    } else {
      fFileHandlerInstances++;
      //HLTDebug("publisher %p: %d references to file handler instance", this, fFileHandlerInstances);
    }
    if (fpFileHandler) {
      if (!fpFileHandler->SetAliInput(pRunLoader)) {
	iResult=-EFAULT;
      }
    } else {
      AliErrorStream() << "can not allocate file handler object" << endl;
      iResult=-ENOMEM;
    }
  } else {
    AliErrorStream() << "can not get runLoader" << endl;
    iResult=-EFAULT;
  }
  return iResult;
}

int AliHLTTPCDigitPublisherComponent::DoDeinit()
{
  // see header file for class documentation
  int iResult=0;
  if (fCurrEvent>=0) {
    fpFileHandler->FreeDigitsTree();
    fCurrEvent=-1;
  }
  //HLTDebug("publisher %p: %d references to file handler instance", this, fFileHandlerInstances);
  if (--fFileHandlerInstances==0 && fpFileHandler!=NULL) {
    try {
      if (fpFileHandler) {
	delete fpFileHandler;
      }
    }
    catch (...) {
      HLTFatal("exeption during object cleanup");
      iResult=-EFAULT;
    }
    fpFileHandler=NULL;
  }
  return iResult;
}

int AliHLTTPCDigitPublisherComponent::GetEvent(const AliHLTComponentEventData& evtData,
					       AliHLTComponentTriggerData& trigData,
					       AliHLTUInt8_t* outputPtr, 
					       AliHLTUInt32_t& size,
					       vector<AliHLTComponentBlockData>& outputBlocks)
{
  // see header file for class documentation
  int iResult=0;
  if (outputPtr==NULL || size==0) {
    HLTError("no target buffer provided");
    return -EFAULT;
  }

  if (fpFileHandler) {
    int event=GetEventCount();
    AliHLTTPCUnpackedRawData* pTgt=reinterpret_cast<AliHLTTPCUnpackedRawData*>(outputPtr);
    if (pTgt) {
      UInt_t nrow=0;
      UInt_t tgtSize=size-sizeof(AliHLTTPCUnpackedRawData);
      if (fCurrEvent>=0 && fCurrEvent!=event) {
	HLTDebug("new event %d, free digit tree for event %d", event, fCurrEvent);
	fpFileHandler->FreeDigitsTree();
      }
      fCurrEvent=event;
      HLTDebug("converting digits for slice %d partition %d", fMinSlice, fMinPart);
      fpFileHandler->Init(fMinSlice,fMinPart);
      AliHLTTPCDigitRowData* pData=fpFileHandler->AliDigits2Memory(nrow, event, reinterpret_cast<Byte_t*>(pTgt->fDigits), &tgtSize);
      if (pData==NULL && tgtSize>0 && tgtSize>fMaxSize) {
	HLTDebug("target buffer too small: %d byte required, %d available", tgtSize+sizeof(AliHLTTPCUnpackedRawData), size);
	// indicate insufficient buffer size, on occasion the frameworks calls
	// again with the corrected buffer 
	fMaxSize=tgtSize;
	iResult=-ENOSPC;
      } else if (pData!=pTgt->fDigits) {
	HLTError("can not read directly into output buffer");
	try {delete pData;}
	catch (...) {/* no action */}
	iResult=-EIO;
      } else {
	size=tgtSize+sizeof(AliHLTTPCUnpackedRawData);
	AliHLTComponentBlockData bd;
	FillBlockData( bd );
	bd.fOffset = 0;
	bd.fSize = size;
	bd.fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification(fMinSlice, fMinSlice, fMinPart, fMinPart);
	outputBlocks.push_back( bd );
	HLTDebug("added AliHLTTPCUnpackedRawData size %d, first row %d nof digits %d", size, pTgt->fDigits->fRow, pTgt->fDigits->fNDigit);
      }
    }
  } else {
    AliErrorStream() << "component not initialized" << endl;
    iResult=-EFAULT;
  }
  return iResult;
}
