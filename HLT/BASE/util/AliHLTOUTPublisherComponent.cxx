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

/** @file   AliHLTOUTPublisherComponent.cxx
    @author Matthias Richter
    @date   2008-06-11
    @brief  A data publisher for data block out of the HLTOUT data
*/

#include <cstdlib>
#include "AliHLTOUTPublisherComponent.h"
#include "AliHLTOUT.h"
#include "TString.h"
#include "AliRawReader.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOUTPublisherComponent)

AliHLTOUTPublisherComponent::AliHLTOUTPublisherComponent()
  :
  fFilterRules(),
  fMaxSize(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTOUTPublisherComponent::~AliHLTOUTPublisherComponent()
{
  // see header file for class documentation
}

const char* AliHLTOUTPublisherComponent::GetComponentID()
{
  // see header file for class documentation
  return "AliHLTOUTPublisher";
}

AliHLTComponentDataType AliHLTOUTPublisherComponent::GetOutputDataType()
{
  // see header file for class documentation
  if (fFilterRules.size()==1) return fFilterRules[0].fDataType;
  if (fFilterRules.size()==0) return kAliHLTAnyDataType;
  return kAliHLTMultipleDataType;
}

int AliHLTOUTPublisherComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  // see header file for class documentation
  tgtList.clear();
  AliHLTComponentBlockDataList::iterator desc=fFilterRules.begin();
  while (desc!=fFilterRules.end()) {
    AliHLTComponentDataTypeList::iterator type=tgtList.begin();
    while (type!=tgtList.end()) {
      if (*type==(*desc).fDataType) break;
      type++;
    }
    if (type==tgtList.end()) tgtList.push_back((*desc).fDataType);
    desc++;
  }
  return tgtList.size();
}

void AliHLTOUTPublisherComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  constBase=fMaxSize;
  inputMultiplier=0.0; // there is no new data, just forwarded descriptors
}

AliHLTComponent* AliHLTOUTPublisherComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTOUTPublisherComponent;
}

int AliHLTOUTPublisherComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;
  TString argument="";
  int bMissingParam=0;
  AliHLTComponentBlockData rule;
  FillBlockData(rule);
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

    // -datatype
    if (argument.CompareTo("-datatype")==0) {
      if ((bMissingParam=(i+2>=argc))) break;

      if (!MatchExactly(rule.fDataType,kAliHLTAnyDataType)) {
	// the data type has already been set, add to list
	// and reset
	fFilterRules.push_back(rule);
	FillBlockData(rule);
      }

      SetDataType(rule.fDataType, argv[i+1], argv[i+2]);
      i+=2;

      // -origin
    } else if (argument.CompareTo("-origin")==0) {
      if ((bMissingParam=(i+1>=argc))) break;

      if (!MatchExactly(rule.fDataType,kAliHLTAnyDataType)) {
	// the data type has already been set, add to list
	// and reset
	fFilterRules.push_back(rule);
	FillBlockData(rule);
      }

      SetDataType(rule.fDataType, NULL, argv[i+1]);
      i+=1;

      // -typeid
    } else if (argument.CompareTo("-typeid")==0) {
      if ((bMissingParam=(i+1>=argc))) break;

      if (!MatchExactly(rule.fDataType,kAliHLTAnyDataType)) {
	// the data type has already been set, add to list
	// and reset
	fFilterRules.push_back(rule);
	FillBlockData(rule);
      }

      SetDataType(rule.fDataType, argv[i+1], NULL);
      i+=1;

      // -dataspec
    } else if (argument.CompareTo("-dataspec")==0) {
      if ((bMissingParam=(++i>=argc))) break;

      if (rule.fSpecification!=kAliHLTVoidDataSpec) {
	// the specification has already been set, add to list
	// and reset
	fFilterRules.push_back(rule);
	FillBlockData(rule);
      }

      TString parameter(argv[i]);
      parameter.Remove(TString::kLeading, ' '); // remove all blanks
      char* pRemnant=NULL;
      rule.fSpecification=strtoul(parameter.Data(), &pRemnant, 0);
      if (pRemnant!=NULL && pRemnant[0]!=0) {
	HLTError("invalid parameter/remnant (%s) for argument %s, number expected", pRemnant, argument.Data());
	iResult=-EINVAL;
      }
    } else {
      HLTError("unknown argument %s", argument.Data());
      iResult=-EINVAL;
      break;
    }
  }
  if (iResult>=0) {
    // add the pending rule or at least the empty default rule
    fFilterRules.push_back(rule);
    FillBlockData(rule);
  }
  return iResult;
}

int AliHLTOUTPublisherComponent::DoDeinit()
{
  // see header file for class documentation
  int iResult=0;
  fFilterRules.clear();
  return iResult;
}

int AliHLTOUTPublisherComponent::GetEvent( const AliHLTComponentEventData& /*evtData*/,
					   AliHLTComponentTriggerData& /*trigData*/,
					   AliHLTUInt8_t* outputPtr, 
					   AliHLTUInt32_t& size,
					   AliHLTComponentBlockDataList& outputBlocks )
{
  // see header file for class documentation
  int iResult=0;

  // process data events only
  if (!IsDataEvent()) return 0;

  unsigned int offset=0;
  AliHLTOUT* pHLTOUT=NULL;
  AliRawReader* pRawReader=GetRawReader();
  if ((pHLTOUT=AliHLTOUT::GetGlobalInstance())!=NULL) {
    // this is the HLTOUT instance set globally by the AliHLTOUT::AliHLTOUTGlobalInstanceGuard
    // used for data input from the HLTOUT to the publishers of a kChain handler
  } else
  if (pRawReader) {
    pRawReader->Reset();
    pHLTOUT=AliHLTOUT::New(pRawReader);
    if (pHLTOUT) iResult=pHLTOUT->Init();
//   } else {
    // this is just a hack and work-around for the missing HLT AliLoader.
    // Because of that the AliRoot framework does not provide the digit tree.
    // The HLTDigits.root file is opened directly in the AliHLTOUTDigitReader.
    // Later, the TTree digit tree object must be fetched here and passed to
    // the HLTOUT instance.
    // Maybe it's obsolete anyhow: 
    // 1. AliHLT reconstruction from digit data is not supported
    // 2. When integrated into HLTOUHandler of type kChain, most likely the
    //    AliRawReaderMemory will be used.
    // The functionality is tested also with the AliHLTOUTDigitReader, for the
    // mentioned reasons I comment this branch.
//     pHLTOUT=AliHLTOUT::New(NULL, GetEventCount());
  } else {
    if (GetEventCount()==0) {
      HLTFatal("can not get RunLoader or RawReader, event processing aborted for current and subsequent events");
    }
    iResult=-ENODEV;
  }
  if (iResult>=0 && pHLTOUT) {
    if (true) { // condition was deprecated but keep for the sake of diff
      for (iResult=pHLTOUT->SelectFirstDataBlock();
	   iResult>=0;
	   iResult=pHLTOUT->SelectNextDataBlock()) {
	AliHLTComponentDataType dt=kAliHLTVoidDataType;
	AliHLTUInt32_t spec=kAliHLTVoidDataSpec;
	pHLTOUT->GetDataBlockDescription(dt, spec);
	if (fFilterRules.size()>0) {
	  // check if the block is selected
	  unsigned int rule=0;
	  for (; rule<fFilterRules.size(); rule++) {
	    if (fFilterRules[rule].fDataType!=dt) continue;
	    if (fFilterRules[rule].fSpecification!=kAliHLTVoidDataSpec &&
		fFilterRules[rule].fSpecification!=spec) continue;
	    break;
	  }
	  // skip the block if none of the filter rules matches
	  if (rule>=fFilterRules.size()) continue;
	}
	const AliHLTUInt8_t* pBuffer=NULL;
	AliHLTUInt32_t bufferSize=0;
	if ((iResult=pHLTOUT->GetDataBuffer(pBuffer, bufferSize))>=0) {
	  if (bufferSize+offset<=size) {
	    memcpy(outputPtr+offset, pBuffer, bufferSize);
	    AliHLTComponentBlockData bd;
	    FillBlockData( bd );
	    bd.fOffset = offset;
	    bd.fSize = bufferSize;
	    bd.fDataType = dt;
	    bd.fSpecification = spec;
	    outputBlocks.push_back( bd );
	  } else {
	    // we keep the loop going in order to collect the full size
	    fMaxSize=offset+bufferSize;
	  }
	  offset+=bufferSize;
	}
      }
      // -ENOENT is not an error but the return value for 'no more data block'
      if (iResult==-ENOENT) iResult=0;

      // indicate too little space in buffer for repeated processing
      if (offset>size) {
	iResult=-ENOSPC;
      }
    } else if (GetEventCount()<5) {
      const char* message="";
      if (GetEventCount()==4) message=", suppressing further messages";
      HLTError("failed initializing HLTOUT%s", message);
    }
    AliHLTOUT::Delete(pHLTOUT);
    pHLTOUT=NULL;
  } else {
    if (GetEventCount()==0) {
      HLTFatal("can not create HLTOUT instance, event processing aborted for current and most likely subsequent events");
    }
    iResult=-ENODEV;
  }

  // finally set the output size
  if (iResult>=0)
    size=offset;
  else
    size=0;

  return iResult;
}
