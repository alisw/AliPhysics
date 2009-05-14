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

/** @file   AliHLTBlockFilterComponent.cxx
    @author Matthias Richter
    @date   
    @brief  A simple data block filter and merger, merges block descriptors
*/

#include <cstdlib>
#include "AliHLTBlockFilterComponent.h"
#include "TString.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTBlockFilterComponent)

AliHLTBlockFilterComponent::AliHLTBlockFilterComponent()
  :
  fFilterRules()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTBlockFilterComponent::~AliHLTBlockFilterComponent()
{
  // see header file for class documentation
}

void AliHLTBlockFilterComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back(kAliHLTAnyDataType);
}

AliHLTComponentDataType AliHLTBlockFilterComponent::GetOutputDataType()
{
  // see header file for class documentation
  if (fFilterRules.size()==1) return fFilterRules[0].fDataType;
  if (fFilterRules.size()==0) return kAliHLTAnyDataType;
  return kAliHLTMultipleDataType;
}

int AliHLTBlockFilterComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
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

void AliHLTBlockFilterComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  constBase=0;
  inputMultiplier=0.0; // there is no new data, just forwarded descriptors
}

int AliHLTBlockFilterComponent::DoInit( int argc, const char** argv )
{
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
  if (iResult>=0 && (rule.fSpecification!=kAliHLTVoidDataSpec || !MatchExactly(rule.fDataType,kAliHLTAnyDataType))) {
    // add the pending rule
    fFilterRules.push_back(rule);
    FillBlockData(rule);
  }
  return iResult;
}

int AliHLTBlockFilterComponent::DoDeinit()
{
  int iResult=0;
  fFilterRules.clear();
  return iResult;
}

int AliHLTBlockFilterComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/,
					 const AliHLTComponentBlockData* /*blocks*/, 
					 AliHLTComponentTriggerData& /*trigData*/,
					 AliHLTUInt8_t* /*outputPtr*/, 
					 AliHLTUInt32_t& size,
					 AliHLTComponentBlockDataList& /*outputBlocks*/ )
{
  // see header file for class documentation
  int iResult=0;
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock();
       pBlock!=NULL; 
       pBlock=GetNextInputBlock()) {
    if (IsSelected(*pBlock)) {
      HLTDebug("block type %s %#x (ptr=%p offset=%d size=%d) selected by filter rules", 
	       DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, 
	       pBlock->fPtr, pBlock->fOffset, pBlock->fSize);
      Forward();
    } else {
      HLTDebug("block type %s %#x discarded by filter rules", DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification);
    }
  }
  size=0;
  return iResult;
}

int AliHLTBlockFilterComponent::IsSelected(const AliHLTComponentBlockData& block)
{
  // see header file for class documentation
  AliHLTComponentBlockDataList::iterator desc=fFilterRules.begin();
  //HLTDebug("check block: %s spec %#x", DataType2Text(block.fDataType, 1).c_str(), block.fSpecification);
  if (desc==fFilterRules.end()) return 1; // no filter rules
  do {
    // match if
    // 1. data types match or filter data type not set
    // 2. data spec match or filter data wpec not set
    // 3. either filter data type or spec is set
    //HLTDebug("check rule : %s spec %#x", DataType2Text((*desc).fDataType, 2).c_str(), block.fSpecification);
    if (((*desc).fDataType==block.fDataType) &&
	((*desc).fSpecification==block.fSpecification || (*desc).fSpecification==kAliHLTVoidDataSpec) &&
	(!MatchExactly((*desc).fDataType,kAliHLTAnyDataType) || (*desc).fSpecification!=kAliHLTVoidDataSpec)) {
      return 1;
    }
  } while (++desc!=fFilterRules.end());
  
  return 0;
}
