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

/** @file   AliHLTBlockDataCollection.cxx
    @author Matthias Richter
    @date   
    @brief  A collection of AliHLTComponentBlockData descriptors providing
            argument parsing and basic selection.
*/

#include <cstdlib>
#include "AliHLTBlockDataCollection.h"
#include "AliHLTComponent.h"
#include "TString.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTBlockDataCollection)

AliHLTBlockDataCollection::AliHLTBlockDataCollection()
  : TObject()
  , AliHLTLogging()
  , fFilterRules()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTBlockDataCollection::~AliHLTBlockDataCollection()
{
  // see header file for class documentation
}

int AliHLTBlockDataCollection::Add(const AliHLTComponentBlockData& block)
{
  // see header file for class documentation
  fFilterRules.push_back(block);
  return fFilterRules.size();
}

int AliHLTBlockDataCollection::ScanArgument( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;
  TString argument="";
  int bMissingParam=0;
  AliHLTComponentBlockData rule;
  AliHLTComponent::FillBlockData(rule);
  int i=0;
  for (; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

    // -datatype
    if (argument.CompareTo("-datatype")==0) {
      if ((bMissingParam=(i+2>=argc))) break;

      if (!MatchExactly(rule.fDataType,kAliHLTAnyDataType)) {
	// the data type has already been set, add to list
	// and reset
	fFilterRules.push_back(rule);
	AliHLTComponent::FillBlockData(rule);
      }

      AliHLTComponent::SetDataType(rule.fDataType, argv[i+1], argv[i+2]);
      i+=2;

      // -origin
    } else if (argument.CompareTo("-origin")==0) {
      if ((bMissingParam=(i+1>=argc))) break;

      if (!MatchExactly(rule.fDataType,kAliHLTAnyDataType)) {
	// the data type has already been set, add to list
	// and reset
	fFilterRules.push_back(rule);
	AliHLTComponent::FillBlockData(rule);
      }

      AliHLTComponent::SetDataType(rule.fDataType, NULL, argv[i+1]);
      i+=1;

      // -typeid
    } else if (argument.CompareTo("-typeid")==0) {
      if ((bMissingParam=(i+1>=argc))) break;

      if (!MatchExactly(rule.fDataType,kAliHLTAnyDataType)) {
	// the data type has already been set, add to list
	// and reset
	fFilterRules.push_back(rule);
	AliHLTComponent::FillBlockData(rule);
      }

      AliHLTComponent::SetDataType(rule.fDataType, argv[i+1], NULL);
      i+=1;

      // -dataspec
    } else if (argument.CompareTo("-dataspec")==0) {
      if ((bMissingParam=(++i>=argc))) break;

      if (rule.fSpecification!=kAliHLTVoidDataSpec) {
	// the specification has already been set, add to list
	// and reset
	fFilterRules.push_back(rule);
	AliHLTComponent::FillBlockData(rule);
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
      // terminate at the first unknown argument
      break;
    }
  }

  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EPROTO;
  }
  if (iResult>=0) {
    if (rule.fSpecification!=kAliHLTVoidDataSpec || !MatchExactly(rule.fDataType,kAliHLTAnyDataType)) {
      // add the pending rule
      fFilterRules.push_back(rule);
      AliHLTComponent::FillBlockData(rule);
    }
    iResult=i;
  }

  return iResult;
}

int AliHLTBlockDataCollection::IsSelected(const AliHLTComponentBlockData& block)
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

int AliHLTBlockDataCollection::IsEmpty()
{
  // see header file for class documentation
  if (fFilterRules.size()==0) return 1;
  return 0;
}
