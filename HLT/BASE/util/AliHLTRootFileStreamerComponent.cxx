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

/** @file   AliHLTRootFileStreamerComponent.cxx
    @author Matthias Richter
    @date   
    @brief  Save objects in a ROOT memory file

                                                                          */

#include "AliHLTRootFileStreamerComponent.h"
#include "TString.h"

/** the global object for component registration */
AliHLTRootFileStreamerComponent gAliHLTRootFileStreamerComponent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTRootFileStreamerComponent)

AliHLTRootFileStreamerComponent::AliHLTRootFileStreamerComponent()
  :
  AliHLTProcessor(),
  fDataType(kAliHLTVoidDataType),
  fSpecification(~(AliHLTUInt32_t)0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

AliHLTRootFileStreamerComponent::AliHLTRootFileStreamerComponent(const AliHLTRootFileStreamerComponent&)
  :
  AliHLTProcessor(),
  fDataType(kAliHLTVoidDataType),
  fSpecification(~(AliHLTUInt32_t)0)
{
  // see header file for class documentation
}

AliHLTRootFileStreamerComponent& AliHLTRootFileStreamerComponent::operator=(const AliHLTRootFileStreamerComponent&)
{
  // see header file for class documentation
  return *this;
}

AliHLTRootFileStreamerComponent::~AliHLTRootFileStreamerComponent()
{
  // see header file for class documentation
}

void AliHLTRootFileStreamerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back(kAliHLTAnyDataType);
}

AliHLTComponentDataType AliHLTRootFileStreamerComponent::GetOutputDataType()
{
  // see header file for class documentation
  return fDataType;
}

void AliHLTRootFileStreamerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  constBase=500;
  inputMultiplier=5.0;
}

int AliHLTRootFileStreamerComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation

  int iResult=0;
  TString argument="";
  int bMissingParam=0;
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

    // -datatype
    if (argument.CompareTo("-datatype")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      memcpy(&fDataType.fID, argv[i], TMath::Min(kAliHLTComponentDataTypefIDsize, (Int_t)strlen(argv[i])));
      if ((bMissingParam=(++i>=argc))) break;
      memcpy(&fDataType.fOrigin, argv[i], TMath::Min(kAliHLTComponentDataTypefOriginSize, (Int_t)strlen(argv[i])));

      // -dataspec
    } else if (argument.CompareTo("-dataspec")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      TString parameter(argv[i]);
      parameter.Remove(TString::kLeading, ' '); // remove all blanks
      if (parameter.IsDigit()) {
	fSpecification=(AliHLTUInt32_t)parameter.Atoi();
      } else if (parameter.BeginsWith("0x") &&
		 parameter.Replace(0,2,"",0).IsHex()) {
	sscanf(parameter.Data(),"%x", &fSpecification);
      } else {
	HLTError("wrong parameter for argument %s, number expected", argument.Data());
	iResult=-EINVAL;
      }
    } else {
      HLTError("unknown argument %s", argument.Data());
      break;
    }
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTRootFileStreamerComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/,
					    AliHLTComponentTriggerData& /*trigData*/ )
{
  // see header file for class documentation
  int iResult=0;
  AliHLTMemoryFile* pFile=CreateMemoryFile(fDataType,fSpecification);
  if (pFile) {
    const TObject* pObj=GetFirstInputObject(kAliHLTAnyDataType);
    int count=0;
    while (pObj && iResult>=0) {
      iResult=Write(pFile, pObj);
      if (iResult) {
	count++;
	HLTDebug("wrote object of class %s, data type %s", pObj->ClassName(), (DataType2Text(GetDataType(pObj)).c_str())); 
      }
      pObj=GetNextInputObject();
    }
    HLTInfo("wrote %d object(s) from %d input blocks to file", count, GetNumberOfInputBlocks());
    iResult=CloseMemoryFile(pFile);
  } else {
    iResult=-ENOMEM;
  }
  return iResult;
}
