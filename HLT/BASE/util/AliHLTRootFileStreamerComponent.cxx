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
//**************************************************************************

/// @file   AliHLTRootFileStreamerComponent.cxx
/// @author Matthias Richter
/// @date   
/// @brief  Save objects in a ROOT memory file
///

#include "AliHLTRootFileStreamerComponent.h"
#include "TString.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTRootFileStreamerComponent)

AliHLTRootFileStreamerComponent::AliHLTRootFileStreamerComponent()
  :
  AliHLTProcessor(),
  fDataType(kAliHLTVoidDataType),
  fSpecification(~(AliHLTUInt32_t)0)
{
  // The RootFileStreamer provides a stand alone component to write incoming
  // TObject like structures into a ROOT memory file. A ROOT memory file is
  // a ROOT file stored in memory instead on disk (AliHLTMemoryFile) The file
  // is published via the output stream. On the receiver side the file can
  // be directly written to disk and appears like a normal root file.
  //
  // Component ID: \b ROOTFileStreamer                                    <br>
  // Library: \b libAliHLTUtil.so                                         <br>
  // Input Data Types: ::kAliHLTAnyDataType                               <br>
  // Output Data Types: according to component arguments,
  //                    ::kAliHLTVoidDataType by default                  <br>
}

AliHLTRootFileStreamerComponent::~AliHLTRootFileStreamerComponent()
{
  // destructor
}

void AliHLTRootFileStreamerComponent::GetInputDataTypes( AliHLTComponentDataTypeList& list)
{
  // overloaded from AliHLTComponent
  list.clear();
  list.push_back(kAliHLTAllDataTypes);
}

AliHLTComponentDataType AliHLTRootFileStreamerComponent::GetOutputDataType()
{
  // overloaded from AliHLTComponent
  return fDataType;
}

void AliHLTRootFileStreamerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // overloaded from AliHLTComponent
  constBase=500;
  inputMultiplier=5.0;
}

int AliHLTRootFileStreamerComponent::DoInit( int argc, const char** argv )
{
  // overloaded from AliHLTComponent: initialization

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
  // overloaded from AliHLTProcessor: event processing
  int iResult=0;
  AliHLTMemoryFile* pFile=CreateMemoryFile(fDataType,fSpecification);
  if (pFile) {
    int count=0;
    for (const TObject* pObj=GetFirstInputObject();
	 pObj && iResult>=0;
	 pObj=GetNextInputObject()) {
      iResult=Write(pFile, pObj);
      if (iResult) {
	count++;
	HLTDebug("wrote object of class %s, data type %s", pObj->ClassName(), (DataType2Text(GetDataType(pObj)).c_str())); 
      }
    }
    HLTInfo("wrote %d object(s) from %d input blocks to file", count, GetNumberOfInputBlocks());
    iResult=CloseMemoryFile(pFile);
  } else {
    iResult=-ENOMEM;
  }
  return iResult;
}
