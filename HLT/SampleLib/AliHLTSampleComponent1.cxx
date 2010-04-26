// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
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

//  @file   AliHLTSampleComponent1.cxx
//  @author Matthias Richter, Timm M. Steinbeck
//  @date   
//  @brief  A sample processing component for the HLT.
//          Component illustrates the basic functionality and component
//          initialization. 

#if __GNUC__== 3
using namespace std;
#endif

#include "AliHLTSampleComponent1.h"
#include "TString.h"
#include "TObjString.h"
#include "TMap.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTSampleComponent1)

/** one global instance used for registration */
AliHLTSampleComponent1 gAliHLTSampleComponent1;

AliHLTSampleComponent1::AliHLTSampleComponent1()
  : AliHLTProcessor()
  , fArgument1(0)
  , fArgument2(0)
{
  // an example component which implements the ALICE HLT processor
  // interface and illustrates the basic interface methods
  //
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  //
  // NOTE: all helper classes should be instantiated in DoInit()
}

AliHLTSampleComponent1::~AliHLTSampleComponent1()
{
  // destructor
  //
  // NOTE: implement proper cleanup in DoDeinit()
}

const char* AliHLTSampleComponent1::GetComponentID()
{ 
  // component property: id
  return "Sample-component1";
}

void AliHLTSampleComponent1::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // component property: list of input data types
    list.push_back(kAliHLTAnyDataType);
}

AliHLTComponentDataType AliHLTSampleComponent1::GetOutputDataType()
{
  // component property: output data type
  return kAliHLTVoidDataType;
}

void AliHLTSampleComponent1::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // component property: output size estimator
  constBase = 0;
  inputMultiplier = 0;
}

void AliHLTSampleComponent1::GetOCDBObjectDescription( TMap* const targetMap)
{
  // Get a list of OCDB object description.
  // The list of objects is provided in a TMap
  // - key: complete OCDB path, e.g. GRP/GRP/Data
  // - value: short description why the object is needed
  // Key and value objects created inside this class go into ownership of
  // target TMap.
  if (!targetMap) return;
  targetMap->Add(new TObjString("HLT/ConfigSample/SampleComponent1"),
		 new TObjString("configuration object"));
}

AliHLTComponent* AliHLTSampleComponent1::Spawn()
{
  // Spawn function, return new class instance
  return new AliHLTSampleComponent1;
}

int AliHLTSampleComponent1::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;

  // init stage 1: default values for all data members
  fArgument1=0;
  fArgument2=0;

  // init stage 2: read configuration object
  // ScanConfigurationArgument() needs to be implemented
  TString cdbPath="HLT/ConfigSample/SampleComponent1";
  iResult=ConfigureFromCDBTObjString(cdbPath);

  // init stage 3: read the component arguments
  if (iResult>=0) {
    iResult=ConfigureFromArgumentString(argc, argv);
  }

  if (iResult>=0) {
    // implement the component initialization
    if (!fArgument1) {
      HLTError("mandatory argument \'-mandatory1\' missing");
      iResult=-EPROTO;
    }
    if (!fArgument2) {
      HLTError("mandatory argument \'-mandatory2\' missing");
      iResult=-EPROTO;
    }
  }

  if (iResult<0) {
    // implement cleanup
  }

  return iResult;
}

int AliHLTSampleComponent1::ScanConfigurationArgument( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;

  TString argument="";
  TString configuration=""; 
  int bMissingParam=0;
  int i=0;
    argument=argv[i];
    if (argument.IsNull()) return 0;

    // -mandatory1
    if (argument.CompareTo("-mandatory1")==0) {
      if (++i>=argc) return -EPROTO;
      HLTInfo("got \'-mandatory1\' argument: %s", argv[i]);
      fArgument1=1;

      // -mandatory2
    } else if (argument.CompareTo("-mandatory2")==0) {
      fArgument2=1;
      HLTInfo("got \'-mandatory2\' argument");

      // -config1
    } else if (argument.CompareTo("-config1")==0) {
      if (++i>=argc) return -EPROTO;
      HLTInfo("got \'%s\' argument: %s", argument.Data(), argv[i]);

      // -config2
    } else if (argument.CompareTo("-config2")==0) {
      HLTInfo("got \'%s\' argument", argument.Data());

    } else {
      // no recognized argument
      i--;
    }

  return i+1;
}

int AliHLTSampleComponent1::DoDeinit()
{
  // see header file for class documentation
  HLTInfo("processing cleanup");
  return 0;
}

int AliHLTSampleComponent1::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
				      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
				      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks ) {
  // see header file for class documentation
  HLTInfo("processing data");
  if (evtData.fStructSize==0 && blocks==NULL && trigData.fStructSize==0 &&
      outputPtr==0 && size==0)
  {
    outputBlocks.clear();
    // this is just to get rid of the warning "unused parameter"
  }
  return 0;
}
int AliHLTSampleComponent1::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // see header file for class documentation
  int iResult=0;
  const char* path="HLT/ConfigSample/SampleComponent1";
  const char* defaultNotify="";
  if (cdbEntry) {
    path=cdbEntry;
    defaultNotify=" (default)";
  }

  HLTInfo("reconfigure from entry %s%s, chain id %s", path, defaultNotify,(chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
  iResult=ConfigureFromCDBTObjString(path);

  return iResult;
}

int AliHLTSampleComponent1::ReadPreprocessorValues(const char* modules)
{
  // see header file for class documentation
  int iResult=0;
  TString detectors(modules!=NULL?modules:"");
  HLTInfo("read preprocessor values for detector(s): %s", detectors.IsNull()?"none":detectors.Data());
  return iResult;
}
