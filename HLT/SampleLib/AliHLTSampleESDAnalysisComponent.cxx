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

//  @file   AliHLTSampleESDAnalysisComponent.cxx
//  @author Matthias Richter
//  @date   2010-04-17
//  @brief  A sample processing component for ESD analysis.
//  @ingroup alihlt_tutorial

#if __GNUC__== 3
using namespace std;
#endif

#include "AliHLTSampleESDAnalysisComponent.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "TMap.h"
#include "TObjString.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTSampleESDAnalysisComponent)

/** one global instance used for registration */
AliHLTSampleESDAnalysisComponent gAliHLTSampleESDAnalysisComponent;

AliHLTSampleESDAnalysisComponent::AliHLTSampleESDAnalysisComponent()
{
  // an example component which implements the ALICE HLT processor
  // interface and does some analysis on the input ESD
  //
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  //
  // NOTE: all helper classes should be instantiated in DoInit()
}

AliHLTSampleESDAnalysisComponent::~AliHLTSampleESDAnalysisComponent()
{
  // destructor
  //
  // NOTE: implement proper cleanup in DoDeinit()
}

const char* AliHLTSampleESDAnalysisComponent::GetComponentID()
{ 
  // component property: id
  return "SampleESDAnalysis";
}

void AliHLTSampleESDAnalysisComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // component property: list of input data types
    list.push_back(kAliHLTDataTypeESDObject);
}

AliHLTComponentDataType AliHLTSampleESDAnalysisComponent::GetOutputDataType()
{
  // component property: output data type
  return kAliHLTDataTypeTObjArray|kAliHLTDataOriginSample;
}

void AliHLTSampleESDAnalysisComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // component property: output size estimator
  constBase = 0;
  inputMultiplier = 0;
}

void AliHLTSampleESDAnalysisComponent::GetOCDBObjectDescription( TMap* const targetMap)
{
  // Get a list of OCDB object description.
  // The list of objects is provided in a TMap
  // - key: complete OCDB path, e.g. GRP/GRP/Data
  // - value: short description why the object is needed
  // Key and value objects created inside this class go into ownership of
  // target TMap.
  if (!targetMap) return;
  targetMap->Add(new TObjString("HLT/ConfigSample/SampleESDAnalysis"),
		new TObjString("configuration object"));
  targetMap->Add(new TObjString("GRP/GRP/Data"),
		new TObjString("GRP object"));
}

AliHLTComponent* AliHLTSampleESDAnalysisComponent::Spawn()
{
  // Spawn function, return new class instance
  return new AliHLTSampleESDAnalysisComponent;
}

int AliHLTSampleESDAnalysisComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;

  // init stage 1: default values for all data members

  // init stage 2: read configuration object
  // ScanConfigurationArgument() needs to be implemented
  TString cdbPath="HLT/ConfigSample/";
  cdbPath+=GetComponentID();
  iResult=ConfigureFromCDBTObjString(cdbPath);

  // init stage 3: read the component arguments
  if (iResult>=0) {
    iResult=ConfigureFromArgumentString(argc, argv);
  }

  if (iResult>=0) {
    // implement the component initialization
  }

  if (iResult<0) {
    // implement cleanup
  }

  return iResult;
}

int AliHLTSampleESDAnalysisComponent::ScanConfigurationArgument(int argc, const char** argv)
{
  // Scan configuration arguments
  // Return the number of processed arguments
  //        -EPROTO if argument format error (e.g. number expected but not found)
  //
  // The AliHLTComponent base class implements a parsing loop for argument strings and
  // arrays of strings which is invoked by ConfigureFromArgumentString/ConfigureFromCDBTObjString
  // The component needs to implement ScanConfigurationArgument in order to decode the arguments.

  int i=0;
  TString argument=argv[i];

  if (argument.IsNull()) return 0;

  // -mandatory1 arg
  if (argument.CompareTo("-mandatory1")==0) {
    if (++i>=argc) return -EINVAL;
    HLTInfo("got \'-mandatory1\' argument: %s", argv[i]);
    return 2; // keyword + 1 argument
  }

  // -optional1 arg
  if (argument.CompareTo("-optional1")==0) {
    if (++i>=argc) return -EINVAL;
    HLTInfo("got \'-optional1\' argument: %s", argv[i]);
    return 2; // keyword + 1 argument
  }

  // -optional2
  if (argument.CompareTo("-optional2")==0) {
    HLTInfo("got \'-optional2\' argument");
    return 1; // only keyword
  }

  return 0;
}

int AliHLTSampleESDAnalysisComponent::DoDeinit()
{
  // component cleanup, delete all instances of helper classes here

  return 0;
}

int AliHLTSampleESDAnalysisComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/,
					      AliHLTComponentTriggerData& /*trigData*/)
{
  // event processing function

  // check if this is a data event, there are a couple of special events
  // which should be ignored for normal processing
  if (!IsDataEvent()) return 0;

  const TObject* obj = GetFirstInputObject(kAliHLTAllDataTypes, "AliESDEvent");

  // input objects are not supposed to be changed by the component, so they
  // are defined const. However, the implementation of AliESDEvent does not
  // support this and we need the const_cast
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(const_cast<TObject*>(obj));
  TObjArray output;
  if (esd != NULL) {
    AliInfoClass(Form("==================== event %3d ================================", GetEventCount()));
    esd->GetStdContent();
    for (Int_t i = 0; i < esd->GetNumberOfTracks(); i++) {
      AliESDtrack* track = esd->GetTrack(i);
      AliInfoClass(Form("-------------------- track %3d --------------------------------", i));
      track->Print("");
      output.Add(track);
    }
  }

  // publish the array of tracks as output
  PushBack(&output, kAliHLTDataTypeTObjArray|kAliHLTDataOriginSample);

  return 0;
}

int AliHLTSampleESDAnalysisComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // reconfigure the component from the specified CDB entry, or default CDB entry
  HLTInfo("reconfigure '%s' from entry %s", chainId, cdbEntry);

  return 0;
}

int AliHLTSampleESDAnalysisComponent::ReadPreprocessorValues(const char* modules)
{
  // read the preprocessor values for the detectors in the modules list
  int iResult=0;
  TString detectors(modules!=NULL?modules:"");
  HLTInfo("read preprocessor values for detector(s): %s", detectors.IsNull()?"none":detectors.Data());
  return iResult;
}
