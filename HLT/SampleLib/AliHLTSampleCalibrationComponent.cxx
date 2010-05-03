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

//  @file   AliHLTSampleCalibrationComponent.cxx
//  @author Matthias Richter
//  @date   2010-04-26
//  @brief  A sample calibration component for the HLT.
//  @ingroup alihlt_tutorial

#if __GNUC__== 3
using namespace std;
#endif

#include "AliHLTSampleCalibrationComponent.h"
#include "AliHLTReadoutList.h"
#include "AliLog.h"
#include "TMap.h"
#include "TObjString.h"
#include "TH1S.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTSampleCalibrationComponent)

/** one global instance used for registration */
AliHLTSampleCalibrationComponent gAliHLTSampleCalibrationComponent;

AliHLTSampleCalibrationComponent::AliHLTSampleCalibrationComponent()
  : AliHLTCalibrationProcessor()
  , fOutputSize(1000)
  , fHisto(NULL)
  , fHistoRange(10000)
{
  // an example component which implements the ALICE HLT calibration
  // processor interface
  //
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  //
  // NOTE: all helper classes should be instantiated in DoInit()
}

AliHLTSampleCalibrationComponent::~AliHLTSampleCalibrationComponent()
{
  // destructor
  //
  // NOTE: implement proper cleanup in DoDeinit()
}

const char* AliHLTSampleCalibrationComponent::GetComponentID()
{ 
  // component property: id
  return "SampleCalibration";
}

void AliHLTSampleCalibrationComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // component property: list of input data types
    list.push_back(kAliHLTAnyDataType);
}

AliHLTComponentDataType AliHLTSampleCalibrationComponent::GetOutputDataType()
{
  // component property: output data type
  return kAliHLTDataTypeFXSCalib|kAliHLTDataOriginSample;
}

void AliHLTSampleCalibrationComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // component property: output size estimator
  constBase = fOutputSize;
  inputMultiplier = 0;
}

void AliHLTSampleCalibrationComponent::GetOCDBObjectDescription( TMap* const targetMap)
{
  // Get a list of OCDB object description.
  // The list of objects is provided in a TMap
  // - key: complete OCDB path, e.g. GRP/GRP/Data
  // - value: short description why the object is needed
  // Key and value objects created inside this class go into ownership of
  // target TMap.
  if (!targetMap) return;
  targetMap->Add(new TObjString("HLT/ConfigSample/SampleCalibration"),
		 new TObjString("configuration object"));
}

AliHLTComponent* AliHLTSampleCalibrationComponent::Spawn()
{
  // Spawn function, return new class instance
  return new AliHLTSampleCalibrationComponent;
}

int AliHLTSampleCalibrationComponent::DoInit( int argc, const char** argv )
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
    fHisto=new TH1S("InputSize", "Input block size", 100, 0, fHistoRange);
    fOutputSize+=EstimateObjectSize(fHisto);
  }

  if (iResult<0) {
    // implement cleanup
  }

  return iResult;
}

int AliHLTSampleCalibrationComponent::ScanConfigurationArgument(int argc, const char** argv)
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
    if (++i>=argc) return -EPROTO;
    HLTInfo("got \'-mandatory1\' argument: %s", argv[i]);
    return 2; // keyword + 1 argument
  }

  // -optional1 arg
  if (argument.CompareTo("-optional1")==0) {
    if (++i>=argc) return -EPROTO;
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

int AliHLTSampleCalibrationComponent::DoDeinit()
{
  // component cleanup, delete all instances of helper classes here
  if (fHisto) delete fHisto;
  fOutputSize=0;

  return 0;
}

int AliHLTSampleCalibrationComponent::ProcessCalibration(const AliHLTComponentEventData& /*evtData*/,
							 AliHLTComponentTriggerData& /*trigData*/)
{
  // event processing function

  // check if this is a data event, there are a couple of special events
  // which should be ignored for normal processing
  if (!IsDataEvent()) return 0;

  // loop over input data blocks
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeDDLRaw);
       pBlock!=NULL; pBlock=GetNextInputBlock()) {
    HLTInfo("block %s specification 0x%x size %d", DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, pBlock->fSize);
    fHisto->Fill(pBlock->fSize);
  }

  // write the histogram out
  // this should not be done for every event, however the call can be implemented
  // like that and the publishing steered by the component argument 
  // '-pushback-period=...'
  if (PushBack(fHisto, kAliHLTDataTypeHistogram)==-ENOSPC) {
    // increase the output size estimator
    // we add the size of the last object, there might be other blocks to
    // be written in addition to the actual object
    fOutputSize+=GetLastObjectSize();
  }

  return 0;
}

int AliHLTSampleCalibrationComponent::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/, 
						     AliHLTComponentTriggerData& /*trigData*/)
{
  // prepare final result and ship to FXS

  AliHLTReadoutList rdList(AliHLTReadoutList::kHLT);
  PushToFXS(fHisto, "HLT", "TestHisto", rdList.Buffer());

  return 0;
}

int AliHLTSampleCalibrationComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // reconfigure the component from the specified CDB entry, or default CDB entry
  HLTInfo("reconfigure '%s' from entry %s", chainId, cdbEntry);

  return 0;
}

int AliHLTSampleCalibrationComponent::ReadPreprocessorValues(const char* modules)
{
  // read the preprocessor values for the detectors in the modules list
  int iResult=0;
  TString detectors(modules!=NULL?modules:"");
  HLTInfo("read preprocessor values for detector(s): %s", detectors.IsNull()?"none":detectors.Data());
  return iResult;
}
