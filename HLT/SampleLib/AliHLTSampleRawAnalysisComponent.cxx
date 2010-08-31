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

/// @file   AliHLTSampleRawAnalysisComponent.cxx
/// @author Matthias Richter
/// @date   2010-08-29
/// @brief  A sample processing component for raw data
/// @ingroup alihlt_tutorial

#include "AliHLTSampleRawAnalysisComponent.h"
#include "AliLog.h"
#include "AliRawReaderMemory.h"
#include "TMap.h"
#include "TObjString.h"

/////////////////////////////////////////////////////////////////////////
// includes files likely to be skipped in a real implementation

// AliHLTDAQ needed to get generic mapping of data origin to DDL numbers
// not since a real component will select only blocks of the origin it
// is made for
#include "AliHLTDAQ.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTSampleRawAnalysisComponent)

/** one global instance used for  registration of
 * AliHLTSampleRawAnalysisComponent
 * Note: there are two ways of component registration
 * - via a global object
 * - via the AliHLTModuleAgent::RegisterComponents function
 * @see @ref alihlt_component_registration
 */
AliHLTSampleRawAnalysisComponent gAliHLTSampleRawAnalysisComponent;

AliHLTSampleRawAnalysisComponent::AliHLTSampleRawAnalysisComponent()
  : AliHLTProcessor()
  , fVerbosity(0)
  , fRawReader(NULL)
{
  // an example component which implements the ALICE HLT processor
  // interface and does some analysis on the input raw data
  //
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  //
  // NOTE: all helper classes should be instantiated in DoInit()
}

AliHLTSampleRawAnalysisComponent::~AliHLTSampleRawAnalysisComponent()
{
  // destructor
  //
  // NOTE: implement proper cleanup in DoDeinit()
}

const char* AliHLTSampleRawAnalysisComponent::GetComponentID()
{ 
  // component property: id
  return "SampleRawAnalysis";
}

void AliHLTSampleRawAnalysisComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // component property: list of input data types
  list.push_back(kAliHLTDataTypeDDLRaw);
}

AliHLTComponentDataType AliHLTSampleRawAnalysisComponent::GetOutputDataType()
{
  // component property: output data type
  return kAliHLTDataTypeHistogram|kAliHLTDataOriginSample;
}

void AliHLTSampleRawAnalysisComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // component property: output size estimator
  constBase = 1000;
  inputMultiplier = 0;
}

void AliHLTSampleRawAnalysisComponent::GetOCDBObjectDescription( TMap* const targetMap)
{
  // Get a list of OCDB object description.
  // The list of objects is provided in a TMap
  // - key: complete OCDB path, e.g. GRP/GRP/Data
  // - value: short description why the object is needed
  // Key and value objects created inside this class go into ownership of
  // target TMap.
  if (!targetMap) return;
  targetMap->Add(new TObjString("HLT/ConfigSample/SampleRawAnalysis"),
		new TObjString("configuration object"));
}

AliHLTComponent* AliHLTSampleRawAnalysisComponent::Spawn()
{
  // Spawn function, return new class instance
  return new AliHLTSampleRawAnalysisComponent;
}

int AliHLTSampleRawAnalysisComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;

  // init stage 1: default values for all data members

  // init stage 2: read configuration object
  // ScanConfigurationArgument() needs to be implemented
  // HLT component configuration objects are located in the HLT/ConfigDET
  // of the OCDB. TObjString configuration objects can be generated with
  // the macro HLT/exa/makeComponentConfigurationObject.C, e.g.
  // aliroot -b -q -l $ALICE_ROOT/HLT/exa/makeComponentConfigurationObject.C'("HLT/ConfigSample/SampleRawAnalysis", "")'
  TString cdbPath="HLT/ConfigSample/";
  cdbPath+=GetComponentID();
  iResult=ConfigureFromCDBTObjString(cdbPath);

  // init stage 3: read the component arguments
  if (iResult>=0) {
    iResult=ConfigureFromArgumentString(argc, argv);
  }

  // implement the component initialization
  do {
    if (iResult<0) break;

    fRawReader=new AliRawReaderMemory;
    if (!fRawReader) {
      iResult=-ENOMEM;
      break;
    }

    // implement further initialization
  } while (0);

  if (iResult<0) {
    // implement cleanup
    if (fRawReader) delete fRawReader;
    fRawReader=NULL;
  }

  return iResult;
}

int AliHLTSampleRawAnalysisComponent::ScanConfigurationArgument(int argc, const char** argv)
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

  // -verbose
  if (argument.CompareTo("-verbose")==0) {
    fVerbosity=1;
    return 1; // only keyword
  }

  return 0;
}

int AliHLTSampleRawAnalysisComponent::DoDeinit()
{
  // component cleanup, delete all instances of helper classes here
  if (fRawReader) delete fRawReader;
  fRawReader=NULL;

  return 0;
}

int AliHLTSampleRawAnalysisComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/,
					      AliHLTComponentTriggerData& /*trigData*/)
{
  // event processing function
  int iResult=0;

  // check if this is a data event, there are a couple of special events
  // which should be ignored for normal processing
  if (!IsDataEvent()) return 0;

  if (fVerbosity>0) AliInfoClass(Form("==================== event %3d ================================", GetEventCount()));

  // loop over the raw input data blocks and set up the rawreader
  for (const AliHLTComponentBlockData* pBlock = GetFirstInputBlock(kAliHLTDataTypeDDLRaw);
       pBlock!=NULL && iResult>=0;
       pBlock=GetNextInputBlock()) {
    // extract DDL id from specification
    int ddlnum=-1;
    for (unsigned pos=0; pos<8*sizeof(AliHLTUInt32_t); pos++) {
      if (pBlock->fSpecification & (0x1<<pos)) {
	if (ddlnum>=0) {
	  // this is just an example, please avoid warnings in every event since those will
	  // seturate the logging system. Consider AliHLTErrorGuard for such cases, e.g.
	  // ALIHLTERRORGUARD(5, "nasty error, first occurence in event %d", event);
	  // this will show the error 5 times and then a summary at the end
	  HLTWarning("Can not uniquely identify DDL number from specification, skipping data block %s 0x%08x",
		     DataType2Text(pBlock->fDataType).c_str(),
		     pBlock->fSpecification);
	  ddlnum=-1;
	  break;
	}
	ddlnum=pos;
      }
    }
    if (ddlnum<0) continue;
    int id=AliHLTDAQ::DdlID(AliHLTDAQ::DetectorName(pBlock->fDataType.fOrigin), ddlnum);
    if (id<0) {
      // this is just an example, please avoid warnings in every event since those will
      // seturate the logging system. Consider AliHLTErrorGuard for such cases, e.g.
      // ALIHLTERRORGUARD(5, "nasty error, first occurence in event %d", event);
      // this will show the error 5 times and then a summary at the end
      HLTWarning("failed to get DDL id of link %d for data block %s, 0x%08x", ddlnum,
		 DataType2Text(pBlock->fDataType).c_str(),
		 pBlock->fSpecification);
      continue;
    }

    // add data block to rawreader
    if(!fRawReader->AddBuffer((UChar_t*) pBlock->fPtr, pBlock->fSize, id)){
      HLTError("Could not add buffer of data block  %s, 0x%08x to rawreader",
	       DataType2Text(pBlock->fDataType).c_str(),
	       pBlock->fSpecification);
    }
  }

  // now use raw reader to do some data processing
  UChar_t* pData=NULL;
  while (fRawReader->ReadNextData(pData)) {
    AliInfoClass(Form("got data for equipment id %d: size %d", fRawReader->GetEquipmentId(), fRawReader->GetEquipmentSize()));
  }

  // clear the rawreader
  fRawReader->ClearBuffers();    

  // publish the array of tracks as output

  return iResult;
}

int AliHLTSampleRawAnalysisComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // reconfigure the component from the specified CDB entry, or default CDB entry
  // function is invoked by the framework if a reconfigure command was received.
  // 
  int iResult=0;
  TString cdbPath;
  if (cdbEntry) {
    cdbPath=cdbEntry;
  } else {
    cdbPath="HLT/ConfigSample/";
    cdbPath+=GetComponentID();
  }
  AliInfoClass(Form("reconfigure '%s' from entry %s%s", chainId, cdbPath.Data(), cdbEntry?"":" (default)"));
  iResult=ConfigureFromCDBTObjString(cdbPath);

  return iResult;
}

int AliHLTSampleRawAnalysisComponent::ReadPreprocessorValues(const char* modules)
{
  // read the preprocessor values for the detectors in the modules list
  // function is invoked by the framework if the pendolino indivates an update
  // of online calibration objects, e.g temperature and pressure measurements.
  int iResult=0;
  TString detectors(modules!=NULL?modules:"");
  AliInfoClass(Form("read preprocessor values for detector(s): %s", detectors.IsNull()?"none":detectors.Data()));
  return iResult;
}
